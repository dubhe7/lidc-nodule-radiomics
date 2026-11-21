# -*- coding: utf-8 -*-
"""
LIDC-IDRI 子集下载器  v0.2
=========================================================
为什么不用 NBIA Data Retriever：
  官方 manifest 是全量 1018 例 / 133 GB，GUI 里逐个取消勾选不现实。
  本脚本走 TCIA 的 NBIA Search REST API，按病例数截断、只取 CT(+SEG)，
  默认 100 例 ≈ 8~10 GB。

端点依据（TCIA 官方 NBIA Search REST API Guide）：
  base  = https://services.cancerimagingarchive.net/nbia-api/services/v1
  getSeries?Collection=LIDC-IDRI&format=json   -> 全量 series 元数据（一次拉完）
  getImage?SeriesInstanceUID=<uid>             -> 该 series 的 ZIP

用法（在你本机的普通终端里跑；沙箱/受限网络里跑不通）：

    # 先测 1 例，确认能通再批量
    python download_lidc.py --ids LIDC-IDRI-0002

    # 下前 100 例
    python download_lidc.py --limit 100

    # 只要 CT 不要 SEG
    python download_lidc.py --limit 100 --modalities CT

输出结构（与现有 LIDC-IDRI-0001 完全一致，可直接被 extract_case.py 读）：
    <out>/LIDC-IDRI-0001/<SeriesUID>/CT_<SeriesUID>/*.dcm
    <out>/LIDC-IDRI-0001/<SeriesUID>/SEG_<SeriesUID>/*.dcm

断点续传：重复执行会跳过已有文件的 series，中断后直接重跑即可。
只依赖标准库。
"""
import os
import io
import time
import json
import argparse
import zipfile
import urllib.request
import urllib.parse

BASE = "https://services.cancerimagingarchive.net/nbia-api/services/v1"
COLLECTION = "LIDC-IDRI"
UA = {"User-Agent": "Mozilla/5.0 (academic research; lidc-subset-downloader/0.3)"}

# 运行期由 detect_channel() 设置；_get 一律经由这个 opener 发请求
OPENER = None


def detect_channel(force_proxy=""):
    """探测可用的网络通道：直连 or 本机代理。

    TCIA 对国内直连典型表现是 TLS 握手阶段被 RST（WinError 10054），
    所以探测必须真的发一次 HTTPS 请求，TCP 连上不算数。
    """
    global OPENER
    test_url = "https://services.cancerimagingarchive.net/"

    # 先读 Windows 系统代理设置（即使没启用，历史端口也有参考价值）
    reg_ports = []
    try:
        import winreg
        k = winreg.OpenKey(winreg.HKEY_CURRENT_USER,
                           r"Software\Microsoft\Windows\CurrentVersion\Internet Settings")
        server, _ = winreg.QueryValueEx(k, "ProxyServer")
        for part in server.replace(";", ",").split(","):
            part = part.strip()
            # 形如 127.0.0.1:10810 或 https=host:port
            if ":" in part:
                hp = part.split("=")[-1]
                if hp.count(":") == 1:
                    reg_ports.append(hp)
    except Exception:
        pass

    if force_proxy:
        candidates = [("代理 " + force_proxy,
                       urllib.request.ProxyHandler(
                           {"http": f"http://{force_proxy}", "https": f"http://{force_proxy}"}))]
    else:
        candidates = [("直连", None)]
        ports = ["7890", "7897", "10809", "10810", "10808", "8118"]
        ports += [p for p in reg_ports if p not in ports]
        for port in ports:
            candidates.append((
                f"本机代理 127.0.0.1:{port}",
                urllib.request.ProxyHandler(
                    {"http": f"http://127.0.0.1:{port}", "https": f"http://127.0.0.1:{port}"}),
            ))

    for name, handler in candidates:
        opener = urllib.request.build_opener(handler) if handler else urllib.request.build_opener()
        try:
            req = urllib.request.Request(test_url, headers=UA)
            with opener.open(req, timeout=15) as r:
                r.read(256)
            print(f"[0] 网络通道 OK：{name}")
            OPENER = opener
            return
        except Exception as e:
            print(f"    通道不可用（{name}）：{type(e).__name__} {str(e)[:80]}")
    raise RuntimeError(
        "所有网络通道都探测失败。\n"
        "  1) 确认你的代理软件（Clash/V2Ray 等）正在运行，并查看它的 HTTP 端口；\n"
        "  2) 用 --proxy 127.0.0.1:端口 强制指定，例如：\n"
        "     python download_lidc.py --ids LIDC-IDRI-0002 --proxy 127.0.0.1:7890")


def _get(url, retries=4, timeout=90):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with OPENER.open(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:
            last = e
            wait = 2 ** i
            print(f"      重试 {i+1}/{retries} ({type(e).__name__}) {wait}s 后继续...")
            time.sleep(wait)
    raise RuntimeError(f"请求失败: {url}\n  {last}")


def fetch_all_series():
    """一次请求拿全量 series（1398 条），按 PatientID 分组"""
    print(f"[1] 拉取 {COLLECTION} 全量 series 列表 ...")
    url = f"{BASE}/getSeries?{urllib.parse.urlencode({'Collection': COLLECTION, 'format': 'json'})}"
    data = json.loads(_get(url, timeout=180).decode("utf-8", "replace"))
    print(f"    返回 {len(data)} 条 series 记录")

    by_pid = {}
    for s in data:
        pid = s.get("PatientID")
        if not pid:
            continue
        by_pid.setdefault(pid, []).append({
            "uid": s["SeriesInstanceUID"],
            "modality": str(s.get("Modality", "")).upper(),
            "desc": s.get("SeriesDescription", ""),
        })
    return dict(sorted(by_pid.items()))


def download_series(pid, ser, out_root, skip_existing=True):
    uid = ser["uid"]
    target_dir = os.path.join(out_root, pid, uid, f"{ser['modality']}_{uid}")

    if skip_existing and os.path.isdir(target_dir):
        if any(f.lower().endswith(".dcm") for f in os.listdir(target_dir)):
            print(f"      跳过（已存在） {ser['modality']:3s} ...{uid[-12:]}")
            return 0

    url = f"{BASE}/getImage?{urllib.parse.urlencode({'SeriesInstanceUID': uid})}"
    blob = _get(url, retries=3, timeout=300)

    if blob[:2] != b"PK":
        # 不是 zip：可能返回了错误页/空内容，报出来而不是静默写坏文件
        head = blob[:200].decode("utf-8", "replace").replace("\n", " ")
        raise RuntimeError(f"返回不是 ZIP（{len(blob)}B）：{head}")

    os.makedirs(target_dir, exist_ok=True)
    n = 0
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        for name in z.namelist():
            if name.endswith("/"):
                continue
            base = os.path.basename(name)
            if not base:
                continue
            with z.open(name) as src, open(os.path.join(target_dir, base), "wb") as dst:
                dst.write(src.read())
            n += 1

    print(f"      {ser['modality']:3s} ...{uid[-12:]}  {n:>3d} 文件  {len(blob)/1024/1024:>6.1f} MB")
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=100, help="下载前 N 例（默认 100）")
    ap.add_argument("--ids", type=str, default="", help="指定病例，逗号分隔，优先于 --limit")
    ap.add_argument("--out", type=str, default=r"E:/lidc_data/lidc_idri", help="输出根目录")
    ap.add_argument("--modalities", type=str, default="CT,SEG",
                    help="要下载的模态，逗号分隔（默认 CT,SEG；只要图像用 CT）")
    ap.add_argument("--proxy", type=str, default="",
                    help="强制指定代理，如 127.0.0.1:7890（默认自动探测）")
    args = ap.parse_args()

    detect_channel(args.proxy)

    want = {m.strip().upper() for m in args.modalities.split(",") if m.strip()}
    os.makedirs(args.out, exist_ok=True)

    groups = fetch_all_series()

    if args.ids:
        pids = [x.strip() for x in args.ids.split(",") if x.strip()]
        print(f"[2] 使用指定病例 {len(pids)} 例")
    else:
        pids = list(groups.keys())[:args.limit]
        print(f"[2] 取前 {len(pids)} 例（共 {len(groups)} 例可选）")

    # 预统计
    planned = 0
    for pid in pids:
        planned += sum(1 for s in groups.get(pid, []) if s["modality"] in want)
    print(f"    计划下载 {planned} 个 series\n")

    t0 = time.time()
    total_files, done_series, failed = 0, 0, []

    for i, pid in enumerate(pids, 1):
        print(f"[{i}/{len(pids)}] {pid}")
        series = [s for s in groups.get(pid, []) if s["modality"] in want]
        if not series:
            print("      无目标模态，跳过")
            continue
        # CT 优先，保证中断时先拿到图像
        series.sort(key=lambda x: (0 if x["modality"] == "CT" else 1))

        for ser in series:
            try:
                total_files += download_series(pid, ser, args.out)
                done_series += 1
            except Exception as e:
                print(f"      !! 失败: {type(e).__name__}: {str(e)[:130]}")
                failed.append((pid, ser["uid"]))

    dt = (time.time() - t0) / 60
    print(f"\n{'='*60}")
    print(f"完成：{done_series}/{planned} 个 series，{total_files} 个文件，耗时 {dt:.1f} 分钟")
    print(f"数据目录：{os.path.abspath(args.out)}")
    if failed:
        print(f"失败 {len(failed)} 个 series（重跑本脚本会自动续传）：")
        for pid, uid in failed[:10]:
            print(f"  {pid}  ...{uid[-12:]}")
    print(f"\n下一步：python extract_case.py \"{os.path.abspath(args.out)}/LIDC-IDRI-0002\"")


if __name__ == "__main__":
    main()
