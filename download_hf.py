# -*- coding: utf-8 -*-
"""
LIDC-IDRI 国内下载器（HuggingFace 国内镜像 hf-mirror.com）  v1.0
=================================================================
为什么用这条路而不是 TCIA：
  TCIA 国内直连被 RST（WinError 10054），需全程挂代理；
  魔搭 ModelScope 有图像（5.8 MB/s，很快）但**只有图像、没有掩膜和标签**；
  hf-mirror.com 上 MedOtter/LIDC-IDRI 是完整镜像：
    images/   CT 图像（NIfTI，int16，HU，原始几何）
    masks/    共识掩膜（uint8 {0,1}，≥3 位医师标注的结节按 50% 共识合成）
    nodules.csv  每个结节的形态评分与 malignancy_score（1-5）
    aux/LIDC-XML-only.zip  原始 4 医师 XML（8.6 MB，可选）
  镜像无需翻墙，实测单连接 0.4 MB/s、多线程可聚合到 1.3+ MB/s。

输出结构（与 extract_case.py 的 DICOM 路线并存，互不干扰）：
    <out>/images/LIDC-IDRI-0001_9720603192.nii.gz
    <out>/masks/LIDC-IDRI-0001_9720603192.nii.gz
    <out>/labels.csv      每例一行：结节数、最大/平均恶性度、是否有共识掩膜
    <out>/nodules.csv     原始结节级标注（含全部形态评分）

用法：
    python download_hf.py --limit 3            # 先试 3 例（约 150 MB）
    python download_hf.py --limit 100          # 按 ID 顺序前 N 例（注意：有入组顺序偏倚！）
    python download_hf.py --limit 100 --workers 8
    python download_hf.py --sample 300         # 【推荐】按恶性度分层随机抽 300 例（全集比例）
    python download_hf.py --sample 300 --dry-run   # 只打印抽样结果，不下载
    python download_hf.py --backfill           # 补"图到了掩膜没到"的半成品，重写 labels.csv
中断后重跑即续传（已存在且大小一致的文件跳过）。
labels.csv 每次按"磁盘上已有的 + 本次选中的"全量重写，不会丢行。
只依赖标准库。
"""
import os
import csv
import sys
import gzip
import time
import random
import argparse
import threading
import urllib.request
import urllib.error

MIRROR = "https://hf-mirror.com"
REPO = "MedOtter/LIDC-IDRI"
RESOLVE = f"{MIRROR}/datasets/{REPO}/resolve/main/"
UA = {"User-Agent": "Mozilla/5.0 (academic research; lidc-cn-downloader/1.0)"}

_print_lock = threading.Lock()


def log(msg):
    with _print_lock:
        print(msg, flush=True)


def fetch(url, retries=4, timeout=180):
    """下载到内存，失败重试（指数退避）。"""
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:
            last = e
            if i + 1 < retries:
                time.sleep(2 ** i)
    raise RuntimeError(f"{type(last).__name__}: {str(last)[:120]}")


def fetch_to_file(url, path, retries=4):
    """流式下载到文件；先写 .part 再改名，避免中断留下半个文件被当成已完成。"""
    tmp = path + ".part"
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            got = 0
            with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "wb") as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
                    got += len(chunk)
            os.replace(tmp, path)
            return got
        except Exception as e:
            last = e
            if os.path.exists(tmp):
                try:
                    os.remove(tmp)
                except OSError:
                    pass
            if attempt + 1 < retries:
                time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"{type(last).__name__}: {str(last)[:120]}")


def load_nodules(out_dir):
    """下载并解析 nodules.csv（结节级标注）。"""
    csv_path = os.path.join(out_dir, "nodules.csv")
    if not os.path.exists(csv_path) or os.path.getsize(csv_path) < 1000:
        log("[1] 下载 nodules.csv ...")
        blob = fetch(RESOLVE + "nodules.csv")
        with open(csv_path, "wb") as f:
            f.write(blob)
    with open(csv_path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    log(f"    nodules.csv: {len(rows)} 个结节记录")
    return rows


def build_case_table(rows):
    """按 case_id 汇总：只保留至少有一个『进默认 GT（≥3 位医师）』结节的病例。
    这类病例才有干净的掩膜和可用的恶性度标签。"""
    cases = {}
    for r in rows:
        cid = r["case_id"]
        c = cases.setdefault(cid, {
            "case_id": cid,
            "patient_id": r["patient_id"],
            "n_nodules": 0,
            "n_gt": 0,
            "malignancies": [],
        })
        c["n_nodules"] += 1
        if r.get("in_default_gt", "").lower() == "true":
            c["n_gt"] += 1
            try:
                c["malignancies"].append(float(r["malignancy_score"]))
            except (KeyError, ValueError):
                pass
    usable = [c for c in cases.values() if c["n_gt"] > 0]
    usable.sort(key=lambda x: x["patient_id"])
    return usable


def gz_ok(path):
    """完整解压一遍验证 gzip 流未截断（HTTP 200 但 body 不全的情况只有这里能抓到）。"""
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1 << 22):
                pass
        return True
    except Exception:
        return False


def label_class(c):
    """按 case 级 max_malignancy 分三类（与论文常用口径一致）。"""
    ms = c["malignancies"]
    if not ms:
        return None
    m = max(ms)
    if m >= 4:
        return "mal"   # 恶性
    if m <= 2:
        return "ben"   # 良性
    return "mid"       # 灰区（=3）


def already_downloaded(c, out_dir):
    """图像与掩膜都在才算完成（只看图像会漏掉"图到掩膜没到"的半成品）。"""
    return (os.path.exists(os.path.join(out_dir, "images", c["case_id"] + ".nii.gz"))
            and os.path.exists(os.path.join(out_dir, "masks", c["case_id"] + ".nii.gz")))


def stable_sample(groups, n_total, seed):
    """每类用独立 RNG 把该类全集洗牌后取前 q。

    关键：随机序列只依赖 nodules.csv 里该类成员（稳定不变），
    不依赖磁盘上已有哪些文件 —— 因此跨次运行抽出的样本完全一致。
    （旧写法 rng.shuffle(rest) 会随磁盘状态改变序列，同一 seed 抽出不同批次。）
    """
    total_classified = sum(len(v) for v in groups.values())
    quotas = {k: round(n_total * len(v) / total_classified) for k, v in groups.items()}
    picked = {}
    for k, v in groups.items():
        pool = sorted(v, key=lambda c: c["case_id"])
        random.Random(f"{seed}:{k}").shuffle(pool)
        picked[k] = pool[:quotas[k]]
    return picked, quotas


def stratified_select(usable, out_dir, n_total, seed=42, verbose=True):
    """按全集类别比例分层随机抽样；已在磁盘上的病例无条件保留。

    返回 (selected 全部选中, to_download 需补下, quotas 各类配额)。
    抽样成员由 stable_sample 决定，与磁盘状态无关，可复现、可写进方法学。
    """
    groups = {"mal": [], "ben": [], "mid": []}
    for c in usable:
        k = label_class(c)
        if k:
            groups[k].append(c)
    picked, quotas = stable_sample(groups, n_total, seed)

    selected, to_download = [], []
    seen = set()
    for k in ("mal", "ben", "mid"):
        have_list = [c for c in groups[k] if already_downloaded(c, out_dir)]
        need = [c for c in picked[k] if c["case_id"] not in {h["case_id"] for h in have_list}]
        for c in have_list + need:
            if c["case_id"] not in seen:
                seen.add(c["case_id"])
                selected.append(c)
        to_download.extend(need)
        if verbose:
            log(f"    [{k}] 全集 {len(groups[k])} 例 | 已在盘 {len(have_list)} | "
                f"本类目标 {max(len(picked[k]), len(have_list))} | 需补下 {len(need)}")
    return selected, to_download, quotas


def download_case(c, out_dir, overwrite=False):
    """下图像 + 掩膜，返回 (ok, 解压后字节数)。"""
    files = [("images", c["case_id"] + ".nii.gz"), ("masks", c["case_id"] + ".nii.gz")]
    total = 0
    for sub, name in files:
        path = os.path.join(out_dir, sub, name)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if os.path.exists(path) and not overwrite and os.path.getsize(path) > 1024:
            total += os.path.getsize(path)
            continue
        total += fetch_to_file(f"{RESOLVE}{sub}/{name}", path)
    return total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=100, help="下载前 N 个可用病例（默认 100；按 ID 顺序，有偏倚）")
    ap.add_argument("--sample", type=int, default=0, help="按恶性度分层随机抽 N 例（推荐，覆盖 --limit）")
    ap.add_argument("--seed", type=int, default=42, help="随机抽样种子（默认 42，保证可复现）")
    ap.add_argument("--dry-run", action="store_true", help="只打印抽样结果与统计，不下载")
    ap.add_argument("--backfill", action="store_true", help="不抽样：只为已有图像、缺掩膜的病例补掩膜，并重写 labels.csv")
    ap.add_argument("--verify", action="store_true", help="解压校验全部文件，删除截断的并重下（建议每次跑完都来一遍）")
    ap.add_argument("--workers", type=int, default=6, help="并发线程数（默认 6；单连接被限速，并发才跑得动）")
    ap.add_argument("--out", type=str, default=r"E:/lidc_data/lidc_hf", help="输出目录")
    ap.add_argument("--overwrite", action="store_true", help="已存在的文件也重新下载")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    log(f"镜像: {MIRROR}   仓库: {REPO}")
    log(f"输出: {os.path.abspath(args.out)}")

    rows = load_nodules(args.out)
    usable = build_case_table(rows)
    log(f"[2] 有『≥3 位医师共识结节』的可用病例: {len(usable)} 例")

    # 完整性校验：删掉截断文件，稍后重下
    repair_ids = set()
    if args.verify:
        removed = 0
        for sub in ("images", "masks"):
            d = os.path.join(args.out, sub)
            if not os.path.isdir(d):
                continue
            for fn in os.listdir(d):
                if not fn.endswith(".nii.gz"):
                    continue
                p = os.path.join(d, fn)
                if not gz_ok(p):
                    repair_ids.add(fn[:-len(".nii.gz")])
                    os.remove(p)
                    removed += 1
                    log(f"  !! 损坏已删 {sub}/{fn} ({os.path.getsize(p) if os.path.exists(p) else 0})")
        log(f"[1] verify: 校验完成，损坏 {removed} 个（将重新下载）")
        if not args.sample:
            args.backfill = True

    if args.backfill:
        targets = [c for c in usable
                   if os.path.exists(os.path.join(args.out, "images", c["case_id"] + ".nii.gz"))
                   or c["case_id"] in repair_ids]
        to_download = [c for c in targets if not already_downloaded(c, args.out)]
        quotas = None
        log(f"[3] backfill: 目标 {len(targets)} 例，其中需补下 {len(to_download)} 例")
    elif args.sample > 0:
        targets, to_download, quotas = stratified_select(usable, args.out, args.sample, args.seed)
    else:
        targets = usable[:args.limit]
        to_download = [c for c in targets if args.overwrite or not already_downloaded(c, args.out)]
        quotas = None

    if not targets:
        log("没有可用病例，退出")
        return

    # 类别统计（labels.csv 与磁盘保持一致的关键：targets 已包含全部在盘病例）
    cnt = {"mal": 0, "ben": 0, "mid": 0}
    for c in targets:
        k = label_class(c)
        if k:
            cnt[k] += 1
    log(f"[3] 选中 {len(targets)} 例：恶性(>=4) {cnt['mal']} | 良性(<=2) {cnt['ben']} | "
        f"灰区(=3) {cnt['mid']} | 恶:良 = {cnt['mal']/max(cnt['ben'],1):.2f} : 1")
    log(f"    本次需新下载 {len(to_download)} 例，预估 {len(to_download)*45.0/1024:.1f} GB\n")
    if args.dry_run:
        log("[dry-run] 不下载，退出")
        return

    # 并发下载（只下载缺口；已存在的病例文件级跳过）
    t0 = time.time()
    done = [0]
    errors = []
    lock = threading.Lock()
    idx = [0]

    def worker():
        while True:
            with lock:
                if idx[0] >= len(to_download):
                    return
                c = to_download[idx[0]]
                idx[0] += 1
                n = idx[0]
            try:
                b = download_case(c, args.out, args.overwrite)
                with lock:
                    done[0] += 1
                log(f"  [{n}/{len(to_download)}] {c['patient_id']}  {b/1024/1024:>6.1f} MB  "
                    f"结节{c['n_gt']}个 恶性度{max(c['malignancies']) if c['malignancies'] else '-'}")
            except Exception as e:
                with lock:
                    errors.append((c["patient_id"], c["case_id"], str(e)))
                log(f"  [{n}/{len(to_download)}] {c['patient_id']}  !! 失败 {str(e)[:90]}")

    threads = [threading.Thread(target=worker, daemon=True) for _ in range(max(1, args.workers))]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    # 写标签表（放下载之后：失败的病例不写入，保证 labels.csv 与磁盘图像一一对应）
    failed_ids = {cid for _, cid, _ in errors}
    label_path = os.path.join(args.out, "labels.csv")
    with open(label_path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case_id", "patient_id", "n_nodules", "n_gt_nodules",
                    "max_malignancy", "mean_malignancy"])
        written = 0
        for c in targets:
            if c["case_id"] in failed_ids or not already_downloaded(c, args.out):
                continue
            ms = c["malignancies"]
            w.writerow([c["case_id"], c["patient_id"], c["n_nodules"], c["n_gt"],
                        max(ms) if ms else "", round(sum(ms) / len(ms), 2) if ms else ""])
            written += 1
    log(f"[4] 标签表: {label_path}（{written} 行，与磁盘图像对齐）")

    dt = (time.time() - t0) / 60
    log("\n" + "=" * 62)
    log(f"完成: 新下载 {done[0]}/{len(to_download)} 例，耗时 {dt:.1f} 分钟，"
        f"平均 {max(len(to_download),1)/max(dt,0.01):.1f} 例/分钟")
    if errors:
        log(f"失败 {len(errors)} 例（重跑本脚本自动续传）:")
        for pid, _, err in errors[:8]:
            log(f"  {pid}: {err}")
    log(f"数据目录: {os.path.abspath(args.out)}")
    log(f"下一步: python extract_nifti.py --root \"{os.path.abspath(args.out)}\"")


if __name__ == "__main__":
    main()
