# -*- coding: utf-8 -*-
"""
LIDC-IDRI 单例影像组学特征提取  v0.1
=========================================================
用途：把「一个病例的 CT + 4 位医师勾画」变成一张可用的特征表。
这一个脚本跑通之后，多例批处理只是外面套一层 for 循环。

用法：
    python extract_case.py "E:/lidc_data/lidc_idri/LIDC-IDRI-0001"

输出（写到 ./out/）：
    <PID>_features.csv   一位医师一行，一阶 + 形状 + GLCM 特征
    <PID>_report.png     三联图：CT+勾画 / HU 直方图 / 观察者一致性 Dice

依赖：numpy, pydicom, matplotlib（本机已具备，不需要 pyradiomics）
"""
import os
import sys
import glob
import json
import math
import time
import numpy as np
import pydicom


def _open_out(path, **kw):
    """打开输出文件；若被 Excel/编辑器占用，退回到带时间戳的新文件而不是崩掉。

    批处理 100 例时这个坑会反复出现——CSV 一被打开，整条流程就断。
    """
    try:
        return open(path, "w", **kw), path
    except PermissionError:
        stem, ext = os.path.splitext(path)
        alt = f"{stem}_{time.strftime('%Y%m%d_%H%M%S')}{ext}"
        print(f"      [!] {os.path.basename(path)} 被占用，改写到 {os.path.basename(alt)}")
        return open(alt, "w", **kw), alt

# ---------------------------------------------------------------- 基础 IO

def _find(patient_dir, modality):
    """按 modality 前缀（CT_/SEG_/SR_/DX_）找 dcm 文件"""
    return sorted(glob.glob(
        os.path.join(patient_dir, "**", f"{modality}_*", "**", "*.dcm"),
        recursive=True))


def load_ct(patient_dir):
    """读 CT 序列 -> (HU 体数据, z 坐标数组, 像素间距 mm)

    关键点：
      1. 必须按 ImagePositionPatient[2] 排序，不能信文件名或 InstanceNumber
      2. HU = 原始像素 * RescaleSlope + RescaleIntercept（LIDC 是 slope=1, intercept=-1024）
    """
    files = _find(patient_dir, "CT")
    if not files:
        raise FileNotFoundError(f"没找到 CT 序列: {patient_dir}")

    slices = []
    for f in files:
        ds = pydicom.dcmread(f, force=True)
        slices.append(ds)
    slices.sort(key=lambda d: float(d.ImagePositionPatient[2]))

    d0 = slices[0]
    slope = float(getattr(d0, "RescaleSlope", 1) or 1)
    inter = float(getattr(d0, "RescaleIntercept", 0) or 0)

    vol = np.stack([s.pixel_array.astype(np.float32) for s in slices], axis=0)
    vol = vol * slope + inter

    z = np.array([float(s.ImagePositionPatient[2]) for s in slices], dtype=np.float64)
    ps = [float(v) for v in d0.PixelSpacing]           # [row_spacing(y), col_spacing(x)]
    dz = float(np.median(np.diff(z))) if len(z) > 1 else float(d0.SliceThickness)
    spacing = (dz, ps[0], ps[1])                        # (z, y, x) mm
    return vol, z, spacing


def load_segs(patient_dir, z_ct):
    """读全部 SEG -> [{annotator, label, masks:{z_index: 2D bool}}]

    LIDC 的 SEG 是 DICOM-SEG：一个文件里 N 个 frame，每 frame 对应一层 CT。
    帧与 CT 层的对应关系靠 PerFrameFunctionalGroupsSequence 里的 ImagePositionPatient[2]，
    不要靠 frame 顺序猜。
    """
    files = _find(patient_dir, "SEG")
    out = []
    for f in files:
        ds = pydicom.dcmread(f, force=True)
        seglabel = "unknown"
        if "SegmentSequence" in ds and len(ds.SegmentSequence):
            seglabel = str(ds.SegmentSequence[0].SegmentLabel)

        # 医师代号藏在 SeriesDescription 或 SegmentLabel 里（形如 "Annotation IL057_127364"）
        annot = seglabel.split("Annotation")[-1].strip() or os.path.basename(os.path.dirname(os.path.dirname(f)))[-8:]

        px = ds.pixel_array                                # (nframes, H, W) uint8
        pf = ds.PerFrameFunctionalGroupsSequence
        masks = {}
        for i in range(px.shape[0]):
            zpos = float(pf[i].PlanePositionSequence[0].ImagePositionPatient[2])
            k = int(np.argmin(np.abs(z_ct - zpos)))        # 最近层
            if abs(z_ct[k] - zpos) > 1e-3:
                continue
            m = px[i] > 0
            if m.any():
                masks[k] = m
        if masks:
            out.append({"annotator": annot, "label": seglabel, "masks": masks})
    return out


# ---------------------------------------------------------------- 特征计算

def _first_order(vals):
    """ROI 内 HU 的一阶统计量"""
    v = vals.astype(np.float64)
    mean = v.mean()
    std = v.std(ddof=1) if v.size > 1 else 0.0
    z = (v - mean) / (std + 1e-9)
    n = v.size

    # 直方图熵（32 bin，基于 ROI 自身范围）
    hist, _ = np.histogram(v, bins=32, range=(v.min(), v.max() + 1e-9))
    p = hist.astype(np.float64) / max(hist.sum(), 1)
    nz = p[p > 0]
    entropy = float(-(nz * np.log2(nz)).sum())
    uniformity = float((p ** 2).sum())

    return {
        "n_voxels": int(n),
        "hu_mean": float(mean),
        "hu_std": float(std),
        "hu_min": float(v.min()),
        "hu_max": float(v.max()),
        "hu_p10": float(np.percentile(v, 10)),
        "hu_p25": float(np.percentile(v, 25)),
        "hu_median": float(np.median(v)),
        "hu_p75": float(np.percentile(v, 75)),
        "hu_p90": float(np.percentile(v, 90)),
        "hu_skewness": float((z ** 3).sum() / n),
        "hu_kurtosis": float((z ** 4).sum() / n - 3.0),
        "hu_entropy": entropy,
        "hu_uniformity": uniformity,
    }


def _shape(mask3d, spacing):
    """形状特征：体积、表面积、球形度、最大径"""
    dz, dy, dx = spacing
    voxel_vol = dz * dy * dx
    n = int(mask3d.sum())
    vol = n * voxel_vol

    # 表面积：统计与背景相邻的体素面（6 邻域）
    p = np.pad(mask3d, 1, mode="constant", constant_values=False)
    c = p[1:-1, 1:-1, 1:-1]
    faces = 0.0
    faces += (c & ~p[:-2, 1:-1, 1:-1]).sum() * (dy * dx)   # z 方向
    faces += (c & ~p[2:, 1:-1, 1:-1]).sum() * (dy * dx)
    faces += (c & ~p[1:-1, :-2, 1:-1]).sum() * (dz * dx)   # y 方向
    faces += (c & ~p[1:-1, 2:, 1:-1]).sum() * (dz * dx)
    faces += (c & ~p[1:-1, 1:-1, :-2]).sum() * (dz * dy)   # x 方向
    faces += (c & ~p[1:-1, 1:-1, 2:]).sum() * (dz * dy)
    surf = float(faces)

    # 球形度：与同体积球体的表面积之比，1 = 完美球
    sphericity = ((36 * math.pi * vol ** 2) ** (1 / 3) / surf) if surf > 0 else 0.0

    # 最大 3D 直径（体素多时随机降采样，O(n^2) 不可接受）
    zz, yy, xx = np.nonzero(mask3d)
    pts = np.column_stack([zz * dz, yy * dy, xx * dx])
    if len(pts) > 1500:
        idx = np.random.RandomState(42).choice(len(pts), 1500, replace=False)
        pts_s = pts[idx]
    else:
        pts_s = pts
    d2 = ((pts_s[:, None, :] - pts_s[None, :, :]) ** 2).sum(-1)
    max_diam = float(np.sqrt(d2.max())) if len(pts_s) > 1 else 0.0

    # 最大单层面积
    per_slice = mask3d.reshape(mask3d.shape[0], -1).sum(1)
    max_area = float(per_slice.max() * dy * dx)

    return {
        "volume_mm3": float(vol),
        "surface_mm2": surf,
        "sphericity": float(sphericity),
        "equiv_diameter_mm": float((6 * vol / math.pi) ** (1 / 3)),
        "max_3d_diameter_mm": max_diam,
        "max_axial_area_mm2": max_area,
        "n_slices": int((per_slice > 0).sum()),
        "elongation": float(max_diam / ((6 * vol / math.pi) ** (1 / 3))) if vol > 0 else 0.0,
    }


def _glcm_feats(vol, mask3d, levels=8):
    """简化 GLCM：轴位面内 4 个方向，量化到 8 级，取平均

    纯手工实现，不依赖 scikit-image —— 8 级灰度量化的 GLCM 已经是
    教科书级影像组学特征，够写进方法学章节。
    """
    vals = vol[mask3d]
    if vals.size < 16:
        return {}
    lo, hi = np.percentile(vals, [1, 99])
    if hi <= lo:
        hi = lo + 1.0
    q = np.clip(((vol - lo) / (hi - lo) * levels).astype(np.int16), 0, levels - 1)

    glcm = np.zeros((levels, levels), dtype=np.float64)
    nz = np.nonzero(mask3d)
    offsets = [(0, 0, 1), (0, 0, -1), (0, 1, 0), (0, -1, 0)]
    for o in offsets:
        nz2 = (nz[0] + o[0], nz[1] + o[1], nz[2] + o[2])
        ok = np.ones(len(nz[0]), dtype=bool)
        for d in range(3):                                   # 越界剔除
            ok &= (nz2[d] >= 0) & (nz2[d] < mask3d.shape[d])
        if not ok.any():
            continue
        a = nz[0][ok], nz[1][ok], nz[2][ok]
        b = nz2[0][ok], nz2[1][ok], nz2[2][ok]
        both = mask3d[b]                                     # 邻居也必须在 ROI 内
        if not both.any():
            continue
        ia = q[a][both]
        ib = q[b][both]
        np.add.at(glcm, (ia, ib), 1.0)

    s = glcm.sum()
    if s == 0:
        return {}
    glcm /= s
    i, j = np.meshgrid(np.arange(levels), np.arange(levels), indexing="ij")
    mu_i = (i * glcm).sum()
    mu_j = (j * glcm).sum()
    sd_i = math.sqrt(((i - mu_i) ** 2 * glcm).sum()) + 1e-9
    sd_j = math.sqrt(((j - mu_j) ** 2 * glcm).sum()) + 1e-9
    nz_p = glcm[glcm > 0]

    return {
        "glcm_contrast": float(((i - j) ** 2 * glcm).sum()),
        "glcm_homogeneity": float((glcm / (1.0 + (i - j) ** 2)).sum()),
        "glcm_energy": float((glcm ** 2).sum()),
        "glcm_correlation": float(((i - mu_i) * (j - mu_j) * glcm).sum() / (sd_i * sd_j)),
        "glcm_entropy": float(-(nz_p * np.log2(nz_p)).sum()),
    }


# ---------------------------------------------------------------- 一致性

def dice_matrix(segs, shape):
    """4 位医师两两 Dice —— LIDC 的经典观察者间一致性指标"""
    names = [s["annotator"] for s in segs]
    n = len(segs)
    d = np.eye(n)
    for a in range(n):
        for b in range(a + 1, n):
            va = vb = inter = 0
            za, zb = segs[a]["masks"], segs[b]["masks"]
            for k in set(za) | set(zb):
                ma = za.get(k, np.zeros(shape[1:], bool))
                mb = zb.get(k, np.zeros(shape[1:], bool))
                va += int(ma.sum())
                vb += int(mb.sum())
                inter += int((ma & mb).sum())
            d[a, b] = d[b, a] = (2 * inter / (va + vb)) if (va + vb) else 0.0
    return names, d


# ---------------------------------------------------------------- 主流程

def process(patient_dir, out_dir=None):
    pid = os.path.basename(patient_dir.rstrip("/\\"))
    out_dir = out_dir or os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
    os.makedirs(out_dir, exist_ok=True)

    print(f"[1/4] 读 CT ...")
    vol, z, spacing = load_ct(patient_dir)
    print(f"      {vol.shape[0]} 层 {vol.shape[1]}x{vol.shape[2]}, "
          f"间距 {spacing}, HU {vol.min():.0f}~{vol.max():.0f}")

    print(f"[2/4] 读 SEG ...")
    segs = load_segs(patient_dir, z)
    print(f"      {len(segs)} 位医师: {[s['annotator'] for s in segs]}")
    if not segs:
        raise RuntimeError("没有任何可用勾画，无法提取特征")

    print(f"[3/4] 提特征 ...")
    rows = []
    for s in segs:
        m3 = np.zeros(vol.shape, dtype=bool)
        for k, m2 in s["masks"].items():
            m3[k] |= m2
        if not m3.any():
            continue
        rec = {"patient_id": pid, "annotator": s["annotator"]}
        rec.update(_first_order(vol[m3]))
        rec.update(_shape(m3, spacing))
        rec.update(_glcm_feats(vol, m3))
        rows.append(rec)
        print(f"      {s['annotator']:>14s} 体素={rec['n_voxels']:>6d} "
              f"体积={rec['volume_mm3']:>9.1f}mm3 球形度={rec['sphericity']:.3f} "
              f"均值={rec['hu_mean']:>7.1f}HU")

    print(f"[4/4] 出图 + 存表 ...")
    names, dmat = dice_matrix(segs, vol.shape)
    _plot(vol, segs, names, dmat, os.path.join(out_dir, f"{pid}_report.png"))

    # CSV（不引 pandas，省一个依赖）
    csv_path = os.path.join(out_dir, f"{pid}_features.csv")
    keys = list(rows[0].keys())
    fh, csv_path = _open_out(csv_path, encoding="utf-8-sig", newline="")
    with fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(str(r.get(k, "")) for k in keys) + "\n")

    json_path = os.path.join(out_dir, f"{pid}_features.json")
    fh, json_path = _open_out(json_path, encoding="utf-8")
    with fh:
        json.dump({"patient_id": pid, "spacing": spacing,
                   "dice": {"names": names, "matrix": dmat.tolist()},
                   "records": rows}, fh, ensure_ascii=False, indent=2)

    print(f"\n完成：\n  {csv_path}\n  {json_path}\n  {out_dir}/{pid}_report.png")
    print(f"\n观察者一致性 Dice 均值（不含对角线）："
          f"{(dmat.sum() - len(names)) / (len(names) ** 2 - len(names)):.3f}")
    return rows, dmat


def _plot(vol, segs, names, dmat, path):
    """三联图：CT+勾画轮廓 / HU 直方图 / Dice 热图"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    # 找勾画面积最大的层做展示
    best_k, best_a = 0, -1
    for s in segs:
        for k, m in s["masks"].items():
            a = int(m.sum())
            if a > best_a:
                best_a, best_k = a, k

    fig, ax = plt.subplots(1, 3, figsize=(17, 5.6))
    ax[0].imshow(vol[best_k], cmap="gray", vmin=-1000, vmax=400)
    colors = ["#ff4d4f", "#36cfc9", "#ffc53d", "#9254de"]
    for i, s in enumerate(segs):
        if best_k in s["masks"]:
            ax[0].contour(s["masks"][best_k].astype(float), levels=[0.5],
                          colors=colors[i % 4], linewidths=1.4)
    ax[0].set_title(f"Axial slice #{best_k}  (contours: {len(segs)} annotators)")
    ax[0].axis("off")

    for i, s in enumerate(segs):
        vals = []
        for k, m in s["masks"].items():
            vals.append(vol[k][m])
        if vals:
            v = np.concatenate(vals)
            ax[1].hist(v, bins=60, range=(-1000, 400), alpha=0.45,
                       label=f"{s['annotator']} (n={v.size})",
                       color=colors[i % 4])
    ax[1].set_xlabel("HU"); ax[1].set_ylabel("voxels")
    ax[1].set_title("Intra-ROI HU distribution"); ax[1].legend(fontsize=8)

    im = ax[2].imshow(dmat, cmap="YlOrRd", vmin=0, vmax=1)
    ax[2].set_xticks(range(len(names)))
    ax[2].set_yticks(range(len(names)))
    short = [n[-10:] if len(n) > 10 else n for n in names]
    ax[2].set_xticklabels(short, rotation=45, ha="right", fontsize=8)
    ax[2].set_yticklabels(short, fontsize=8)
    for i in range(len(names)):
        for j in range(len(names)):
            ax[2].text(j, i, f"{dmat[i, j]:.2f}", ha="center", va="center", fontsize=9)
    ax[2].set_title("Inter-annotator Dice")
    fig.colorbar(im, ax=ax[2], fraction=0.046)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else r"E:/lidc_data/lidc_idri/LIDC-IDRI-0001"
    np.random.seed(42)
    process(target)
