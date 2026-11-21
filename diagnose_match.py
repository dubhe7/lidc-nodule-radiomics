# -*- coding: utf-8 -*-
"""匹配质量核验：只看掩膜连通域 vs GT 质心的距离分布，不提特征。

用途：判定 --match-tol-mm 默认值 5.0 是否过紧，以及未匹配病例的成因
（容差边缘 vs 掩膜碎裂 vs 真缺失）。
"""
import os
import numpy as np
import pandas as pd
import SimpleITK as sitk
from scipy import ndimage as ndi
from scipy.optimize import linear_sum_assignment
from multiprocessing import Pool

ROOT = r"E:/lidc_data/lidc_hf"
MIN_VOX, MIN_MM = 20, 5.0


def one(cid):
    img = sitk.ReadImage(os.path.join(ROOT, "images", cid + ".nii.gz"))
    msk = sitk.ReadImage(os.path.join(ROOT, "masks", cid + ".nii.gz"))
    sp = np.array([img.GetSpacing()[2], img.GetSpacing()[1], img.GetSpacing()[0]])
    m = sitk.GetArrayFromImage(msk) > 0
    lab, n = ndi.label(m, structure=np.ones((3, 3, 3), bool))
    comps = []
    for k in range(1, n + 1):
        mm = lab == k
        s = int(mm.sum())
        if s < MIN_VOX:
            continue
        comps.append((ndi.center_of_mass(mm), s))
    return cid, sp, comps


def main():
    nd = pd.read_csv(os.path.join(ROOT, "nodules.csv"), encoding="utf-8-sig")
    nd = nd[nd["in_default_gt"].astype(str).str.lower() == "true"].copy()
    nd["diameter"] = pd.to_numeric(nd["diameter_mm"], errors="coerce")
    nd = nd[nd["diameter"] >= MIN_MM]
    gt = {c: g[["nodule_id", "centroid_z", "centroid_y", "centroid_x", "diameter_mm"]].values
          for c, g in nd.groupby("case_id")}
    # 只核验已下载到盘的病例（nodules.csv 覆盖全库 703 例，本队列是分层抽样 302 例）
    cases = sorted([c for c in gt if os.path.exists(
        os.path.join(ROOT, "images", c + ".nii.gz"))])
    print(f"GT 病例 {len(gt)}，已下载 {len(cases)}")

    with Pool(6) as p:
        res = p.map(one, cases, chunksize=1)

    rows, unmatched = [], []
    for cid, sp, comps in res:
        g = gt[cid]
        if not comps or not len(g):
            unmatched.append((cid, "no_comp" if not comps else "no_gt", len(comps), len(g)))
            continue
        C = np.array([c for c, _ in comps]) * sp
        G = np.array([[float(r[1]), float(r[2]), float(r[3])] for r in g]) * sp
        d = np.sqrt(((C[:, None, :] - G[None, :, :]) ** 2).sum(-1))
        ri, ci = linear_sum_assignment(d)
        for i, j in zip(ri, ci):
            rows.append({"case_id": cid, "nodule_id": int(g[j][0]), "dist_mm": d[i, j],
                         "vox": comps[i][1]})
        # GT 侧被 assignment 分配到但距离过大 = 真未匹配，记录最近的组件距离
        for j in range(len(g)):
            ii = np.where(ci == j)[0]
            dd = d[ii[0], j] if len(ii) else float("inf")
            if dd > 5.0:
                unmatched.append((cid, int(g[j][0]), round(float(dd), 2),
                                  f"diam{g[j][4]}mm"))

    df = pd.DataFrame(rows)
    print(f"病例 {len(cases)} | 配对 {len(df)}")
    print("距离 mm 分位：", {q: round(df.dist_mm.quantile(q), 2) for q in (.5, .9, .99, 1.0)})
    for t in (5, 6, 8, 10):
        print(f"  tol={t}mm -> 命中 {int((df.dist_mm <= t).sum())} 对，"
              f"距离>{t}mm 的 {int((df.dist_mm > t).sum())} 对")
    print("\n超出 5mm 的配对（含真未匹配）：")
    for u in unmatched:
        print("  ", u)


if __name__ == "__main__":
    main()
