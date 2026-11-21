# -*- coding: utf-8 -*-
"""
build_model.py — LIDC-IDRI 结节级良恶性分类，产出第一个可信 AUC

设计要点（每条都是被坑过才写下来的）：
  1) 特征筛选必须发生在每一折训练集内部。先在全表上选特征再交叉验证，
     等于把测试折的信息提前泄露，AUC 会虚高 0.05~0.15。
  2) 交叉验证按 patient_id 分组。同一病人的多枚结节若被劈到训练/测试两侧，
     CV 就变成"记住病人"而不是"学会结节"。
  3) 灰区（malignancy == 3）不进二分类，单独统计后剔除。
  4) 给三个参照：直径、体积、全特征不筛选。没有参照的 AUC 无法判断是否有价值。

用法：
    python build_model.py                        # 默认：全结节、5折x5次重复、100次置换
    python build_model.py --min-diameter 6       # 只看 >=6mm
    python build_model.py --n-perm 0             # 跳过置换检验（快）
"""
import os
import json
import argparse
import warnings

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score, roc_curve, brier_score_loss

warnings.filterwarnings("ignore")

RNG = np.random.default_rng(20261115)
C_GRID = np.logspace(-3, 1, 9)      # L1 正则强度候选


# ---------------------------------------------------------------- 数据
def load_data(csv_path, min_diam=0.0, max_diam=0.0, endpoint="mal"):
    df = pd.read_csv(csv_path, encoding="utf-8-sig")
    m = pd.to_numeric(df["malignancy_score"], errors="coerce")
    d = pd.to_numeric(df["diameter_mm"], errors="coerce")
    ok = m.notna() & (m >= 1) & (m <= 5)
    if min_diam > 0:
        ok &= d >= min_diam
    if max_diam > 0:
        ok &= d < max_diam
    n_all = len(df)
    df = df[ok].reset_index(drop=True)
    n_grey = int((pd.to_numeric(df.malignancy_score) == 3).sum())
    if endpoint == "mal":                     # 原终点：灰区不进二分类
        df = df[pd.to_numeric(df.malignancy_score) != 3].reset_index(drop=True)
    # grey 终点：全部保留，灰区本身就是正类


    score = pd.to_numeric(df.malignancy_score)
    if endpoint == "grey":
        y = (score == 3).astype(int).to_numpy()
    else:
        y = (score >= 4).astype(int).to_numpy()

    feats = [c for c in df.columns if c.startswith("original_")]
    X = df[feats].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    return df, X, y, feats, dict(n_all=n_all, n_grey=n_grey)


# ---------------------------------------------------------------- 内层选特征
def inner_select(Xtr, ytr, gtr, seed):
    """在训练折内部用 3 折分组 CV 选 L1 强度，返回 (best_C, mask)。"""
    n_min = min(np.bincount(ytr))
    n_grp = len(np.unique(gtr))
    if n_min < 3 or n_grp < 6:
        # 折内样本太薄，退化为弱正则、全特征
        return 1.0, np.ones(Xtr.shape[1], dtype=bool)

    skf = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=seed)
    scores = []
    for C in C_GRID:
        aucs = []
        for itr, iva in skf.split(Xtr, ytr, gtr):
            if len(np.unique(ytr[itr])) < 2 or len(np.unique(ytr[iva])) < 2:
                continue
            sc = StandardScaler().fit(Xtr[itr])
            lr = LogisticRegression(penalty="l1", solver="liblinear", C=C,
                                    class_weight="balanced", max_iter=3000)
            lr.fit(sc.transform(Xtr[itr]), ytr[itr])
            p = lr.predict_proba(sc.transform(Xtr[iva]))[:, 1]
            aucs.append(roc_auc_score(ytr[iva], p))
        scores.append(np.mean(aucs) if aucs else 0.5)

    best_C = float(C_GRID[int(np.argmax(scores))])
    sc = StandardScaler().fit(Xtr)
    lr = LogisticRegression(penalty="l1", solver="liblinear", C=best_C,
                            class_weight="balanced", max_iter=3000)
    lr.fit(sc.transform(Xtr), ytr)
    mask = np.abs(lr.coef_[0]) > 1e-8
    if mask.sum() == 0:                       # LASSO 全清零 -> 保住系数最大的一项
        mask[np.argmax(np.abs(lr.coef_[0]))] = True
    return best_C, mask


# ---------------------------------------------------------------- 外层 CV
def run_cv(X, y, groups, feats, repeats=5, seed=42, select=True, fixed_mask=None,
           verbose=True):
    n = len(y)
    oof_stack, sel_count = [], np.zeros(X.shape[1])
    fold_rows = []
    for rep in range(repeats):
        skf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed + rep)
        oof = np.full(n, np.nan)
        for k, (tr, te) in enumerate(skf.split(X, y, groups), 1):
            if select:
                C, mask = inner_select(X[tr], y[tr], groups[tr], seed + rep * 10 + k)
                sel_count[mask] += 1
            else:
                C, mask = 1.0, (fixed_mask if fixed_mask is not None
                                else np.ones(X.shape[1], dtype=bool))
            sc = StandardScaler().fit(X[tr][:, mask])
            lr = LogisticRegression(penalty="l2", solver="lbfgs", C=1.0,
                                    class_weight="balanced", max_iter=3000)
            lr.fit(sc.transform(X[tr][:, mask]), y[tr])
            oof[te] = lr.predict_proba(sc.transform(X[te][:, mask]))[:, 1]
            fold_rows.append(dict(rep=rep, fold=k, C=C, n_feat=int(mask.sum()),
                                  n_te=len(te)))
        if np.isnan(oof).any():
            if verbose:
                print(f"  [warn] rep{rep} 有未覆盖样本")
        else:
            oof_stack.append(oof)
    aucs = [roc_auc_score(y, o) for o in oof_stack]
    return dict(aucs=np.array(aucs), oof=oof_stack, sel_count=sel_count,
                folds=pd.DataFrame(fold_rows))


def boot_ci(y, p, groups, n_boot=2000, seed=7):
    """按病人整块重抽样的聚类 bootstrap；单类重抽样跳过。"""
    rng = np.random.default_rng(seed)
    uniq = np.unique(groups)
    idx_by = {g: np.where(groups == g)[0] for g in uniq}
    out = []
    for _ in range(n_boot):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([idx_by[g] for g in pick])
        yy = y[idx]
        if len(np.unique(yy)) < 2:
            continue
        out.append(roc_auc_score(yy, p[idx]))
    out = np.array(out)
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def paired_delta(y, p_new, p_old, groups, n_boot=2000, seed=11):
    """同一批病人重抽样下算 AUC 差值的分布 -> 增量价值是否真的 > 0。"""
    rng = np.random.default_rng(seed)
    uniq = np.unique(groups)
    idx_by = {g: np.where(groups == g)[0] for g in uniq}
    out = []
    for _ in range(n_boot):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([idx_by[g] for g in pick])
        if len(np.unique(y[idx])) < 2:
            continue
        out.append(roc_auc_score(y[idx], p_new[idx]) - roc_auc_score(y[idx], p_old[idx]))
    out = np.array(out)
    return float(out.mean()), float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def univariate_auc(x, y, groups, name):
    x = np.asarray(x, dtype=float)
    if name.startswith("log_"):
        x = np.log1p(np.clip(x, 0, None))
    auc = roc_auc_score(y, x)
    lo, hi = boot_ci(y, x, groups)
    return dict(name=name, auc=float(auc), ci=(lo, hi))


# ---------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="features/nodules_features.csv")
    ap.add_argument("--out", default="out_model")
    ap.add_argument("--min-diameter", type=float, default=0.0)
    ap.add_argument("--max-diameter", type=float, default=0.0)
    ap.add_argument("--add-size", action="store_true",
                    help="把 log(直径) 作为特征强制放进模型，检验纹理是否带来增量")
    ap.add_argument("--endpoint", choices=["mal", "grey"], default="mal",
                    help="mal=良恶性分类（默认，原逻辑）；grey=灰区分歧终点（3分 vs 确信）")
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--n-perm", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    df, X, y, feats, info = load_data(args.csv, args.min_diameter, args.max_diameter,
                                      endpoint=args.endpoint)
    if args.add_size:
        X = np.column_stack([X, np.log1p(pd.to_numeric(df.diameter_mm).to_numpy())])
        feats = feats + ["SIZE_log_diameter"]
    groups = df["patient_id"].astype(str).to_numpy()
    n, n_mal, n_ben = len(y), int(y.sum()), int((1 - y).sum())
    n_pat = df["patient_id"].nunique()
    tag = "灰区(3分)" if args.endpoint == "grey" else "恶性(4-5分)"
    print(f"样本: {n} 结节 / {n_pat} 病人 | 正类[{tag}] {n_mal} : 负类 {n_ben} "
          f"| 灰区已剔 {info['n_grey']}（原始 {info['n_all']} 行）")
    print(f"EPV(10) 可容纳特征上限: {min(n_mal, n_ben)//10}")

    # ---- 参照模型（无筛选，单变量）
    base = []
    for col, tag in [("diameter_mm", "log_diameter"), ("gt_volume_mm3", "log_volume")]:
        if col in df.columns:
            base.append(univariate_auc(df[col].to_numpy(), y, groups, tag))
    # 全特征不筛选（反例）
    naive = run_cv(X, y, groups, feats, repeats=1, seed=args.seed, select=False)
    base.append(dict(name="all107_noSelect", auc=float(naive["aucs"][0]),
                     ci=boot_ci(y, naive["oof"][0], groups)))

    # ---- 主模型：LASSO 嵌套
    print("\n主模型运行中（嵌套 CV）...")
    cv = run_cv(X, y, groups, feats, repeats=args.repeats, seed=args.seed, select=True)
    auc_mean, auc_sd = float(cv["aucs"].mean()), float(cv["aucs"].std())
    aucs_txt = " / ".join(f"{a:.3f}" for a in cv["aucs"])
    pooled = np.concatenate(cv["oof"]) if cv["oof"] else np.array([])
    y_pool = np.tile(y, len(cv["oof"]))
    g_pool = np.tile(groups, len(cv["oof"]))
    lo, hi = boot_ci(y_pool, pooled, g_pool) if len(pooled) else (np.nan, np.nan)
    brier = float(brier_score_loss(y_pool, pooled)) if len(pooled) else np.nan

    delta = None
    if args.add_size and len(pooled):
        p_size = np.tile(np.log1p(pd.to_numeric(df.diameter_mm).to_numpy()),
                         len(cv["oof"]))
        delta = paired_delta(y_pool, pooled, p_size, g_pool)

    n_folds = max(1, args.repeats * 5)
    freq = pd.DataFrame({"feature": feats, "selected": cv["sel_count"] / n_folds})
    freq = freq[freq.selected > 0].sort_values("selected", ascending=False)
    freq.to_csv(os.path.join(args.out, "selected_features.csv"),
                index=False, encoding="utf-8-sig")

    # ---- 置换检验
    p_perm = np.nan
    perm_aucs = []
    if args.n_perm > 0:
        print(f"置换检验 {args.n_perm} 次（这步最慢）...")
        rng = np.random.default_rng(args.seed + 999)
        for i in range(args.n_perm):
            yp = rng.permutation(y)
            if len(np.unique(yp)) < 2:
                continue
            r = run_cv(X, yp, groups, feats, repeats=1, seed=args.seed, select=True,
                       verbose=False)
            perm_aucs.append(float(r["aucs"][0]))
            if (i + 1) % 20 == 0:
                print(f"  {i+1}/{args.n_perm} 完成，当前零分布均值 {np.mean(perm_aucs):.3f}")
        if perm_aucs:
            p_perm = (1 + sum(a >= auc_mean for a in perm_aucs)) / (1 + len(perm_aucs))

    # ---- 病例级敏感性分析（同一病例取结节预测最大值）
    case_auc = np.nan
    case_n = 0
    if cv["oof"]:
        p0 = cv["oof"][0]
        d = pd.DataFrame({"case": df.case_id.values, "pat": groups, "p": p0, "y": y})
        cg = d.groupby("case").agg(p=("p", "max"), y=("y", "max"), pat=("pat", "first"))
        case_n = len(cg)
        case_auc = float(roc_auc_score(cg.y, cg.p))

    # ---- 落盘
    if cv["oof"]:
        pd.DataFrame({"case_id": df.case_id, "patient_id": groups,
                      "malignancy_score": df.malignancy_score, "y": y,
                      "oof_prob": cv["oof"][0]}).to_csv(
            os.path.join(args.out, "oof_predictions.csv"), index=False,
            encoding="utf-8-sig")

    lines = []
    A = lines.append
    A("=" * 66)
    A("LIDC-IDRI 结节级良恶性分类 · 首个 AUC")
    A("=" * 66)
    A(f"样本        {n} 结节 / {n_pat} 病人 | 恶 {n_mal} : 良 {n_ben} "
      f"= {n_mal/max(n_ben,1):.2f} : 1")
    A(f"灰区剔除    {info['n_grey']} 例（malignancy = 3）")
    A(f"直径门槛    {'无' if args.min_diameter <= 0 else f'>= {args.min_diameter} mm'}"
      f"{'' if args.max_diameter <= 0 else f' 且 < {args.max_diameter} mm'}")
    if args.add_size:
        A("模型输入    107 纹理/形状特征 + log(直径)（强制纳入，用于检验增量）")
    A("")
    A("【主模型】LASSO-L1 逻辑回归，嵌套 CV（外层 5 折 × "
      f"{args.repeats} 次重复，按 patient_id 分组）")
    A(f"  AUC       {auc_mean:.3f} ± {auc_sd:.3f}   各次: {aucs_txt}")
    A(f"  95% CI    {lo:.3f} – {hi:.3f}（病人级聚类 bootstrap，2000 次）")
    A(f"  Brier     {brier:.3f}（0.25 = 无信息基线）")
    A(f"  每折入选   {cv['folds'].n_feat.mean():.1f} ± {cv['folds'].n_feat.std():.1f} 个特征")
    if not np.isnan(p_perm):
        A(f"  置换检验  p = {p_perm:.3f}（{len(perm_aucs)} 次，零分布 "
          f"{np.mean(perm_aucs):.3f} ± {np.std(perm_aucs):.3f}）")
    if delta:
        A(f"  增量(纹理-直径) ΔAUC = {delta[0]:+.4f}  95% CI {delta[1]:+.4f} ~ {delta[2]:+.4f}"
          f"  -> {'有增量' if delta[1] > 0 else '无增量（CI 跨 0）'}")
    A("")
    A("【参照】")
    for b in base:
        A(f"  {b['name']:<16s} AUC {b['auc']:.3f}  ({b['ci'][0]:.3f} – {b['ci'][1]:.3f})")
    A("")
    A("【病例级敏感性】同一病例取各结节预测最大值")
    A(f"  {case_n} 例 → AUC {case_auc:.3f}" if case_n else "  跳过")
    A("")
    A("【入选最稳定的特征（入选折数占比）】")
    for _, r in freq.head(12).iterrows():
        A(f"  {r.feature:<44s} {r.selected*100:5.0f}%")
    A("")
    ok_ci = (not np.isnan(lo)) and lo > 0.60
    ok_perm = np.isnan(p_perm) or p_perm < 0.05
    verdict = ("可用：CI 下限 > 0.60 且优于单变量基线" if (ok_ci and ok_perm and auc_mean > max(
        b["auc"] for b in base)) else "可疑：未同时满足 CI 下限 > 0.60、置换 p < 0.05、优于基线")
    A(f"【判定】{verdict}")
    A("=" * 66)

    report = "\n".join(lines)
    print("\n" + report)
    with open(os.path.join(args.out, "report.txt"), "w", encoding="utf-8") as f:
        f.write(report + "\n")

    json.dump(dict(n=n, n_pat=n_pat, n_mal=n_mal, n_ben=n_ben,
                   auc_mean=auc_mean, auc_sd=auc_sd, ci=[lo, hi],
                   aucs=cv["aucs"].tolist(), brier=brier, p_perm=p_perm,
                   delta_over_size=delta,
                   case_auc=case_auc, case_n=case_n,
                   baselines=base,
                   mean_selected=float(cv["folds"].n_feat.mean())),
              open(os.path.join(args.out, "metrics.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)

    # ---- ROC 图
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5.6, 5.2), dpi=140)
        if len(pooled):
            fpr, tpr, _ = roc_curve(y_pool, pooled)
            ax.plot(fpr, tpr, lw=2, color="#c0392b",
                    label=f"LASSO-LR  AUC={auc_mean:.3f} (95%CI {lo:.2f}-{hi:.2f})")
        for b in base:
            ax.plot([], [], " ", label=f"{b['name']}: {b['auc']:.3f}")
        ax.plot([0, 1], [0, 1], "--", lw=1, color="#7f8c8d")
        ax.set_xlabel("1 - Specificity")
        ax.set_ylabel("Sensitivity")
        ax.set_title(f"LIDC-IDRI nodule-level  n={n} ({n_mal} mal / {n_ben} ben)")
        ax.legend(loc="lower right", fontsize=8, frameon=False)
        ax.set_aspect("equal")
        fig.tight_layout()
        fig.savefig(os.path.join(args.out, "roc.png"))
        print(f"ROC 图: {os.path.join(args.out, 'roc.png')}")
    except Exception as e:
        print("绘图跳过:", e)


if __name__ == "__main__":
    main()
