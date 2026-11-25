# Zenodo DOI 发布操作指南（当前状态版，约 10 分钟）

> 更新于 2026-11-25 20:00。仓库侧动作（push、v1.0.0 release）已由 AI 代办完成，**只剩 Zenodo 账号侧的两步**——因为需要 GitHub 账号本人 OAuth 登录，任何人代办不了。

## 当前已完成（不用再动）

- [x] 代码与文档全部 push 到 `https://github.com/dubhe7/lidc-nodule-radiomics`（最新 commit `130e2c2`）
- [x] `CITATION.cff` / `RELEASE_NOTES_v1.0.0.md` / README 披露段 已在仓库中
- [x] **v1.0.0 release 已发布**：https://github.com/dubhe7/lidc-nodule-radiomics/releases/tag/v1.0.0 （含复现性修复，快照为 `130e2c2`）

## 第 1 步：登录 Zenodo 并开启仓库归档（一次性，3 分钟）

1. 打开 https://zenodo.org → 右上角 **Log in** → 选 **Log in with GitHub**（dubhe7 账号 OAuth 授权）。
2. 首次登录会要求补全资料（姓名、邮箱、机构可填 桂林医学院 / Guilin Medical University），保存。
3. 进 https://zenodo.org/account/settings/github/ → 页面列出你的仓库 → 找到 `dubhe7/lidc-nodule-radiomics` → 把开关拨到 **ON**。

开启后 Zenodo 会自动抓取**已有的 v1.0.0 release**（无需再打一次 release）。

## 第 2 步：拿 DOI（1–5 分钟）

1. 打开 https://zenodo.org/me/uploads （或 My dashboard），等出现 `lidc-nodule-radiomics` 条目，状态 `Published`。
2. 点进去，右侧 DOI 区有两个号：
   - **版本 DOI**（本版专用，形如 `10.5281/zenodo.XXXXXXX`）
   - **Concept DOI**（跨版本通用，README/CITATION 里用它）
3. 把 **Concept DOI** 复制发我，我回填 `CITATION.cff` 与 `README.md` 两处占位并 push。

## 验收标准

- [ ] Zenodo dashboard 出现条目，状态 Published，有 DOI
- [ ] Zenodo 条目文件列表里能看到代码文件、README、CITATION.cff
- [ ] DOI 链接可在无痕窗口打开
- [ ] DOI 已回填进仓库（AI 执行）

## 署名说明

当前 CITATION.cff 署名为 `dubhe7`。若日后需要真名（中文核心投稿、简历），改 CFF 重打 release 即可——版本 DOI 会变，**Concept DOI 不变**，引用连续性不受影响。唯一规则：简历/投稿中出现的名字必须与 DOI 上的一致。
