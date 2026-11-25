# Zenodo DOI 发布操作指南（作者执行，约 15 分钟）

> 分工：以下 1–5 步只能由 GitHub 账号本人操作（需要登录态）。做之前先把新文件 push 上去（命令见文末）。

## 前置：先把本批新文件推上仓库

```bash
cd "E:/新建文件夹 (4)/2026-11-20-20-27-11/output/lidc-radiomics-repo"
git add README.md CITATION.cff RELEASE_NOTES_v1.0.0.md ZENODO_STEPS.md
git commit -m "docs: add citation metadata, v1.0.0 release notes, AI disclosure, zenodo guide"
git push
```

## 1. Zenodo 关联 GitHub（一次性）

1. 浏览器打开 https://zenodo.org ，点右上角 **Sign up / Log in** → 选 **Log in with GitHub**（用 dubhe7 账号 OAuth 授权）。
2. 登录后进 https://zenodo.org/account/settings/github/ ，找到 `dubhe7/lidc-nodule-radiomics`，点 **Enable**（开启后该仓库的每个 release 都会自动归档）。

## 2. 打 release 触发归档

1. 打开 https://github.com/dubhe7/lidc-nodule-radiomics/releases → **Draft a new release**。
2. Choose a tag：输入 `v1.0.0`（main 分支最新提交上新建）。
3. Title：`v1.0.0 — 首次公开归档`。
4. Describe：把 `RELEASE_NOTES_v1.0.0.md` 全文粘贴进去。
5. **Publish release**。Zenodo 会在几分钟内自动抓取仓库快照。

## 3. 拿 DOI

1. 回 https://zenodo.org/uploads/ （或 My dashboard），等状态变红条目出现（通常 1–5 分钟）。
2. 点进去，右侧 **Cite all versions?** / DOI 区拿到两个号：
   - **版本 DOI**（v1.0.0 专用）
   - **概念 DOI**（代表所有版本，推荐写进 README/CITATION.cff 的是它）
3. 复制 DOI，形如 `10.5281/zenodo.XXXXXXX`。

## 4. 回填 DOI（两处）

1. `CITATION.cff`：取消 `doi:` 行注释，填入真实 DOI。
2. `README.md` 引用段：把 `10.5281/zenodo.XXXXXXX` 占位替换为真实 DOI。
3. 再 commit + push 一次（不影响已归档快照，README 属于"活文档"）。

## 5. 验收标准

- [ ] Zenodo dashboard 出现条目，状态为 published，有 DOI；
- [ ] GitHub Releases 页显示 v1.0.0 + Release notes；
- [ ] Zenodo 条目里能看到 README、CITATION.cff、代码文件（不含 results 图以外的多余大文件即可）；
- [ ] DOI 链接在新标签页打开能访问。

## 署名决定（唯一待你拍板的项）

当前 CITATION.cff 署名为 `dubhe7`（假名）。若中文核心投稿阶段需要真名，届时改 CITATION.cff 重新打 release 即可（Zenodo 会生成新版本 DOI，概念 DOI 不变）。两种都可以，没有对错，只有一致性问题：**简历/投稿里出现的名字必须和 DOI 上的一致**。
