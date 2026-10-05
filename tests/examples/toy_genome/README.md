# 小型合成参考基因组

`toy_uprobe.fa` 和 `toy_uprobe.gtf` 是相互匹配的合成数据，仅用于功能测试，不能用于真实生物学分析。固定随机种子：20261004。FASTA 使用染色体名 `1`、`2`，共 20,000 bp；GTF 使用标准的 1-based inclusive 坐标、九列制表符分隔以及 gene_id、gene_name、transcript_id 属性。

| 基因 | 染色体 | 链 | 外显子数 | 标准转录本长度 |
| --- | --- | --- | --- | --- |
| ToyGeneA | 1 | + | 3 | 2700 bp |
| ToyGeneB | 1 | - | 3 | 2700 bp |
| ToyGeneC | 2 | + | 2 | 2400 bp |

每个基因包含两个 UTR 注释；未添加 CDS，因此后端的 CDS/exon 回退逻辑会选择 exon。`targets.txt` 包含三个可使用的 gene_name。

## 已部署的 WSL 数据

Ubuntu WSL 路径：`/home/qzhang/uprobe_all/genomes/public/toy_uprobe`。

已在 `/home/qzhang/uprobe_all/cathe/genomes.yaml` 注册 `toy_uprobe`；原配置备份为 `genomes.yaml.before_toy_uprobe.bak`。索引配置只启用 Bowtie2，并关闭 Jellyfish。已生成基因组和转录本 Bowtie2 索引，分别位于 `bowtie2_genome/toy_uprobe`、`bowtie2_transcript/toy_uprobe`。

前端刷新 Genome 页面后选择 `toy_uprobe`。创建任务时选择相同基因组，基因目标填写 `ToyGeneA`、`ToyGeneB`、`ToyGeneC`。DNA 区间测试可使用 `1:1001-4900` 或 `2:1001-4200`。本地副本也可通过 Genome 页上传到另一个私人测试基因组；上传后后端会自动登记 WSL 文件路径。

## 验证范围

使用正在运行后端的 Python 环境验证了基因名称、8 个外显子和3条转录本提取，并生成了 FASTA 索引。未提交完整探针设计任务；完整任务仍取决于所选协议和工具依赖。

现有后端使用 GTF start 直接作为 Python 切片起点，每个外显子少取1 bp，实际输出转录本长度为2697、2697、2398 bp，详见 `validation.json`。数据保持标准 GTF 格式，未为此修改坐标或后端代码。

`prepare_wsl.py` 可在后端环境、后端项目目录下生成并登记数据。为避免覆盖已有数据，目标目录已存在时会退出。
