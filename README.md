# bulk RNA-seq Snakemake Workflow

该仓库提供一个完整的 bulk RNA-seq 分析流程，执行顺序如下：

1. FastQC（原始数据）
2. fastp
3. FastQC（清洗后）
4. MultiQC（阶段1）
5. STAR
6. SAMtools（索引）
7. Salmon
8. RSeQC + Qualimap
9. MultiQC（阶段2）
10. DESeq2（PCA + DEG）

## 目录结构

- `Snakefile`：主流程定义
- `config/config.yaml`：流程参数
- `config/samples.tsv`：样本 FASTQ 列表
- `config/metadata.tsv`：差异分析分组信息
- `config/tx2gene.tsv`：转录本到基因映射
- `envs/*.yaml`：conda 环境
- `scripts/run_deseq2.R`：DESeq2 分析脚本

## 使用方法

```bash
# 1) 创建 conda 环境并 dry-run
snakemake -n -p --use-conda

# 2) 正式运行
snakemake --cores 16 --use-conda
```

## 配置说明

请至少修改以下路径：

- `star_index`
- `salmon_index`
- `gtf`
- `rseqc_bed`
- `samples`
- `metadata`
- `tx2gene`

`contrasts` 的格式必须为：`因子列:实验组:对照组`，例如：

```yaml
contrasts:
  - "condition:treatment:control"
```

## 主要输出

- `results/multiqc_qc/multiqc_report.html`
- `results/multiqc_final/multiqc_report.html`
- `results/deseq2/PCA_plot.pdf`
- `results/deseq2/*_DEG.tsv`

