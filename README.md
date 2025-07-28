# mtDNApipe

MitoVariantPipe 是一个用于从二代测序（NGS）数据中检测线粒体DNA（mtDNA）低频异质性（low-frequency heteroplasmy）的生物信息学分析流程。

该流程整合了质量控制、序列比对、Indel重比对以及一套定制的脚本，用于处理测序错误、PCR扩增偏好，并最终实现对低至0.1%频率的突变进行高精度检测。

## 流程概览

```
Raw FASTQ -> Quality Filtering (Custom Script) -> Adapter/Quality Trimming (fastp) -> Overlap Correction (Custom C++) -> Alignment (BWA) -> BAM Processing (Picard & GATK) -> mtDNA Read Extraction -> Custom Deduplication -> Variant Calling
```

## 1. 先决条件 (Prerequisites)

### 软件依赖
本流程依赖于一系列标准的生物信息学工具。我们强烈建议使用 `Conda` 来管理这些依赖。

-   Python (>= 3.7)
-   BWA (0.7.17)
-   Fastp
-   Samtools
-   Picard
-   GATK (3.8) - **注意**: 本流程使用 GATK 3.x 版本的语法。
-   Perl
-   A C++ compiler (e.g., g++)

### 参考基因组


## 2. 安装 (Installation)

我们推荐使用 Conda 来创建独立的运行环境，这可以避免与系统其他软件的冲突。

1.  **克隆仓库**
    ```bash
    git clone https://github.com/your-username/MitoVariantPipe.git
    cd MitoVariantPipe
    ```

2.  **创建并激活 Conda 环境**
    使用项目提供的 `environment.yml` 文件来自动安装所有软件依赖。
    ```bash
    conda env create -f environment.yml
    conda activate mitovariantpipe
    ```
    这个过程可能需要几分钟。激活环境后，所有需要的工具（如 `bwa`, `fastp`）都将在你的 `PATH` 中。

3.  **编译 C++ 脚本**
    本流程包含一个用于校正 Read Overlap 的 C++ 脚本，需要手动编译。
    ```bash
    cd scripts/overlap_corrector/
    make
    cd ../../ # 返回项目根目录
    ```
    编译成功后，会在 `scripts/overlap_corrector/` 目录下生成一个名为 `overlap_corrector` 的可执行文件。

4.  **准备参考文件**
    -   将你的线粒体参考基因组 `hg19_mt.fa`（或类似文件）放置在一个方便访问的位置。
    -   检查 `ref/chrM_refAllele.txt` 文件是否存在。

至此，安装完成！

## 3. 使用方法 (Usage)

所有分析步骤都已封装在 `run_pipeline.sh` 脚本中。

### 命令格式
```bash
bash run_pipeline.sh -i <sample_id> \
                     -1 <path/to/read1.fq.gz> \
                     -2 <path/to/read2.fq.gz> \
                     -r <path/to/ref_genome.fa> \
                     -o <path/to/output_directory>
```

### 参数说明
-   `-i` : **[必需]** 样本ID，例如 `Sample01`。
-   `-1` : **[必需]** Read 1 的 FASTQ 文件路径 (`.fq.gz`)。
-   `-2` : **[必需]** Read 2 的 FASTQ 文件路径 (`.fq.gz`)。
-   `-r` : **[必需]** 线粒体参考基因组的 FASTA 文件路径。
-   `-o` : **[必需]** 输出目录的路径。所有结果将保存在 `<output_directory>/<sample_id>/` 下。
-   `-t` : **[可选]** BWA比对时使用的线程数 (默认为 8)。
-   `-h` : 显示帮助信息。

### 运行测试示例
为了验证安装是否成功，你可以运行 `example/` 目录下的测试数据。
```bash
# 确保你已经激活了 conda 环境
# conda activate mitovariantpipe

# 运行测试脚本
# 注意：你需要将 /path/to/your/hg19_mt.fa 替换为真实路径
bash example/run_test.sh /path/to/your/hg19_mt.fa
```
`run_test.sh` 的内容如下：
```bash
#!/bin/bash
# A script to test the pipeline with example data.
# Usage: bash run_test.sh <path_to_your_reference_genome>

REF_GENOME=$1
if [ -z "$REF_GENOME" ]; then
    echo "Error: Please provide the path to your reference genome."
    exit 1
fi

bash ../run_pipeline.sh -i test_sample \
                        -1 test_R1.fq.gz \
                        -2 test_R2.fq.gz \
                        -r $REF_GENOME \
                        -o ../test_output
```

## 4. 输入与输出文件

### 输入文件
1.  **双端测序数据**: `_R1.fq.gz` 和 `_R2.fq.gz`。
2.  **线粒体参考基因组**: FASTA 格式。

### 主要输出文件
所有输出文件均位于 `<output_dir>/<sample_id>/` 目录下。

-   `*_analysis.log`: 整个流程的运行日志。
-   `*.realn.bam`: GATK Indel 重比对后的最终 BAM 文件。
-   `*.mt.no.softclip.bam`: 移除了软剪切（soft-clipped）reads 的线粒体 BAM 文件，用于下游分析。
-   `*.conrmdup.v5.end10.txt`: **最终的突变结果文件**。这是一个文本文件，包含了检测到的低频突变位点、参考碱基、突变碱基、频率、覆盖深度等信息。
-   其他中间文件：流程会生成大量中间文件（如 `.raw.bam`, `.sort.picard.bam` 等），用于调试和追溯。

## 5. 持续更新与贡献
本流程目前处于持续开发和优化阶段。我们欢迎任何形式的反馈、问题报告或功能建议。

如果您在使用中遇到问题，或者有改进建议，请在 GitHub 仓库中提交一个 [Issue](https://github.com/your-username/MitoVariantPipe/issues)。

## 6. 许可证 (License)
本项目采用 [MIT License](./LICENSE)。
