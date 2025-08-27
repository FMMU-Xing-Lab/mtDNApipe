#!/bin/bash
set -e # 如果任何命令失败，则立即退出

# --- Default parameters ---
THREADS=8
MIN_FREQ=0.1
MIN_READS=5

# --- Help Message ---
usage() {
    echo "Usage: $0 -i <sample_id> -1 <read1.fq.gz> -2 <read2.fq.gz> -r <ref.fa> -o <output_dir>"
    echo "Options:"
    echo "  -i  Sample ID (e.g., sample01)"
    echo "  -1  Path to Read 1 FASTQ file"
    echo "  -2  Path to Read 2 FASTQ file"
    echo "  -r  Path to the mitochondrial reference genome (e.g., hg19_mt.fa)"
    echo "  -o  Output directory"
    echo "  -t  Number of threads (default: $THREADS)"
    echo "  -h  Display this help message"
    exit 1
}

# --- Parse command-line arguments ---
while getopts "i:1:2:r:o:t:h" opt; do
    case ${opt} in
        i) SAMPLE_ID=$OPTARG ;;
        1) R1_FQ=$OPTARG ;;
        2) R2_FQ=$OPTARG ;;
        r) REF_GENOME=$OPTARG ;;
        o) OUTPUT_DIR=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        h) usage ;;
        \?) usage ;;
    esac
done

# --- Check for mandatory arguments ---
if [ -z "$SAMPLE_ID" ] || [ -z "$R1_FQ" ] || [ -z "$R2_FQ" ] || [ -z "$REF_GENOME" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing mandatory arguments."
    usage
fi

# --- Get base directory of the pipeline ---
PIPELINE_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
SCRIPTS_DIR="$PIPELINE_DIR/scripts"

# --- Setup Paths, Status Directory and Log File ---
SAMPLE_OUT_DIR="$OUTPUT_DIR/$SAMPLE_ID"
STATUS_DIR="$SAMPLE_OUT_DIR/status_flags" # 用于存放成功标记文件的目录
LOG_FILE="$SAMPLE_OUT_DIR/${SAMPLE_ID}_analysis.log"

# --- Create necessary directories ---
mkdir -p "$SAMPLE_OUT_DIR"
mkdir -p "$STATUS_DIR" # 创建标记文件目录

# --- Start Logging ---
# 只有在日志文件不存在时才添加起始标题，避免重复运行时多次写入
if [ ! -f "$LOG_FILE" ]; then
    echo "--- Starting Pipeline for Sample: $SAMPLE_ID ---" | tee -a "$LOG_FILE"
    date | tee -a "$LOG_FILE"
fi
echo "---" | tee -a "$LOG_FILE"
echo "Pipeline resumed or started for Sample: $SAMPLE_ID" | tee -a "$LOG_FILE"
date | tee -a "$LOG_FILE"
echo "---" | tee -a "$LOG_FILE"

# --- Step 1: Custom Python Filtering ---
STEP1_FLAG="$STATUS_DIR/01_filter_fastq.success"
if [ -f "$STEP1_FLAG" ]; then
    echo "[1/9] Filtering raw FASTQ files... SKIPPED (already done)" | tee -a "$LOG_FILE"
else
    echo "[1/9] Filtering raw FASTQ files..." | tee -a "$LOG_FILE"
    python3 $SCRIPTS_DIR/filter_fastq.py \
        --r1 "$R1_FQ" --r2 "$R2_FQ" \
        --o1 "$SAMPLE_OUT_DIR/${SAMPLE_ID}.noReads_R1.fq.gz" \
        --o2 "$SAMPLE_OUT_DIR/${SAMPLE_ID}.noReads_R2.fq.gz" \
        --log "$SAMPLE_OUT_DIR/${SAMPLE_ID}_noReads.log" \
        --window 10 --threshold 30
    touch "$STEP1_FLAG" # 成功后创建标记文件
fi

# --- Step 2: fastp Trimming ---
STEP2_FLAG="$STATUS_DIR/02_fastp_trim.success"
if [ -f "$STEP2_FLAG" ]; then
    echo "[2/9] Trimming with fastp... SKIPPED (already done)" | tee -a "$LOG_FILE"
else
    echo "[2/9] Trimming with fastp..." | tee -a "$LOG_FILE"
    fastp -i "$SAMPLE_OUT_DIR/${SAMPLE_ID}.noReads_R1.fq.gz" \
        -o "$SAMPLE_OUT_DIR/${SAMPLE_ID}.fastp_R1.good.fq.gz" \
        -I "$SAMPLE_OUT_DIR/${SAMPLE_ID}.noReads_R2.fq.gz" \
        -O "$SAMPLE_OUT_DIR/${SAMPLE_ID}.fastp_R2.good.fq.gz" \
        -h "$SAMPLE_OUT_DIR/${SAMPLE_ID}_fastp.html" \
        --detect_adapter_for_pe --length_required 50 -e 30 &>> "$LOG_FILE"
    touch "$STEP2_FLAG"
fi

# # --- Step 3: C++ Overlap Correction ---
STEP3_FLAG="$STATUS_DIR/03_Overlap_correction.success"
if [ -f "$STEP3_FLAG" ]; then
    echo "[3/9]  Correcting overlaps... SKIPPED (already done)" | tee -a "$LOG_FILE"
else
    echo "[3/9]  Correcting overlaps..." | tee -a "$LOG_FILE"
    gunzip  "$SAMPLE_OUT_DIR/${SAMPLE_ID}.fastp_R1.good.fq.gz"
    gunzip  "$SAMPLE_OUT_DIR/${SAMPLE_ID}.fastp_R2.good.fq.gz"
    $SCRIPTS_DIR/overlap_corrector/overlap_corrector \
    "$SAMPLE_OUT_DIR/${SAMPLE_ID}.fastp_R1.good.fq" \
    "$SAMPLE_OUT_DIR/${SAMPLE_ID}.fastp_R2.good.fq" \
    "$SAMPLE_OUT_DIR/${SAMPLE_ID}_R1.good.fq" \
    "$SAMPLE_OUT_DIR/${SAMPLE_ID}_R2.good.fq"

    gzip "$SAMPLE_OUT_DIR/${SAMPLE_ID}.fastp_R1.good.fq"
    gzip "$SAMPLE_OUT_DIR/${SAMPLE_ID}.fastp_R2.good.fq"
    gzip "$SAMPLE_OUT_DIR/${SAMPLE_ID}_R1.good.fq"
    gzip "$SAMPLE_OUT_DIR/${SAMPLE_ID}_R2.good.fq"

    touch "$STEP3_FLAG"
fi

# --- Step 4: BWA Alignment ---
STEP4_FLAG="$STATUS_DIR/04_BWA_alignment.success"
if [ -f "$STEP4_FLAG" ]; then
    echo "[4/9] Aligning with BWA-MEM... SKIPPED (already done)" | tee -a "$LOG_FILE"
else
    echo "[4/9] Aligning with BWA-MEM..." | tee -a "$LOG_FILE"
    # 使用正确的BWA路径和RG字符串（使用$'\t'转义tab）
    bwa mem -t "$THREADS" -M -R "@RG\tID:${SAMPLE_ID}\tLB:mtDNA\tPL:ILLUMINA\tSM:${SAMPLE_ID}" \
        "$REF_GENOME" \
        "$SAMPLE_OUT_DIR/${SAMPLE_ID}_R1.good.fq.gz" \
        "$SAMPLE_OUT_DIR/${SAMPLE_ID}_R2.good.fq.gz" 2>> "$LOG_FILE" | \
        samtools view -bS - -o "$SAMPLE_OUT_DIR/${SAMPLE_ID}.raw.bam"
    echo "Step 4 finished ~" | tee -a "$LOG_FILE"  # 可选：添加日志
    touch "$STEP4_FLAG"
fi

# --- Step 5: Picard & GATK Pre-processing ---
STEP5_FLAG="$STATUS_DIR/05_Picard_Pre-processing.success"
if [ -f "$STEP5_FLAG" ]; then
    echo "[5/9] Sorting, Indexing, and Realigning BAM... SKIPPED (already done)" | tee -a "$LOG_FILE"
else
    echo "[5/9] Sorting, Indexing, and Realigning BAM..." | tee -a "$LOG_FILE"
    PICARD_TMP_DIR="$SAMPLE_OUT_DIR/tmp"
    mkdir -p "$PICARD_TMP_DIR"
    java -Xmx10g -Djava.io.tmpdir="$PICARD_TMP_DIR" -jar ${which picard} SortSam.jar \
        I="$SAMPLE_OUT_DIR/${SAMPLE_ID}.raw.bam" O="$SAMPLE_OUT_DIR/${SAMPLE_ID}.sort.picard.bam" \
        SO=coordinate VALIDATION_STRINGENCY=SILENT 2>> "$LOG_FILE"

    java -Xmx10g -Djava.io.tmpdir="$PICARD_TMP_DIR" -jar ${which picard} BuildBamIndex.jar \
        INPUT="$SAMPLE_OUT_DIR/${SAMPLE_ID}.sort.picard.bam" VALIDATION_STRINGENCY=SILENT 2>> "$LOG_FILE"

    java -Xmx10g -Djava.io.tmpdir="$PICARD_TMP_DIR" -jar ${which GATK} GenomeAnalysisTK.jar \
        -R "$REF_GENOME" -T RealignerTargetCreator \
        -o "$SAMPLE_OUT_DIR/${SAMPLE_ID}.realn.intervals" \
        -I "$SAMPLE_OUT_DIR/${SAMPLE_ID}.sort.picard.bam" 2>> "$LOG_FILE"

    java -Xmx10g -Djava.io.tmpdir="$PICARD_TMP_DIR" -jar ${which GATK} GenomeAnalysisTK.jar \
        -R "$REF_GENOME" -T IndelRealigner \
        --maxReadsForRealignment 100000 \
        -targetIntervals "$SAMPLE_OUT_DIR/${SAMPLE_ID}.realn.intervals" \
        -o "$SAMPLE_OUT_DIR/${SAMPLE_ID}.realn.bam" \
        -I "$SAMPLE_OUT_DIR/${SAMPLE_ID}.sort.picard.bam" 2>> "$LOG_FILE"
    rm -rf "$PICARD_TMP_DIR" # Clean up temp dir

    touch "$STEP5_FLAG"
fi

# --- Step 6: Mitochondria-specific steps ---
STEP6_FLAG="$STATUS_DIR/06_Mitochondria-specific.success"
if [ -f "$STEP6_FLAG" ]; then
    echo "[6/9] Filtering for mitochondrial reads and removing soft-clips... SKIPPED (already done)" | tee -a "$LOG_FILE"
else
    echo "[6/9] Filtering for mitochondrial reads and removing soft-clips..." | tee -a "$LOG_FILE"
    samtools index "$SAMPLE_OUT_DIR/${SAMPLE_ID}.realn.bam"
    perl $SCRIPTS_DIR/fetch_mt.pl "$SAMPLE_OUT_DIR/${SAMPLE_ID}.realn.bam" | \
        samtools view -bS - -o "$SAMPLE_OUT_DIR/${SAMPLE_ID}.mt.bam" 2>> "$LOG_FILE"
    samtools index "$SAMPLE_OUT_DIR/${SAMPLE_ID}.mt.bam"

    python3 $SCRIPTS_DIR/To_softclip_bam.py \
        "$SAMPLE_OUT_DIR/${SAMPLE_ID}.mt.bam" \
        "$SAMPLE_OUT_DIR/${SAMPLE_ID}.mt.softclip.bam" \
        "$SAMPLE_OUT_DIR/${SAMPLE_ID}.mt.no.softclip.bam"
    samtools index "$SAMPLE_OUT_DIR/${SAMPLE_ID}.mt.no.softclip.bam"
    touch "$STEP6_FLAG"
fi

# --- Step 7: Consensus-based Deduplication ---
STEP7_FLAG="$STATUS_DIR/07_Consensus-based_Deduplication.success"
if [ -f "$STEP7_FLAG" ]; then
    echo "[7/9] Performing consensus-based deduplication... SKIPPED (already done)" | tee -a "$LOG_FILE"
else
    echo "[7/9] Performing consensus-based deduplication..." | tee -a "$LOG_FILE"
    python3 $SCRIPTS_DIR/con_rmdup_v1.5.py \
        "$SAMPLE_ID" \
        "$SAMPLE_OUT_DIR/${SAMPLE_ID}.mt.no.softclip.bam" \
        "$OUTPUT_DIR" 30 4
    touch "$STEP7_FLAG"
fi

# --- Step 8: Convert Binary to Text ---
STEP8_FLAG="$STATUS_DIR/08_convert_binary.success"
if [ -f "$STEP8_FLAG" ]; then
    echo "[8/9] Converting binary genotype files to text... SKIPPED (already done)" | tee -a "$LOG_FILE"
else
    echo "[8/9] Converting binary genotype files to text..." | tee -a "$LOG_FILE"
    python3 $SCRIPTS_DIR/convert_binary_to_text_v1.py \
        --input_bin "$SAMPLE_OUT_DIR/${SAMPLE_ID}.RawGenotypes.Total" \
        --output_txt "$SAMPLE_OUT_DIR/${SAMPLE_ID}.RawGenotypes.Total.txt" \
        --ref_file "$PIPELINE_DIR/ref/chrM_refAllele.txt" \
        --file_type Total

    python3 $SCRIPTS_DIR/convert_binary_to_text_v1.py \
        --input_bin "$SAMPLE_OUT_DIR/${SAMPLE_ID}.RawGenotypes.cordup" \
        --output_txt "$SAMPLE_OUT_DIR/${SAMPLE_ID}.RawGenotypes.cordup.txt" \
        --ref_file "$PIPELINE_DIR/ref/chrM_refAllele.txt" \
        --file_type cordup
    touch "$STEP8_FLAG"
fi

# --- Step 9: Call Mutations ---
STEP9_FLAG="$STATUS_DIR/09_call_mutations.success"
if [ -f "$STEP9_FLAG" ]; then
    echo "[9/9] Calling low-frequency mutations... SKIPPED (already done)" | tee -a "$LOG_FILE"
else
    echo "[9/9] Calling low-frequency mutations..." | tee -a "$LOG_FILE"
    python3 $SCRIPTS_DIR/call_mut_v1.6.py \
        "$SAMPLE_ID" "$OUTPUT_DIR" \
        "$SAMPLE_OUT_DIR/${SAMPLE_ID}.conrmdup.v5.end10.txt" 10 \
        --min_freq $MIN_FREQ --min_reads $MIN_READS
    touch "$STEP9_FLAG"
fi

echo "--- Pipeline Finished Successfully for Sample: $SAMPLE_ID ---" | tee -a "$LOG_FILE"
date | tee -a "$LOG_FILE"
