import pysam
import os
import sys
import numpy as np
import pandas as pd
from collections import defaultdict
import gc
import multiprocessing as mp
from multiprocessing import Pool
import struct
import array
import time
import ctypes
from tqdm import tqdm

# 定义DNA碱基
dna_letters = ['A', 'C', 'G', 'T', 'N']

# 获取命令行参数
sample = sys.argv[1]
bam_file = sys.argv[2]
out_dir = sys.argv[3]
BaseQ_thld_hi = int(sys.argv[4])
num_processes = int(sys.argv[5]) if len(sys.argv) > 5 else max(1, mp.cpu_count() - 2)

# 创建输出目录
os.makedirs(f"{out_dir}/{sample}", exist_ok=True)

# 输出文件路径
out_genotypeTotal_file = f"{out_dir}/{sample}/{sample}.RawGenotypes.Total"
out_genotypeVerySensitive_file = f"{out_dir}/{sample}/{sample}.RawGenotypes.cordup"
out_totalCts_file = f"{out_dir}/{sample}/{sample}.QualifiedTotalCts"
out_family_size_dist = f"{out_dir}/{sample}/{sample}.FamilySizeDistribution.txt"
out_molecule_length_dist = f"{out_dir}/{sample}/{sample}.MoleculeLengthDistribution.txt"
mito_ref_file = r"/mnt/data2/lishengjing/20250509_con_rmdup/chrM_refAllele.txt"

# 创建临时目录存储中间结果
temp_dir = f"{out_dir}/{sample}/temp"
os.makedirs(temp_dir, exist_ok=True)

# 读取参考基因组
mito_ref = pd.read_table(mito_ref_file, names=["pos", "base"])
max_bp = mito_ref.shape[0]


# 轻量级Read表示
class LightRead:
    __slots__ = ('seq', 'qual', 'is_reverse', 'aligned_pairs', 'query_name')

    def __init__(self, read):
        self.seq = read.seq
        self.qual = read.query_qualities
        self.is_reverse = read.is_reverse
        self.aligned_pairs = np.array(read.get_aligned_pairs(matches_only=True), dtype=np.int32)
        self.query_name = read.query_name


# 共享内存数据结构
shared_family_size = mp.Array(ctypes.c_uint, 101)  # 0-100
shared_mol_length = mp.Array(ctypes.c_uint, 1001)  # 0-1000


# ================== 核心处理函数 ==================

def process_molecule_group(args):
    """处理一组分子的核心函数，优化内存和速度"""
    group_idx, molecule_data_list, BaseQ_thld_hi, max_bp, temp_dir = args

    # 初始化结果矩阵
    TotalMoleculeCtsMatrix = np.zeros((max_bp, 2), dtype=np.uint32)

    # 初始化统计数组（本地）
    family_size_counts = array.array('I', [0] * 101)  # 支持最大family size为100
    molecule_length_counts = array.array('I', [0] * 1001)  # 支持最大长度1000bp

    # 创建临时文件
    with open(f"{temp_dir}/total_{group_idx}.bin", "wb") as total_out, \
            open(f"{temp_dir}/cordup_{group_idx}.bin", "wb") as cordup_out:

        # 处理每个分子
        for molecule_id, read_pairs in molecule_data_list:
            # 初始化基因型矩阵
            SG_Genotypes = np.zeros((max_bp, 5), dtype=np.uint8)
            DB_Genotypes = np.zeros((max_bp, 5), dtype=np.uint8)
            Strand_mtx = np.zeros((max_bp, 2), dtype=np.uint8)

            # 统计分子长度
            start_pos, end_pos = map(int, molecule_id.split('_'))
            molecule_length = abs(end_pos - start_pos)
            if molecule_length < 1000:
                molecule_length_counts[molecule_length] += 1
            else:
                molecule_length_counts[1000] += 1

            # 处理每个read对
            for read1, read2 in read_pairs:
                # 获取序列和质量信息
                seq1, seq2 = read1.seq, read2.seq
                qual1, qual2 = read1.qual, read2.qual

                # 获取对齐位置
                pos_array1 = read1.aligned_pairs
                pos_array2 = read2.aligned_pairs

                # 找出重叠和非重叠区域
                if len(pos_array1) > 0 and len(pos_array2) > 0:
                    ref_positions1 = pos_array1[:, 1]
                    ref_positions2 = pos_array2[:, 1]
                    overlap_positions = np.intersect1d(ref_positions1, ref_positions2)
                else:
                    overlap_positions = np.array([])

                # 处理非重叠区域 (read1特有)
                if len(pos_array1) > 0:
                    read1_unique_pos = pos_array1[~np.isin(pos_array1[:, 1], overlap_positions)]
                    for row in read1_unique_pos:
                        read_pos, ref_pos = row
                        if qual1[read_pos] > BaseQ_thld_hi:
                            base = seq1[read_pos]
                            if base in dna_letters:
                                base_idx = dna_letters.index(base)
                                SG_Genotypes[ref_pos, base_idx] = min(SG_Genotypes[ref_pos, base_idx] + 1, 255)
                                Strand_mtx[ref_pos, int(read1.is_reverse)] = min(
                                    Strand_mtx[ref_pos, int(read1.is_reverse)] + 1, 255)
                        else:
                            SG_Genotypes[ref_pos, 4] = min(SG_Genotypes[ref_pos, 4] + 1, 255)

                # 处理非重叠区域 (read2特有)
                if len(pos_array2) > 0:
                    read2_unique_pos = pos_array2[~np.isin(pos_array2[:, 1], overlap_positions)]
                    for row in read2_unique_pos:
                        read_pos, ref_pos = row
                        if qual2[read_pos] > BaseQ_thld_hi:
                            base = seq2[read_pos]
                            if base in dna_letters:
                                base_idx = dna_letters.index(base)
                                SG_Genotypes[ref_pos, base_idx] = min(SG_Genotypes[ref_pos, base_idx] + 1, 255)
                                Strand_mtx[ref_pos, int(read2.is_reverse)] = min(
                                    Strand_mtx[ref_pos, int(read2.is_reverse)] + 1, 255)
                        else:
                            SG_Genotypes[ref_pos, 4] = min(SG_Genotypes[ref_pos, 4] + 1, 255)

                # 处理重叠区域
                if len(overlap_positions) > 0:
                    # 获取重叠区域的位置
                    read1_overlap = pos_array1[np.isin(pos_array1[:, 1], overlap_positions)]
                    read2_overlap = pos_array2[np.isin(pos_array2[:, 1], overlap_positions)]

                    # 按参考位置排序
                    read1_overlap = read1_overlap[read1_overlap[:, 1].argsort()]
                    read2_overlap = read2_overlap[read2_overlap[:, 1].argsort()]

                    # 确保长度一致
                    min_len = min(len(read1_overlap), len(read2_overlap))
                    if min_len > 0:
                        for i in range(min_len):
                            read1_pos, ref_pos = read1_overlap[i]
                            read2_pos, _ = read2_overlap[i]

                            base1, base2 = seq1[read1_pos], seq2[read2_pos]
                            qual1_val, qual2_val = qual1[read1_pos], qual2[read2_pos]

                            # 双链一致
                            if base1 == base2:
                                if qual1_val > BaseQ_thld_hi or qual2_val > BaseQ_thld_hi:
                                    if base1 in dna_letters:
                                        base_idx = dna_letters.index(base1)
                                        DB_Genotypes[ref_pos, base_idx] = min(DB_Genotypes[ref_pos, base_idx] + 1, 255)
                                        Strand_mtx[ref_pos, 0] = min(Strand_mtx[ref_pos, 0] + 1, 255)
                                        Strand_mtx[ref_pos, 1] = min(Strand_mtx[ref_pos, 1] + 1, 255)
                                else:
                                    DB_Genotypes[ref_pos, 4] = min(DB_Genotypes[ref_pos, 4] + 1, 255)
                            # 双链不一致
                            else:
                                if qual1_val > BaseQ_thld_hi and qual2_val > BaseQ_thld_hi:
                                    DB_Genotypes[ref_pos, 4] = min(DB_Genotypes[ref_pos, 4] + 1, 255)
                                elif qual1_val > BaseQ_thld_hi:
                                    if base1 in dna_letters:
                                        base_idx = dna_letters.index(base1)
                                        DB_Genotypes[ref_pos, base_idx] = min(DB_Genotypes[ref_pos, base_idx] + 1, 255)
                                        Strand_mtx[ref_pos, 0] = min(Strand_mtx[ref_pos, 0] + 1, 255)
                                        Strand_mtx[ref_pos, 1] = min(Strand_mtx[ref_pos, 1] + 1, 255)
                                elif qual2_val > BaseQ_thld_hi:
                                    if base2 in dna_letters:
                                        base_idx = dna_letters.index(base2)
                                        DB_Genotypes[ref_pos, base_idx] = min(DB_Genotypes[ref_pos, base_idx] + 1, 255)
                                        Strand_mtx[ref_pos, 0] = min(Strand_mtx[ref_pos, 0] + 1, 255)
                                        Strand_mtx[ref_pos, 1] = min(Strand_mtx[ref_pos, 1] + 1, 255)
                                else:
                                    DB_Genotypes[ref_pos, 4] = min(DB_Genotypes[ref_pos, 4] + 1, 255)

            # 识别突变
            covered_positions = np.where(np.sum(SG_Genotypes + DB_Genotypes, axis=1) > 0)[0]
            for pos_idx in covered_positions:
                total_genotype = SG_Genotypes[pos_idx] + DB_Genotypes[pos_idx]
                FamSize = sum(total_genotype[:4])

                # 更新family size统计
                if FamSize < 100:
                    family_size_counts[FamSize] += 1
                else:
                    family_size_counts[100] += 1

                if FamSize > 0:
                    CallIdx = np.argmax(total_genotype[:4])
                    Call = dna_letters[CallIdx]
                    Ref = mito_ref["base"].iloc[pos_idx].upper()

                    # 位置校正
                    pos_raw = pos_idx + 1
                    pos_cor = pos_raw - 16569 if pos_raw > 16569 else pos_raw

                    # 计算统计量
                    GT_Cts = total_genotype[CallIdx]
                    SG_Cts = SG_Genotypes[pos_idx, CallIdx]
                    DB_Cts = DB_Genotypes[pos_idx, CallIdx]
                    CSS = GT_Cts / FamSize if FamSize > 0 else 0

                    # 判断是否记录到严格输出
                    record_to_strict = False
                    is_overlap = DB_Genotypes[pos_idx].sum() > 0

                    # 应用新规则
                    if not is_overlap:
                        if FamSize == 1:
                            record_to_strict = True
                        elif FamSize == 2:
                            if GT_Cts == 2:
                                record_to_strict = True
                            else:
                                Call = 'N'
                        else:
                            if CSS > 0.66:
                                record_to_strict = True
                    else:
                        if FamSize == 1:
                            record_to_strict = True
                        elif FamSize == 2:
                            if GT_Cts == 2:
                                record_to_strict = True
                            else:
                                Call = 'N'
                        else:
                            if CSS > 0.66:
                                record_to_strict = True

                    # 更新总覆盖计数
                    TotalMoleculeCtsMatrix[pos_idx, 0] += 1
                    if record_to_strict:
                        TotalMoleculeCtsMatrix[pos_idx, 1] += 1

                    # 写入突变记录（仅当是突变时）
                    if Call != Ref and Call != 'N':
                        # 获取分子起始结束位置
                        start, end = map(int, molecule_id.split('_'))

                        # 获取链支持信息
                        forward = int(Strand_mtx[pos_idx, 0])
                        reverse = int(Strand_mtx[pos_idx, 1])

                        # 获取碱基索引
                        base_idx = dna_letters.index(Call)

                        # 修复：使用正确的字段数和格式字符串
                        variant_data = struct.pack('IIIIdIIIIII',
                                                   start,
                                                   end,
                                                   pos_cor,
                                                   base_idx,
                                                   CSS,
                                                   FamSize,
                                                   GT_Cts,
                                                   DB_Cts,
                                                   SG_Cts,
                                                   forward,
                                                   reverse)

                        # 写入到总文件
                        total_out.write(variant_data)

                        # 如果符合严格标准，写入到严格文件
                        if record_to_strict:
                            cordup_out.write(variant_data)

    # 保存统计结果
    np.save(f"{temp_dir}/cts_matrix_{group_idx}.npy", TotalMoleculeCtsMatrix)
    with open(f"{temp_dir}/family_size_{group_idx}.bin", "wb") as f:
        family_size_counts.tofile(f)
    with open(f"{temp_dir}/mol_length_{group_idx}.bin", "wb") as f:
        molecule_length_counts.tofile(f)

    return group_idx


# ================== 主函数 ==================

def main():
    print(f"Starting processing for sample {sample} with {num_processes} processes...")
    start_time = time.time()

    # 步骤1: 构建分子字典（使用轻量级Read对象）
    print("Building molecule dictionary...")
    bam_input = pysam.AlignmentFile(bam_file, "rb")
    MoleculeDict = {}
    read_cache = {}
    read_count = 0

    for read in bam_input:
        read_count += 1
        if read_count % 1000000 == 0:
            print(f"Processed {read_count} reads")

        # 缓存未配对的reads
        if read.query_name not in read_cache:
            read_cache[read.query_name] = LightRead(read)
            continue

        # 配对reads
        read1 = read_cache.pop(read.query_name)
        read2 = LightRead(read)

        # 识别正向和反向reads
        if read1.is_reverse and not read2.is_reverse:
            fwd_read, rev_read = read2, read1
        elif not read1.is_reverse and read2.is_reverse:
            fwd_read, rev_read = read1, read2
        else:
            continue

        # 创建分子ID
        if fwd_read.aligned_pairs.size > 0:
            start_pos = fwd_read.aligned_pairs[0, 1]
        else:
            continue

        if rev_read.aligned_pairs.size > 0:
            end_pos = rev_read.aligned_pairs[-1, 1]
        else:
            continue

        molecule_id = f"{start_pos}_{end_pos}"

        # 存储分子信息
        if molecule_id not in MoleculeDict:
            MoleculeDict[molecule_id] = []
        MoleculeDict[molecule_id].append((fwd_read, rev_read))

    bam_input.close()
    print(f"Identified {len(MoleculeDict)} molecules from {read_count} reads")
    print(f"Molecule dictionary built in {time.time() - start_time:.2f} seconds")

    # 步骤2: 将分子分成组用于并行处理
    print("Preparing data for parallel processing...")
    molecules_grouped = []
    group_size = max(1000, len(MoleculeDict) // (num_processes * 10))
    current_group = []
    group_idx = 0

    # 转换为列表并删除原始字典以释放内存
    molecule_list = list(MoleculeDict.items())
    del MoleculeDict
    gc.collect()

    for i, (mol_id, read_pairs) in enumerate(molecule_list):
        current_group.append((mol_id, read_pairs))

        if len(current_group) >= group_size:
            molecules_grouped.append((group_idx, current_group, BaseQ_thld_hi, max_bp, temp_dir))
            current_group = []
            group_idx += 1

    if current_group:
        molecules_grouped.append((group_idx, current_group, BaseQ_thld_hi, max_bp, temp_dir))
        group_idx += 1

    print(f"Split into {len(molecules_grouped)} groups for parallel processing")

    # 步骤3: 并行处理分子组
    print("Starting parallel processing...")
    process_start = time.time()

    with Pool(processes=num_processes) as pool:
        results = list(tqdm(pool.imap(process_molecule_group, molecules_grouped),
                            total=len(molecules_grouped),
                            desc="Processing groups"))

    print(f"Parallel processing completed in {time.time() - process_start:.2f} seconds")

    # 步骤4: 合并结果
    print("Merging results...")
    merge_start = time.time()

    # 初始化合并矩阵
    TotalMoleculeCtsMatrix_combined = np.zeros((max_bp, 2), dtype=np.uint32)
    family_size_dist_combined = array.array('I', [0] * 101)
    molecule_length_dist_combined = array.array('I', [0] * 1001)

    # 合并覆盖矩阵
    for idx in results:
        cts_matrix = np.load(f"{temp_dir}/cts_matrix_{idx}.npy")
        TotalMoleculeCtsMatrix_combined += cts_matrix

        # 合并family size分布
        with open(f"{temp_dir}/family_size_{idx}.bin", "rb") as f:
            family_size_chunk = array.array('I')
            family_size_chunk.fromfile(f, 101)
            for i in range(101):
                family_size_dist_combined[i] += family_size_chunk[i]

        # 合并分子长度分布
        with open(f"{temp_dir}/mol_length_{idx}.bin", "rb") as f:
            mol_length_chunk = array.array('I')
            mol_length_chunk.fromfile(f, 1001)
            for i in range(1001):
                molecule_length_dist_combined[i] += mol_length_chunk[i]

    # 环状基因组处理
    if max_bp == 17119:
        for i in range(550):
            for col in range(2):
                TotalMoleculeCtsMatrix_combined[i][col] += TotalMoleculeCtsMatrix_combined[16569 + i][col]
        TotalMoleculeCtsMatrix_combined = TotalMoleculeCtsMatrix_combined[:16569]

    # 合并突变文件
    with open(out_genotypeTotal_file, "wb") as total_out, \
            open(out_genotypeVerySensitive_file, "wb") as cordup_out:

        for idx in results:
            # 合并Total文件
            with open(f"{temp_dir}/total_{idx}.bin", "rb") as f:
                while True:
                    # 使用正确的记录大小
                    record_size = struct.calcsize('IIIIdIIIIII')
                    data = f.read(record_size)
                    if not data:
                        break
                    total_out.write(data)

            # 合并cordup文件
            with open(f"{temp_dir}/cordup_{idx}.bin", "rb") as f:
                while True:
                    record_size = struct.calcsize('IIIIdIIIIII')
                    data = f.read(record_size)
                    if not data:
                        break
                    cordup_out.write(data)

    # 步骤5: 输出最终结果
    # 输出覆盖深度文件
    with open(out_totalCts_file, "w") as out_totalCts:
        for pos in range(len(TotalMoleculeCtsMatrix_combined)):
            out_totalCts.write(
                f"{pos + 1}\t{TotalMoleculeCtsMatrix_combined[pos][0]}\t{TotalMoleculeCtsMatrix_combined[pos][1]}\n")

    # 输出family size分布
    with open(out_family_size_dist, "w") as f_out:
        f_out.write("FamilySize\tCount\n")
        for size, count in enumerate(family_size_dist_combined):
            if count > 0:
                f_out.write(f"{size}\t{count}\n")

    # 输出molecule长度分布
    with open(out_molecule_length_dist, "w") as f_out:
        f_out.write("MoleculeLength\tCount\n")
        for length, count in enumerate(molecule_length_dist_combined):
            if count > 0 and length < 1000:
                f_out.write(f"{length}\t{count}\n")
        if molecule_length_dist_combined[1000] > 0:
            f_out.write(f"1000+\t{molecule_length_dist_combined[1000]}\n")

    # 清理临时文件
    for idx in results:
        for suffix in [f"cts_matrix_{idx}.npy",
                       f"family_size_{idx}.bin",
                       f"mol_length_{idx}.bin",
                       f"total_{idx}.bin",
                       f"cordup_{idx}.bin"]:
            try:
                os.remove(f"{temp_dir}/{suffix}")
            except:
                pass

    print(f"Processing completed in {time.time() - start_time:.2f} seconds")


if __name__ == "__main__":
    main()
