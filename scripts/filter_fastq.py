#!/usr/bin/env python3
import gzip
import argparse
from itertools import zip_longest
import sys

def check_window_quality(qual_str, window_size, threshold):
    """
    检查质量字符串是否在任一滑窗中平均质量低于阈值
    :param qual_str: 质量字符串（ASCII字符）
    :param window_size: 滑窗大小
    :param threshold: 质量阈值（整数值）
    :return: 布尔值，True表示质量合格，False表示不合格
    """
    n = len(qual_str)
    # 如果序列长度小于窗口大小，则检查整个序列
    if n <= window_size:
        # 计算整个序列的平均质量
        total_qual = sum(ord(q) - 33 for q in qual_str)
        if n == 0:
            return True  # 避免除零错误
        average_qual = total_qual / n
        return average_qual >= threshold
    
    # 计算所有质量值的列表
    quals = [ord(q) - 33 for q in qual_str]
    # 使用累加和快速计算窗口平均质量
    cumsum = [0]
    for q_val in quals:
        cumsum.append(cumsum[-1] + q_val)
    
    # 滑动窗口检查
    for i in range(len(qual_str) - window_size + 1):
        start = i
        end = i + window_size
        total = cumsum[end] - cumsum[start]
        avg_qual = total / window_size
        if avg_qual < threshold:
            return False
    return True

def main():
    parser = argparse.ArgumentParser(description='Filter paired-end reads by sliding window quality')
    parser.add_argument('--r1', required=True, help='Input R1 FASTQ gzipped file')
    parser.add_argument('--r2', required=True, help='Input R2 FASTQ gzipped file')
    parser.add_argument('--o1', required=True, help='Output filtered R1 FASTQ gzipped file')
    parser.add_argument('--o2', required=True, help='Output filtered R2 FASTQ gzipped file')
    parser.add_argument('--log', required=True, help='Log file path')
    parser.add_argument('--window', type=int, required=True, help='Sliding window size')
    parser.add_argument('--threshold', type=int, required=True, help='Quality score threshold (Phred+33)')
    
    args = parser.parse_args()
    
    # 初始化计数器
    r1_total = r2_total = 0
    r1_self_fail = r2_self_fail = 0
    r1_pair_fail = r2_pair_fail = 0
    
    # 打开输入输出文件
    try:
        with gzip.open(args.r1, 'rt') as r1_in, \
             gzip.open(args.r2, 'rt') as r2_in, \
             gzip.open(args.o1, 'wt') as r1_out, \
             gzip.open(args.o2, 'wt') as r2_out:
            
            while True:
                # 读取R1的四行记录
                r1_header = r1_in.readline().rstrip()
                if not r1_header: 
                    break
                r1_seq = r1_in.readline().rstrip()
                r1_plus = r1_in.readline().rstrip()  # 通常为"+"
                r1_qual = r1_in.readline().rstrip()
                r1_total += 1
                
                # 读取R2的四行记录
                r2_header = r2_in.readline().rstrip()
                r2_seq = r2_in.readline().rstrip()
                r2_plus = r2_in.readline().rstrip()
                r2_qual = r2_in.readline().rstrip()
                r2_total += 1
                
                # 检查R1和R2的质量
                r1_ok = check_window_quality(r1_qual, args.window, args.threshold)
                r2_ok = check_window_quality(r2_qual, args.window, args.threshold)
                
                # 根据质量决定是否保留
                if r1_ok and r2_ok:
                    # 双方都合格，写入输出文件
                    r1_out.write(f"{r1_header}\n{r1_seq}\n{r1_plus}\n{r1_qual}\n")
                    r2_out.write(f"{r2_header}\n{r2_seq}\n{r2_plus}\n{r2_qual}\n")
                else:
                    # 记录失败原因
                    if not r1_ok:
                        r1_self_fail += 1
                    else:
                        r1_pair_fail += 1
                    
                    if not r2_ok:
                        r2_self_fail += 1
                    else:
                        r2_pair_fail += 1
                        
    except Exception as e:
        sys.stderr.write(f"Error processing files: {str(e)}\n")
        sys.exit(1)
    
    # 计算最终保留的reads数
    r1_kept = r1_total - (r1_self_fail + r1_pair_fail)
    r2_kept = r2_total - (r2_self_fail + r2_pair_fail)
    
    # 写入日志文件
    with open(args.log, 'w') as log_file:
        log_file.write(f"Total R1 reads: {r1_total}\n")
        log_file.write(f"Total R2 reads: {r2_total}\n")
        log_file.write(f"R1 failed due to self low quality: {r1_self_fail}\n")
        log_file.write(f"R1 failed due to pair's low quality: {r1_pair_fail}\n")
        log_file.write(f"R2 failed due to self low quality: {r2_self_fail}\n")
        log_file.write(f"R2 failed due to pair's low quality: {r2_pair_fail}\n")
        log_file.write(f"Final kept R1 reads: {r1_kept}\n")
        log_file.write(f"Final kept R2 reads: {r2_kept}\n")

if __name__ == '__main__':
    main()
