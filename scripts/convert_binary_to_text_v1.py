import sys
import struct
import os
import pandas as pd
from tqdm import tqdm
import argparse

# 定义DNA碱基
dna_letters = ['A', 'C', 'G', 'T']

def main():
    parser = argparse.ArgumentParser(description="Convert custom binary genotype format to a text file.")
    parser.add_argument('--input_bin', required=True, help='Path to the input binary file.')
    parser.add_argument('--output_txt', required=True, help='Path to the output text file.')
    parser.add_argument('--ref_file', required=True, help='Path to the mitochondrial reference allele file.')
    parser.add_argument('--file_type', required=True, choices=['Total', 'cordup'], help="Type of file being processed, used for description.")
    
    args = parser.parse_args()

    # 读取参考基因组
    try:
        mito_ref = pd.read_csv(args.ref_file, sep='\t', names=["pos", "base"])
        ref_dict = dict(zip(mito_ref['pos'], mito_ref['base']))
    except FileNotFoundError:
        print(f"Error: Reference file not found at {args.ref_file}")
        sys.exit(1)

    # 检查输入文件
    if not os.path.exists(args.input_bin):
        print(f"Error: Input binary file not found at {args.input_bin}")
        sys.exit(1)
        
    # 获取文件大小用于进度条
    file_size = os.path.getsize(args.input_bin)

    # 处理线粒体环状基因组
    def adjust_position(pos):
        return pos + 16569 if pos <= 0 else pos

    # 二进制格式 (11个字段)
    record_format = 'IIIIdIIIIII'
    record_size = struct.calcsize(record_format)

    # 打开输入输出文件
    with open(args.input_bin, 'rb') as bin_file, \
            open(args.output_txt, 'w') as text_file:

        # 写入表头
        header = "\t".join([
            "Molecule", "Start", "Position", "Variant", "CalledBase",
            "RefBase", "FamSize", "GT_Cts", "CSS", "DB_Cts",
            "SG_Cts", "ForwardStrand", "ReverseStrand"
        ])
        text_file.write(header + "\n")

        # 进度条
        pbar = tqdm(total=file_size, unit='B', unit_scale=True, desc=f"Converting {args.file_type}")

        while True:
            # 读取二进制记录
            data = bin_file.read(record_size)
            if not data:
                break

            pbar.update(len(data))

            # 解包数据
            try:
                unpacked = struct.unpack(record_format, data)
            except struct.error:
                print(f"Warning: Could not unpack final {len(data)} bytes. File might be truncated.")
                break

            start, end, raw_position, base_idx, css, fam_size, gt_cts, db_cts, sg_cts, forward, reverse = unpacked

            # 处理环状基因组位置
            position = adjust_position(raw_position)

            # 获取碱基信息
            called_base = dna_letters[base_idx]
            ref_base = ref_dict.get(position, 'N')

            # 构建变异描述
            variant = f"chrM:{position}{ref_base}>{called_base}"
            molecule_id = f"{start}_{end}"

            # 构建输出行
            output_fields = [
                molecule_id, str(start), str(position), variant, called_base, ref_base,
                str(fam_size), str(gt_cts), f"{css:.6f}", str(db_cts),
                str(sg_cts), str(forward), str(reverse)
            ]

            text_file.write("\t".join(output_fields) + "\n")

        pbar.close()

    print(f"Conversion complete. Output: {args.output_txt}")


if __name__ == "__main__":
    main()
