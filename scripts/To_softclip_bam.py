#!/usr/bin/env python3
import pysam
import argparse


def extract_softclip_reads(input_bam, output_with_softclip, output_without_softclip):
    """
    实际可用的解决方案：可能会有少量误判，但效率较高
    """
    # 打开输入BAM文件
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        # 创建输出BAM文件
        with pysam.AlignmentFile(output_with_softclip, "wb", template=infile) as softclip_out, \
                pysam.AlignmentFile(output_without_softclip, "wb", template=infile) as non_softclip_out:

            # 第一步：收集所有带有softclip的reads的query_name
            softclip_names = set()
            for read in infile:
                if has_softclip(read):
                    softclip_names.add(read.query_name)

            # 重置文件指针
            infile.reset()

            # 第二步：处理所有reads
            for read in infile:
                if read.query_name in softclip_names:
                    softclip_out.write(read)
                else:
                    # 检查配对read是否有softclip
                    if not read.mate_is_unmapped:
                        # 由于我们无法高效获取配对read的softclip状态，这里做一个假设：
                        # 如果当前read没有softclip，且配对read有softclip，那么当前read应该被写入softclip
                        # 但由于我们无法高效检查，我们只能接受这种不完美的情况
                        # 或者我们可以选择将所有配对read也写入softclip（不推荐）
                        pass

                    # 默认情况下，写入non_softclip文件
                    non_softclip_out.write(read)

                    # 如果需要更精确，可以在这里再次扫描文件来检查配对read
                    # 但这会非常慢，不推荐


def has_softclip(read):
    """检查read是否有softclip"""
    if read.cigartuples is None:
        return False
    return any(op == 4 for op, _ in read.cigartuples)  # 4代表softclip


def main():
    parser = argparse.ArgumentParser(description='Extract reads with softclip and their mates from a BAM file')
    parser.add_argument('input_bam', help='Input BAM file')
    parser.add_argument('output_with_softclip', help='Output BAM file for reads with softclip')
    parser.add_argument('output_without_softclip', help='Output BAM file for reads without softclip')

    args = parser.parse_args()

    # 检查输入文件是否存在
    try:
        with open(args.input_bam, 'rb') as f:
            pass
    except IOError:
        print(f"Error: Input file {args.input_bam} does not exist")
        return

    # 运行提取函数
    extract_softclip_reads(args.input_bam, args.output_with_softclip, args.output_without_softclip)
    print(f"Processing complete. Results saved to {args.output_with_softclip} and {args.output_without_softclip}")


if __name__ == "__main__":
    main()

