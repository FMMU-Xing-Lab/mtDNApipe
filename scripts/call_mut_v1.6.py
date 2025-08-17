import argparse
import pandas as pd
from collections import defaultdict

'''
usage:python3 call_mut_v1.5.py <sample> <output_dir> <output_file> <Terminal bp> --min_freq <threshold> --min_reads <threshold>
'''

def main():
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description='Mitochondrial Mutation Filtering Pipeline')
    parser.add_argument('sample_name', help='Sample name')
    parser.add_argument('output_dir', help='Output file')
    parser.add_argument('output_file', help='Output file')
    parser.add_argument('n', type=int, help='Boundary distance filter value')
    parser.add_argument('--min_freq', type=float, default=0.1,
                        help='Minimum mutation frequency threshold (%%), default=0.1')
    parser.add_argument('--min_reads', type=int, default=5,
                        help='Minimum supporting reads threshold, default=5')
    args = parser.parse_args()

    sample_name = args.sample_name
    output_dir = args.output_dir
    output_file = args.output_file
    n = args.n
    min_freq = args.min_freq
    min_reads = args.min_reads
    cutoff = 16569  # 线粒体基因组调整阈值

    # 定义重复区域（根据题目要求）
    repeat_regions = [
        (302, 316),  # Region 1
        (513, 526),  # Region 2
        (566, 573),  # Region 3
        (16181, 16194)  # Region 4
    ]

    qualified_file = f"{output_dir}/{sample_name}/{sample_name}.QualifiedTotalCts"
    input_files = {
        'Total': f"{output_dir}/{sample_name}/{sample_name}.RawGenotypes.Total.txt",
        'cordup': f"{output_dir}/{sample_name}/{sample_name}.RawGenotypes.cordup.txt"
    }

    # 读取覆盖深度文件
    try:
        coverage_df = pd.read_csv(qualified_file, sep='\t',
                                  names=["Position", "Total", "cordup"])
        coverage_dict = coverage_df.set_index('Position').to_dict('index')
    except FileNotFoundError:
        print(f"Error: Qualified total counts file {qualified_file} not found")
        sys.exit(1)

    # 处理每种严格类型
    for strict_type, input_file in input_files.items():
        try:
            # 读取原始突变文件
            df = pd.read_csv(input_file, sep='\t',header=0,
                             names=["Molecule", "Start", "Position", "Variant", "CalledBase",
                                    "RefBase", "FamSize", "GT_Cts", "CSS", "DB_Cts",
                                    "SG_Cts", "ForwardStrand", "ReverseStrand"])
        except FileNotFoundError:
            print(f"Warning: {strict_type} input file {input_file} not found, skipping")
            continue

        # 初始化排除统计字典
        excluded_depth = defaultdict(int)
        #print(df['Molecule'])
        # ===================== 功能1：处理分子出现次数限制 =====================
        # 解析分子坐标
        df[['fragment_start', 'fragment_end']] = df['Molecule'].str.split('_', expand=True).astype(int)

        # 计算分子长度
        df['fragment_length'] = df['fragment_end'] - df['fragment_start']

        # 统计分子出现次数
        molecule_counts = df.groupby('Molecule').size().reset_index(name='count')

        # 合并长度信息
        molecule_counts[['fs', 'fe']] = molecule_counts['Molecule'].str.split('_', expand=True).astype(int)
        molecule_counts['length'] = molecule_counts['fe'] - molecule_counts['fs']

        # 计算允许的最大出现次数（ceil(length/50)）
        molecule_counts['max_allowed'] = (molecule_counts['length'] + 49) // 50

        # 筛选违规分子
        invalid_molecules = molecule_counts[molecule_counts['count'] > molecule_counts['max_allowed']]['Molecule']

        # 记录被排除的突变并统计深度
        excluded_by_molecule = df[df['Molecule'].isin(invalid_molecules)]
        for pos in excluded_by_molecule['Position']:
            excluded_depth[pos] += 1

        # 过滤违规分子
        df = df[~df['Molecule'].isin(invalid_molecules)]

        # ===================== 功能2：处理边界过滤 =====================
        # 坐标调整
        mask = ((df['fragment_start'] > cutoff) & (df['fragment_end'] > cutoff)) | \
               ((df['fragment_start'] <= cutoff) & (df['fragment_end'] > cutoff))
        df['adjusted_start'] = df['fragment_start'].where(~mask, df['fragment_start'] - cutoff)
        df['adjusted_end'] = df['fragment_end'].where(~mask, df['fragment_end'] - cutoff)

        # 应用边界过滤
        boundary_mask = ~((df['Position'] - df['adjusted_start'] < n) |
                          (df['adjusted_end'] - df['Position'] < n))

        # 记录被排除的突变
        excluded_by_boundary = df[~boundary_mask]
        for pos in excluded_by_boundary['Position']:
            excluded_depth[pos] += 1

        # 应用边界过滤条件
        filtered_df = df[boundary_mask & (df['CalledBase'] != df['RefBase'])]

        # ===================== 功能3：双链支持过滤 =====================
        # 根据FamSize应用不同的双链支持条件
        strand_mask = (
            # 当FamSize >= 3时：要求双链支持（Forward和Reverse都等于1）
                ((filtered_df['FamSize'] >= 3) &
                 (filtered_df['ForwardStrand'] >= 1) &
                 (filtered_df['ReverseStrand'] >= 1)) |

                # 当FamSize < 3时：不要求双链支持（总是True）
                (filtered_df['FamSize'] < 3)
        )

        # 记录被双链支持过滤排除的突变
        excluded_by_strand = filtered_df[~strand_mask]
        for pos in excluded_by_strand['Position']:
            excluded_depth[pos] += 1

        # 应用双链支持过滤
        final_filtered_df = filtered_df[strand_mask]

        # ===================== 计算突变频率 =====================
        if final_filtered_df.empty:
            print(f"No valid mutations found for {strict_type}")
            continue

        # 计算突变计数
        mutation_counts = final_filtered_df.groupby(['Position', 'Variant']).size().reset_index(name='Count')
        #mutation_counts[['ref', 'alt']] = mutation_counts['Variant'].str.split('_', expand=True).iloc[:, 1:3]
        mutation_counts['ref'] = mutation_counts['Variant'].str.extract(r'(\w)>')[0]
        mutation_counts['alt'] = mutation_counts['Variant'].str.extract(r'>(\w)')[0]
        
        # 调整覆盖深度（减去排除的突变数）
        mutation_counts['TotalDepth'] = mutation_counts['Position'].apply(
            lambda pos: coverage_dict.get(pos, {}).get(strict_type, 0) - excluded_depth.get(pos, 0)
        )

        # 计算突变频率
        mutation_counts['Frequency(%)'] = mutation_counts.apply(
            lambda row: round(row['Count'] / row['TotalDepth'] * 100, 4) if row['TotalDepth'] > 0 else 0.0,
            axis=1
        )

        # 记录各种过滤前的突变计数
        initial_mutations = len(mutation_counts)

        # ===================== 应用新的过滤条件 =====================
        # 1. 过滤支持reads数不足的突变
        read_filter_mask = mutation_counts['Count'] >= min_reads
        filtered_by_reads = mutation_counts[~read_filter_mask]
        mutation_counts = mutation_counts[read_filter_mask]

        # 2. 过滤频率过低的突变
        freq_filter_mask = mutation_counts['Frequency(%)'] >= min_freq
        filtered_by_freq = mutation_counts[~freq_filter_mask]
        mutation_counts = mutation_counts[freq_filter_mask]

        '''
        # 3. 过滤氧化损伤突变 (VAF < 10% 且 突变类型为C>A/G>T)
        ox_mask = (
                          ((mutation_counts['ref'] == 'C') & (mutation_counts['alt'] == 'A')) |
                          ((mutation_counts['ref'] == 'G') & (mutation_counts['alt'] == 'T'))
                  ) & (mutation_counts['Frequency(%)'] < 10.0)
        filtered_by_ox = mutation_counts[ox_mask]
        mutation_counts = mutation_counts[~ox_mask]
        '''

        # 4. 过滤重复区域的突变
        in_repeat_region = False
        repeat_filter_mask = mutation_counts['Position'].apply(
            lambda pos: any(start <= pos <= end for start, end in repeat_regions)
        )
        filtered_by_repeat = mutation_counts[repeat_filter_mask]
        mutation_counts = mutation_counts[~repeat_filter_mask]

        # ===================== 保存结果 =====================
        # 确保结果不为空
        if not mutation_counts.empty:
            #output_path = f"{output_dir}/{sample_name}/{sample_name}.{strict_type}.filtered_mutations_v2.txt"
            mutation_counts[['Variant', 'Position', 'ref', 'alt', 'Count', 'TotalDepth', 'Frequency(%)']] \
                .to_csv(output_file, sep='\t', index=False)

        print(f"\n===== {strict_type} Results Summary =====")
        print(f"Total valid mutations (after initial filters): {initial_mutations}")
        print(f"Additional filtering summary:")
        print(f"  - Excluded by min reads (<{min_reads}): {len(filtered_by_reads)}")
        print(f"  - Excluded by min frequency (<{min_freq}%): {len(filtered_by_freq)}")
        #print(f"  - Excluded as oxidative damage: {len(filtered_by_ox)}")
        print(f"  - Excluded in repeat regions: {len(filtered_by_repeat)}")
        print(f"Final mutations: {len(mutation_counts)}")

        if not mutation_counts.empty:
            print(f"Results saved to: {output_file}")


if __name__ == "__main__":
    main()
