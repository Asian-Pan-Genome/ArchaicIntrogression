import sys
import pysam
import subprocess
import os
import numpy as np
import pandas as pd
import multiprocessing
from functools import partial

def process_sv_row(row, raw_bam, one_end_window_size, mean_depth_global=None):
    """
    处理单个SV行的函数（作为并行任务的工作单元）。
    注意：此函数中的 samtools 命令使用单线程（-@1），
    因为并行化是在SV级别（多进程）进行的。
    """
    # 从 row 对象中获取信息
    chrom = row['chrom']
    pos = row['pos']
    svlen = row['SVLEN']
    svtype = row['SVTYPE']
    block_id = row['block_id']
    ref_is_intro = (row['introgressed_0'] / (row['introgressed_0'] + row['introgressed_1'])) > (row['non_introgressed_0'] / (row['non_introgressed_0'] + row['non_introgressed_1']))

    # 1. 定义区域
    if svtype == 'INS' or svtype == 'INV':
        start, end = max(pos - one_end_window_size, 1), pos + abs(svlen) + one_end_window_size
    elif svtype == 'DEL':
        start, end = max(pos - one_end_window_size, 1), pos + one_end_window_size
    
    region = f"{chrom}:{start}-{end}"
    tmp_bam = f"./tmp/{raw_bam}.{block_id}.{region}.bam"
    depth_file = f"./tmp/{raw_bam}.{block_id}.{region}.depth.tsv"

    try:
        # 2. 运行 samtools 命令
        # 使用 -@1 是因为我们已经在进程级别上并行了，避免线程过度竞争
        subprocess.run(f"samtools view {raw_bam} {region} -Sb -o {tmp_bam} -@8", shell=True, check=True)
        subprocess.run(f"samtools index {tmp_bam} -@8", shell=True, check=True)
        
        # 将 samtools depth 的输出重定向到文件
        with open(depth_file, "w") as f_out:
            subprocess.run(f"samtools depth -a -r {region} {tmp_bam} -@8", shell=True, check=True, stdout=f_out)

        # 3. 读取深度文件并计算
        df_depth = pd.read_csv(depth_file, sep="\t", header=None, names=['chrom', 'pos', 'depth'])
        #mean_depth = df_depth['depth'].mean()
        mean_depth = mean_depth_global
        
        if svtype == 'DEL':
            mean_depth_region = df_depth[(df_depth['pos'] >= pos) & (df_depth['pos'] <= pos + abs(svlen))]['depth'].mean()
        elif svtype == 'INS':
            mean_depth_region = df_depth[(df_depth['pos'] >= pos - 10) & (df_depth['pos'] <= pos + 10)]['depth'].mean()
        elif svtype == 'INV':
            mean_depth_region_1 = df_depth[(df_depth['pos'] >= pos - 10) & (df_depth['pos'] <= pos + 10)]['depth'].mean()
            mean_depth_region_2 = df_depth[(df_depth['pos'] >= pos + abs(svlen) - 10) & (df_depth['pos'] <= pos + abs(svlen) + 10)]['depth'].mean()
            mean_depth_region = (mean_depth_region_1 + mean_depth_region_2) / 2
        
        depth_ratio = mean_depth_region / mean_depth if mean_depth != 0 else 0
        if ref_is_intro == False:
            status = "Confirmed" if depth_ratio < 0.5 else "NotConfirmed"
        else:
            status = "Confirmed" if depth_ratio >= 0.5 else "NotConfirmed"

        # 4. 返回结果字典
        return {
            'chrom': chrom,
            'pos': pos,
            'block_id': block_id,
            'fisher_pvalue': row['fisher_pvalue'],
            'SVTYPE': svtype,
            'SVLEN': svlen,
            'mean_depth': mean_depth,
            'mean_depth_region': mean_depth_region,
            'depth_ratio': depth_ratio,
            'status': status
        }
    except subprocess.CalledProcessError as e:
        print(f"Error processing {block_id}: {e.stderr.decode()}", file=sys.stderr)
        return None # 错误处理，返回None
    finally:
        # 5. 清理临时文件
        for f in [tmp_bam, f"{tmp_bam}.bai", depth_file]:
            if os.path.exists(f):
                os.remove(f)


def main():
    if len(sys.argv) != 7:
        print(f"Usage: python {sys.argv[0]} <input.freq_test.tsv> <input.bam> <one_end_window_size> <input_cpu> <mean_depth_global> <output.prefix>")
        sys.exit(1)

    # --- 数据准备（这部分不需要并行） ---
    input_pvalue_file = sys.argv[1]
    raw_bam = sys.argv[2]
    one_end_window_size = int(sys.argv[3])
    cpu_count = int(sys.argv[4])
    mean_depth_global = float(sys.argv[5])
    output_prefix = sys.argv[-1]
    
    print("Reading and filtering p-value file...")
    df_pvalue = pd.read_csv(input_pvalue_file, sep="\t", header=0)
    df_pvalue = df_pvalue.dropna(subset=['fisher_pvalue'])
    df_pvalue = df_pvalue[df_pvalue['fisher_pvalue'] < 5e-8]
    df_ins = df_pvalue[(df_pvalue['SVTYPE'] == 'INS') & (df_pvalue['SVLEN'] >= 50)]
    df_del = df_pvalue[(df_pvalue['SVTYPE'] == 'DEL') & (df_pvalue['SVLEN'] <= -50)]
    df_inv = df_pvalue[(df_pvalue['SVTYPE'] == 'INV') & (abs(df_pvalue['SVLEN']) >= 50)]
    df_sv = pd.concat([df_ins, df_del, df_inv], ignore_index=True)

    if not os.path.exists(f"./tmp"):
        os.makedirs(f"./tmp")

    # --- 并行处理 ---
    print(f"Processing {len(df_sv)} SVs in parallel...")
    
    # 将DataFrame的每一行转换为一个独立的处理任务
    tasks = [row for _, row in df_sv.iterrows()]

    # 设置进程数，通常设置为CPU核心数或略少
    # 比如 os.cpu_count() - 2，给系统留一些资源
    num_processes = min(cpu_count, len(tasks)) 
    print(f"Using {num_processes} processes.")
    
    # 使用 with 语句确保进程池被正确关闭
    with multiprocessing.Pool(processes=num_processes) as pool:
        # 使用 functools.partial 传递固定的参数 (raw_bam, one_end_window_size)
        # pool.map 只接受一个可迭代的参数，所以我们将固定的参数“冻结”到函数中
        worker_func = partial(process_sv_row, raw_bam=raw_bam, one_end_window_size=one_end_window_size, mean_depth_global=mean_depth_global)
        results = pool.map(worker_func, tasks)

    # --- 结果汇总 ---
    # 过滤掉处理失败的结果（即返回None的）
    output_data = [res for res in results if res is not None]
    
    print("Writing results to output file...")
    output_df = pd.DataFrame(output_data)
    output_df.to_csv(f'{output_prefix}.sv.confirmed.tsv', sep="\t", index=False, header=True)
    print("Done!")

if __name__ == "__main__":
    main()
