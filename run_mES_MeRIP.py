import subprocess
import glob
import os
import sys
import pandas as pd
from multiprocessing import Pool


def run_command(command):
    print('=======================')
    print('Running: {0}'.format(command))
    p = subprocess.Popen(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    if p.returncode != 0:
        print_format_log(err.decode('utf-8'), log_type='ERROR')
        sys.exit(p.returncode)
    else:
        print_format_log(out.decode('utf-8'), log_type='OUTPUT')
        print(f'Success: {command}')


def print_format_log(logstring, log_type='OUTPUT'):
    print('-----------------------')
    print(f'{log_type}: {logstring}')
    print('-----------------------')


def main():
    run_command('Rscript get_ensembdb_mES.R')
    df = pd.read_csv('hek293t_m6ace/m6ace_hek293t_GT_extend.csv')
    contigs = df['Chr'].unique()
    commands = [f'Rscript convert_genomic2transcriptomic.R {contig}' for contig in contigs]
    os.makedirs('hek293t_m6ace/tmp', exist_ok=True)
    with Pool(50) as pool:
        pool.map(run_command, commands)
        pool.close()
        pool.join()

    files = glob.glob('hek293t_m6ace/tmp/tmp_*.csv')
    df = pd.concat([pd.read_csv(f) for f in files])
    df.to_csv('hek293t_m6ace/Hek293T_m6aceSeq_results_annotated.csv', index=False)
    run_command('python get_final_benchmark_hek293t.py')
    run_command('bedtools merge -i hek293t_m6ace/Hek293T_m6aceSeq_annotated_m6aceincluded.bed > hek293t_m6ace/Hek293T.bed')

if __name__ == '__main__':
    main()