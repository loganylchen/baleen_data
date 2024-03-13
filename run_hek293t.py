import glob
import os
from multiprocessing import Pool

from baleen_helper import *






def main():
    hek293t_df = get_select_hek293t_df()
    output_file = 'HEK293T/M6A.csv'
    final_output = 'HEK293T/M6A_annotated.csv'
    raw_bed = 'HEK293T/M6A.bed'
    final_bed = 'HEK293T/M6A_merged.bed'
    tmp_dir = 'tmp_hek293t'
    hek293t_df.to_csv(output_file, index=False)
    run_command('Rscript get_ensembdb_human.R')
    contigs = hek293t_df['Chr'].unique()
    os.makedirs(tmp_dir, exist_ok=True)
    commands = [f'Rscript convert_genomic2transcriptomic.R {contig} {output_file} {tmp_dir}' for contig in contigs]

    with Pool(50) as pool:
        pool.map(run_command, commands)
        pool.close()
        pool.join()

    files = glob.glob(f'{tmp_dir}/tmp_*.csv')
    df = pd.concat([pd.read_csv(f) for f in files])
    df = df.dropna(subset=['tx_id']).sort_values(['seq_name', 'seq_start'])
    df['seq_start'] = df['seq_start'].astype(int)
    df['seq_end'] = df['seq_end'].astype(int)
    df['seq_strand'] = df['seq_strand'].astype(str)
    df['seq_name'] = df['seq_name'].astype(str)
    hek293t_df['Start'] = hek293t_df['Start'].astype(int)
    hek293t_df['End'] = hek293t_df['End'].astype(int)
    hek293t_df['Strand'] = hek293t_df['Strand'].astype(str)
    hek293t_df['Chr'] = hek293t_df['Chr'].astype(str)
    df.merge(hek293t_df,how='left',left_on=['seq_name','seq_start','seq_end','seq_strand'],right_on=['Chr','Start','End','Strand'])
    df.to_csv(final_output, index=False)
    df['start'] = df['start'] - 6
    df['end'] = df['end'] + 7
    df.loc[df['start'] < 0, 'start'] = 0
    df.sort_values(['tx_id', 'start', 'end']).loc[:, ['tx_id', 'start', 'end'
                                                                    ]].to_csv(raw_bed, sep='\t', index=False, header=False)
    run_command(f'bedtools merge -i {raw_bed} > {final_bed}')

if __name__ == '__main__':
    main()







