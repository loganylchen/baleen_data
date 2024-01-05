import pandas as pd
import numpy as np


print("Loading data...")
hek293t_m6ace_df = pd.read_csv('hek293t_m6ace/Hek293T_m6aceSeq_results_annotated.csv', index_col=0)
print("Drop the items that without tx_id. (Sorted by seq_name and seq_start)")
hek293t_m6ace_df = hek293t_m6ace_df.dropna(subset=['tx_id']).sort_values(['seq_name', 'seq_start'])

print("Loading data...")

float_columns = [  'WT rep1 RML', 'WT rep2 RML',
       'WT rep3 RML', 'Mettl3-KO rep1 RML', 'Mettl3-KO rep2 RML',
       'Mettl3-KO rep3 RML', 'WT DESeq2 padj',
       'Mettl3-KO DESeq2 padj'
                   ]
int_columns = ['Start', 'End']
chr_columns = ['Chr',  'Strand',  'Gene', 'Annotation', 'Motif']
dtypes = dict()
for c in float_columns:
    dtypes[c]= np.float32
for c in int_columns:
    dtypes[c]= np.int64
for c in chr_columns:
    dtypes[c]= str

hek293t_m6ace_raw_df = pd.read_csv('hek293t_m6ace/Hek293T_m6aceSeq_results.csv',dtype=dtypes).dropna(axis=1,how='all').sort_values(['Chr', 'Start'])
hek293t_m6ace_raw_df['m6a'] = 0
hek293t_m6ace_raw_df.loc[(hek293t_m6ace_raw_df['WT DESeq2 padj']<0.05)&((hek293t_m6ace_raw_df['Mettl3-KO DESeq2 padj']>0.05)|(hek293t_m6ace_raw_df['Mettl3-KO DESeq2 padj'].isna())),'m6a'] = 1


print("Merging annotated data and raw data with m6a results")
merged_data = pd.merge(hek293t_m6ace_df, hek293t_m6ace_raw_df, left_on=['seq_name', 'seq_start'],
                       right_on=['Chr', 'Start'], how='left')
merged_data.dropna(axis=1,how='all',inplace=True)
merged_data.to_csv('hek293t_m6ace/Hek293T_m6aceSeq_annotated_m6aceincluded.tsv', sep='\t', index=False)

print("get the bed region!!! expand 6bp upstream and 7bp downstream")

merged_data['start'] = merged_data['start'] - 6
merged_data['end'] = merged_data['end'] + 7
merged_data.loc[merged_data['start'] < 0, 'start'] = 0

merged_data.sort_values(['tx_id', 'start', 'end']).loc[:, ['tx_id', 'start', 'end'
                                                           ]].to_csv(
    'hek293t_m6ace/Hek293T_m6aceSeq_annotated_m6aceincluded.bed', sep='\t', index=False, header=False)
