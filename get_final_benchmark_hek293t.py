import pandas as pd
import numpy as np


print("Loading data...")
hek293t_m6ace_df = pd.read_csv('hek293t_m6ace/Hek293T_m6aceSeq_results_annotated.csv', index_col=0)
print("Drop the items that without tx_id. (Sorted by seq_name and seq_start)")
hek293t_m6ace_df = hek293t_m6ace_df.dropna(subset=['tx_id']).sort_values(['seq_name', 'seq_start'])
hek293t_m6ace_df['seq_start'] = hek293t_m6ace_df['seq_start'].astype(int)
hek293t_m6ace_df['seq_end'] = hek293t_m6ace_df['seq_end'].astype(int)
hek293t_m6ace_df['seq_strand'] = hek293t_m6ace_df['seq_strand'].astype(str)
hek293t_m6ace_df['seq_name'] = hek293t_m6ace_df['seq_name'].astype(str)

hek293t_m6ace_GT_extend_df = pd.read_csv('hek293t_m6ace/m6ace_hek293t_GT_extend.csv')
hek293t_m6ace_GT_extend_df['Start'] = hek293t_m6ace_GT_extend_df['Start'].astype(int)
hek293t_m6ace_GT_extend_df['End'] = hek293t_m6ace_GT_extend_df['End'].astype(int)
hek293t_m6ace_GT_extend_df['Strand'] = hek293t_m6ace_GT_extend_df['Strand'].astype(str)
hek293t_m6ace_GT_extend_df['Chr'] = hek293t_m6ace_GT_extend_df['Chr'].astype(str)
hek293t_m6ace_df.merge(hek293t_m6ace_GT_extend_df,how='left',left_on=['seq_name','seq_start','seq_end','seq_strand'],right_on=['Chr','Start','End','Strand']).to_csv(
    'hek293t_m6ace/Hek293T_m6aceSeq_results_annotated_withMotif.csv',index=False)

print("get the bed region!!! expand 6bp upstream and 7bp downstream")

hek293t_m6ace_df['start'] = hek293t_m6ace_df['start'] - 6
hek293t_m6ace_df['end'] = hek293t_m6ace_df['end'] + 7
hek293t_m6ace_df.loc[hek293t_m6ace_df['start'] < 0, 'start'] = 0

hek293t_m6ace_df.sort_values(['tx_id', 'start', 'end']).loc[:, ['tx_id', 'start', 'end'
                                                           ]].to_csv(
    'hek293t_m6ace/Hek293T_m6aceSeq_annotated_m6aceincluded.bed', sep='\t', index=False, header=False)
