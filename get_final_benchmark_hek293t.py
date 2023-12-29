import pandas as pd

print("Loading data...")
hek293t_m6ace_df = pd.read_csv('hek293t_m6ace/Hek293T_m6aceSeq_results_annotated.csv', index_col=0)
print("Drop the items that without tx_id. (Sorted by seq_name and seq_start)")
hek293t_m6ace_df = hek293t_m6ace_df.dropna(subset=['tx_id']).sort_values(['seq_name', 'seq_start'])

print("Loading data...")
hek293t_m6ace_raw_df = pd.read_csv('hek293t_m6ace/Hek293T_m6aceSeq_results.csv').sort_values(['Chr', 'Start'])

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
