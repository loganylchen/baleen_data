import subprocess
import sys
import pandas as pd


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

def get_hek293t_m6ace_df():
    raw_hek293t_df = pd.read_excel('HEK293T/m6ace.xlsx')
    raw_hek293t_df['Chr'] = raw_hek293t_df['Chr'].str.replace('chr', '')
    raw_hek293t_df['WT RML average'] = raw_hek293t_df.loc[:, ['WT rep1 RML', 'WT rep2 RML', 'WT rep3 RML']].mean(axis=1)
    raw_hek293t_df['Mettl3-KO RML average'] = raw_hek293t_df.loc[:, ['Mettl3-KO rep1 RML', 'Mettl3-KO rep2 RML', 'Mettl3-KO rep3 RML']].mean(axis=1)
    raw_hek293t_df['WT-KO'] = raw_hek293t_df['WT RML average'] / raw_hek293t_df['Mettl3-KO RML average']
    raw_hek293t_df['mettl3-m6a'] = 0
    raw_hek293t_df.loc[raw_hek293t_df['WT-KO'] >= 4, 'mettl3-m6a'] = 1
    return raw_hek293t_df

def get_hek293t_glory_df():
    raw_glory_df = pd.read_excel('HEK293T/glory.xlsx')
    raw_glory_df['Chr'] = raw_glory_df['Chr'].str.replace('chr', '')
    return raw_glory_df

def get_select_hek293t_df():
    m6ace_seq_df = get_hek293t_m6ace_df()
    glory_df = get_hek293t_glory_df()
    df = m6ace_seq_df.merge(glory_df, left_on=['Chr', 'End', 'Strand'], right_on=['Chr', 'Sites', 'Strand'], how='left')
    df = df.loc[((df['mettl3-m6a'] == 1) & (~df['Cluster_info'].isna())) | ((df['mettl3-m6a'] == 0) & (df['Cluster_info'].isna()))]
    return df




