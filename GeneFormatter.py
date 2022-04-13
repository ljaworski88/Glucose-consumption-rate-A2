import openpyxl
import pandas as pd
import getpass as gp

wb = openpyxl.load_workbook('C:/Users/' + gp.getuser() + '/Dropbox/Lab/Experiments/PCR/a2-gcr(gluc_ox-days)-NP/'
                                                         'cleaned_data/PCR Analysis.xlsx')
df = pd.DataFrame()
genes = ['RPL4', 'HPRT I', 'ACTB', 'YWHAZ', 'TBP', 'HMBS', '18s', 'GAPDH', 'TIMP 1', 'TIMP 2', 'TIMP 3', 'MT 1',
         'MMP 1', 'MMP 3', 'MMP 9', 'MMP 13', 'Agg', 'Col IA1', 'Col IIA1', 'T']
# split the data and arrange it by day, the number are based off of the data's position in the excel spreadsheet
for ws in wb:
    for col in ws.iter_cols():
        col2 = [val.value for val in col]
        imp_col_val = col2[0:21], col2[21:42], col2[42:63]
        if isinstance(imp_col_val[0][0], int):
            temp_df = pd.DataFrame({'Detector': genes, 'Cq': imp_col_val[0][1:]})
            temp_df['Sample'] = imp_col_val[0][0]
            temp_df['Exp_Grp'] = ws.title
            temp_df['Day'] = 1
            df = df.append(temp_df, ignore_index=True)
        if isinstance(imp_col_val[1][0], int):
            temp_df = pd.DataFrame({'Detector': genes, 'Cq': imp_col_val[1][1:]})
            temp_df['Sample'] = imp_col_val[1][0]
            temp_df['Exp_Grp'] = ws.title
            temp_df['Day'] = 5
            df = df.append(temp_df, ignore_index=True)
        if isinstance(imp_col_val[2][0], int):
            temp_df = pd.DataFrame({'Detector': genes, 'Cq': imp_col_val[2][1:]})
            temp_df['Sample'] = imp_col_val[2][0]
            temp_df['Exp_Grp'] = ws.title
            temp_df['Day'] = 10
            df = df.append(temp_df, ignore_index=True)
df.loc[df['Detector'] == 'HPRT I', 'Detector'] = 'HPRT_I'
df.loc[df['Detector'] == 'TIMP 1', 'Detector'] = 'TIMP_1'
df.loc[df['Detector'] == 'TIMP 2', 'Detector'] = 'TIMP_2'
df.loc[df['Detector'] == 'TIMP 3', 'Detector'] = 'TIMP_3'
df.loc[df['Detector'] == 'MT 1', 'Detector'] = 'MT_1'
df.loc[df['Detector'] == 'MMP 1', 'Detector'] = 'MMP_1'
df.loc[df['Detector'] == 'MMP 3', 'Detector'] = 'MMP_3'
df.loc[df['Detector'] == 'MMP 9', 'Detector'] = 'MMP_9'
df.loc[df['Detector'] == 'MMP 13', 'Detector'] = 'MMP_13'
df.loc[df['Detector'] == 'Col IA1', 'Detector'] = 'Col_IA1'
df.loc[df['Detector'] == 'Col IIA1', 'Detector'] = 'Col_IIA1'

for id_num in set(df['Sample']):
    df.loc[df['Sample'] == id_num, 'Sample'] = 'S' + str(id_num)
# MMP 9 is being excluded due to poor primer performance and multiple products
df = df[df['Detector'] != 'MMP_9']
# Undefined is being changed to a cq of 40 for the purpose of statistics
df['Cq'].replace(to_replace='Undetermined', value=40.0, inplace=True)
df[['Sample', 'Detector', 'Cq']].to_csv('C:/Users/' + gp.getuser() + '/Dropbox/Lab/Experiments'
                                        '/PCR/a2-gcr(gluc_ox-days)-NP/cleaned_data/PCR_cts.txt',
                                        sep='\t', index=False)
pheno_df = df[['Sample', 'Exp_Grp', 'Day']]
pheno2 = pd.DataFrame()
for id_num in set(pheno_df['Sample']):
    pheno2 = pheno2.append(pheno_df[pheno_df['Sample'] == id_num].iloc[1], ignore_index=True)
    pheno2['Day'] = pheno2['Day'].astype(int)
pheno2.to_csv('C:/Users/' + gp.getuser() + '/Dropbox/Lab/Experiments/PCR/a2-gcr(gluc_ox-days)-NP/'
              'cleaned_data/PCR_phenodata.txt',
              sep='\t', index=False)
