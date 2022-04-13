import openpyxl
import pandas as pd

wb = openpyxl.load_workbook('C:/Users/Student2/Dropbox/Lab/Experiments/PCR/a2-housekeeping_genes-NP/Cleaned Data/'
                            'PCR Analysis for Genorme.xlsx')
df = pd.DataFrame()
genes = ['RPL4', 'HPRT I', 'ACTB', 'YWHAZ', 'TBP', 'HMBS', '18s', 'GAPDH']
for ws in wb:
    for col in ws.iter_cols():
        col2 = [val.value for val in col]
        imp_col_val = col2[1:10], col2[11:20], col2[21:30]
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
for id_num in set(df['Sample']):
    df.loc[df['Sample'] == id_num, 'Sample'] = 'S' + str(id_num)
df[['Sample', 'Detector', 'Cq']].to_csv('C:/Users/Student2/Dropbox/Lab/Experiments/PCR/a2-housekeeping_genes-NP/'
                                        'Cleaned Data/GeneNorm_Cqs.txt',
                                        sep='\t', index=False)
pheno_df = df[['Sample', 'Exp_Grp', 'Day']]
pheno2 = pd.DataFrame()
for id_num in set(pheno_df['Sample']):
    pheno2 = pheno2.append(pheno_df[pheno_df['Sample'] == id_num].iloc[1], ignore_index=True)
    pheno2['Day'] = pheno2['Day'].astype(int)
pheno2.to_csv('C:/Users/Student2/Dropbox/Lab/Experiments/PCR/a2-housekeeping_genes-NP/'
              'Cleaned Data/GeneNorm_PhenoData.txt',
              sep='\t', index=False)
