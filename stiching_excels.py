import pandas as pd

df1 = pd.read_excel('lung_cancer_1000_sequences.xlsx')
df2 = pd.read_excel('lung_protein_sequences.xlsx')

merged_df = pd.concat([df1, df2], ignore_index=True)

merged_df.to_excel('combined_file.xlsx', index=False)
