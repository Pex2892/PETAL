import pandas as pd

df = pd.read_csv('analysis_all_depths.csv', delimiter='\t')

df = df.astype({'depth': 'int32', 'occurrences': 'int32'})

list_index = []
for i, v in df.iterrows():
    if 'hsa:8517' == v['starting_gene_id'] or 'hsa:8517' in str(v['starting_isoforms_id']):
        list_index.append(i)

# print(len(list_index), list_index)

# print(df.iloc[list_index, :])

df_rs = df.loc[list_index, ['depth', 'relation', 'subtype', 'reference_pathway', 'fullpath', 'occurrences']].sort_values(by=['depth', 'occurrences'], ascending=[True, False])
df_rs.to_csv('filter_hsa04014.csv', sep='\t', index=False)