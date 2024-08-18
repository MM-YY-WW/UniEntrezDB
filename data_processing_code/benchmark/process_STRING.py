from utils import *
import pandas as pd

# map = pd.read_csv("Benchmark_Evaluation/data/Protein_Protein_Interaction/STRING/data/protein.STRING_all_connected.sequences.dictionary.tsv", sep='\t', names =['Ensembl_ID', "Protein_Seq"], header=None)
ensembl_dict = load_json("ID_mapping_dicts/other_to_entrez/ensembl_entrez.json")
entrez = []


# for esb in tqdm(map["Ensembl_ID"]):
#     esb=esb.split('.')[1]
#     if esb in ensembl_dict.keys():
#         entrez.append(ensembl_dict[esb])
#     else:
#         entrez.append('-1')
# pd.Series(entrez).to_csv("Benchmark_Evaluation/data/Protein_Protein_Interaction/STRING/entrez.txt", sep='\t', header=False, index=False)
# map['EntrezID'] = entrez
# map.to_csv("Benchmark_Evaluation/data/Protein_Protein_Interaction/STRING/map.csv", sep='\t', header=True, index=False)

PPI_df = pd.read_csv("Benchmark_Evaluation/data/Protein_Protein_Interaction/STRING/data/9606.protein.actions.all_connected.txt", sep='\t', dtype=str)
final_ids = list(pd.read_csv("Benchmark_Evaluation/input_embeddings/all_map.csv", sep='\t', dtype=str)['EntrezID'])
id_a = []
for a in list(PPI_df['item_id_a']):
    if a[5:] in ensembl_dict.keys():
        id_a.append(ensembl_dict[a[5:]])
    else:
        id_a.append('-1') 
id_b = []
for b in list(PPI_df['item_id_b']):
    if b[5:] in ensembl_dict.keys():
        id_b.append(ensembl_dict[b[5:]])
    else:
        id_b.append('-1') 

PPI_df['item_id_a'] = id_a
PPI_df['item_id_b'] = id_b
PPI_df = PPI_df[PPI_df['item_id_a'] != '-1' ]
print(len(PPI_df))
PPI_df = PPI_df[PPI_df['item_id_b'] != '-1' ]
print(len(PPI_df))
PPI_df = PPI_df[PPI_df['item_id_a'].isin(final_ids)]
print(len(PPI_df))
PPI_df = PPI_df[PPI_df['item_id_b'].isin(final_ids)]
print(len(PPI_df))
PPI_df.to_csv("Benchmark_Evaluation/data/Protein_Protein_Interaction/STRING/data/input_all_cols.txt", sep='\t', header=True, index=False)
another = PPI_df[['item_id_a', 'item_id_b', 'mode']]
another.to_csv("Benchmark_Evaluation/data/Protein_Protein_Interaction/STRING/data/input_id_only.csv", sep='\t', header=True, index=False)



# pc = pd.read_csv("Benchmark_Evaluation/data/Protein_Protein_Interaction/STRING/data/input_id_only.csv", sep='\t', dtype=str)

# final_id = list(set(pc['item_id_a'] + pc['item_id_b']))

# train = random.sample(final_id, k=int(len(final_id)*0.3))
# test = list(set(final_id) - set(train))

# # print(len(pc), int(len(pc)*0.2))

# # pc_train = pc[pc['gene1_entrezID'].astype(str).isin(list(train))]
# # pc_train = pc_train[pc_train['gene2_entrezID'].astype(str).isin(list(train))]
# # print(len(pc_train))

# pc_test = pc[pc['item_id_a'].astype(str).isin(list(test))]
# pc_test = pc_test[pc_test['item_id_b'].astype(str).isin(list(test))]
# print(len(pc_test))
# pc_train = pc.loc[~pc.index.isin(pc_test.index)]
# print(len(pc_train))
# pc_test.to_csv("Benchmark_Evaluation/data/Protein_Protein_Interaction/STRING/test_pairs.csv", sep='\t', header=True, index=False)
# pc_train.to_csv("Benchmark_Evaluation/data/Protein_Protein_Interaction/STRING/train_pairs.csv", sep='\t', header=True, index=False)


a=1