import pandas as pd
from utils import *
from tqdm import tqdm
from collections import Counter
from itertools import combinations, permutations, product
import random

# pathway=pd.read_csv("Gene2Vec_Eval/PathwayCommons/PC14.All.hgnc.txt",sep='\t')
# name_id_dict = load_json('R_Annotation/name_entrez_id.json')
# anno_dict = load_json('Go_annotation/result/temp/full_annotation_dict.json')
# gene2vec_emb = np.loadtxt('go_bert/Gene2Vec_baseline/gene2vec_dim_200_iter_9.txt', dtype=str)
# gene2vec_names = [ge[0] for ge in gene2vec_emb]
# gene2vec_ids = [name_id_dict[gn] for gn in gene2vec_names if gn in name_id_dict.keys()]
# new_GGI_dict = {}
# count=0
# valid_labels = ['controls-expression-of', 'interacts-with', 'controls-phosphorylation-of', 'controls-state-change-of', 'in-complex-with', 'catalysis-precedes', 'controls-transport-of', 'consumption-controlled-by']

# for i in tqdm(range(len(pathway))):
#     obj = pathway.iloc[i]
#     gene1 = obj['PARTICIPANT_A']
#     gene2 = obj['PARTICIPANT_B']
#     label = obj['INTERACTION_TYPE']
#     if label in valid_labels and 'CHEBI' not in gene1 and 'CHEBI' not in gene2:  
#         if gene1 in name_id_dict.keys() and gene2 in name_id_dict.keys():
#             id1 = name_id_dict[gene1]
#             id2 = name_id_dict[gene2]
#             ids = id1 + "_" + id2
#             if ids not in new_GGI_dict.keys():
#                 new_GGI_dict[ids] = [label]
#             else:
#                 new_GGI_dict[ids].append(label)
# # all_gene_name = pd.Series(list(set(all_gene_name)))
# # all_gene_name.to_csv("Gene2Vec_Eval/PathwayCommons/all_gene_name.csv", index=False, header=False)
# save_json(new_GGI_dict, "Gene2Vec_Eval/PathwayCommons/new_GGI_dict.json", "w")
# new_GGI_dict = load_json("Gene2Vec_Eval/PathwayCommons/new_GGI_dict.json")
# repeat = {key: value for key, value in new_GGI_dict.items() if len(value)>1}
# unrepeat = {key: value for key, value in new_GGI_dict.items() if len(value)==1}
# gene1_list = []
# gene2_list = []
# label_list = []
# for key, value in tqdm(unrepeat.items()):
#     ids = key.split("_")
#     gene1_list.append(ids[0])
#     gene2_list.append(ids[1])
#     label_list.append(value[0])
# result_dataframe = pd.DataFrame()
# result_dataframe['gene1_entrezID'] = gene1_listcombinations_of_two = list(combinations(keys, 2))
# result_dataframe['interaction_type'] = label_list
# result_dataframe['gene2_entrezID'] = gene2_list
# result_dataframe.to_csv("Gene2Vec_Eval/PathwayCommons/PathwayCommons_GGI.csv", sep='\t', index=False)

# result_dataframe = pd.read_csv('Gene2Vec_Eval/PathwayCommons/PathwayCommons_GGI.csv', sep='\t')
# # genetype_dict = load_json("/home/yuwei/projects_backup/Gene_Ontology_PT/R_Annotation/entrezid_genetype_full.json")
# all_genes = set(list(result_dataframe['gene1_entrezID']) + list(result_dataframe['gene2_entrezID']))

# # genes_type = {str(i): genetype_dict[str(i)] for i in all_genes if str(i) in genetype_dict.keys()}
# # count = Counter(genes_type.values())
# #protein coding': 9972, 'ncRNA': 919, 'pseudo': 718, 'snoRNA': 157, 'other': 95, 'snRNA': 23, 'rRNA': 21, 'tRNA': 10, 'unknown': 4, 'scRNA': 3})
# keys = ['protein coding', 'ncRNA', 'pseudo', 'snoRNA', 'other', 'snRNA', 'rRNA', 'tRNA', 'unknown', 'scRNA']
#permutations_of_two = list(permutations(keys, 2))
# combinations_of_two = list(product(keys, repeat=2))
# pair_keys = ['_'.join(i) for i in combinations_of_two]
# count_pair_dict = {pk:0 for pk in pair_keys}
# for i in range(len(result_dataframe)):
#     object = result_dataframe.iloc[i]
#     id1 = object['gene1_entrezID']
#     id2 = object['gene2_entrezID']
#     if str(id1) in genes_type.keys() and str(id2) in genes_type.keys():
#         count_pair_dict['_'.join([genes_type[str(id1)], genes_type[str(id2)]]) ] +=1 


        
# gene_names = pd.read_csv("Benchmark_Evaluation/data/Functional_Gene_Interaction_Prediction/PathwayCommons/all_gene_name.csv", sep='\t', names=['Gene_Name'])
# name_id_dict = load_json("ID_mapping_dicts/name_to_id/name_entrez_id.json")
# gene_ids = []
# for gn in gene_names['Gene_Name']:
#     if gn in name_id_dict.keys():
#         gene_ids.append(str(name_id_dict[gn]))
#     else:
#         gene_ids.append("-1")

# gene_names['EntrezID'] = gene_ids

# entrez_refseq = load_json("ID_mapping_dicts/entrez_to_other/entrez_refseq.json")
# refseq = []
# for e in tqdm(gene_ids):
#     if str(e) in entrez_refseq.keys():
#         refseq.append(entrez_refseq[str(e)])
#     else:
#         refseq.append('-1')
#         print(e)
# gene_names['RefSeq_ID'] = refseq

# gene_names.to_csv("Benchmark_Evaluation/data/Functional_Gene_Interaction_Prediction/PathwayCommons/map.csv", sep='\t', header=True, index=False)
pc = pd.read_csv("Benchmark_Evaluation/data/Functional_Gene_Interaction_Prediction/PathwayCommons/PathwayCommons_GGI.csv", sep='\t', dtype=str)
# goa_id = pd.read_csv("Benchmark_Evaluation/input_embeddings/GOA_entrez.txt", sep='\t', header=None, dtype=str)
# gene2vec_id = pd.read_csv("Benchmark_Evaluation/input_embeddings/gene2vec_id.csv", sep='\t', header=None, dtype=str)
# final_id = set(goa_id[0]) & set(gene2vec_id[0])
all_map = pd.read_csv("Benchmark_Evaluation/input_embeddings/all_map.csv", sep='\t', dtype=str)
# all_map = all_map[all_map['EntrezID'].isin(list(final_id))]
# all_map.to_csv("Benchmark_Evaluation/input_embeddings/all_map.csv", sep='\t', header=True, index=False)
# pd.Series(list(final_id)).to_csv("Benchmark_Evaluation/input_embeddings/final_id_in_all_benchmark.txt", sep='\t', header=None, index=False)
final_id = list(set(all_map['EntrezID']))
pc = pc[pc['gene1_entrezID'].astype(str).isin(list(final_id))]
pc = pc[pc['gene2_entrezID'].astype(str).isin(list(final_id))]

train = random.sample(final_id, k=int(len(final_id)*0.55))
test = list(set(final_id) - set(train))

# print(len(pc), int(len(pc)*0.2))

# pc_train = pc[pc['gene1_entrezID'].astype(str).isin(list(train))]
# pc_train = pc_train[pc_train['gene2_entrezID'].astype(str).isin(list(train))]
# print(len(pc_train))

pc_test = pc[pc['gene1_entrezID'].astype(str).isin(list(test))]
pc_test = pc_test[pc_test['gene2_entrezID'].astype(str).isin(list(test))]
print(len(pc_test))
pc_train = pc.loc[~pc.index.isin(pc_test.index)]
print(len(pc_train))
pc_test.to_csv("Benchmark_Evaluation/data/Functional_Gene_Interaction_Prediction/PathwayCommons/test_pairs.csv", sep='\t', header=True, index=False)
pc_train.to_csv("Benchmark_Evaluation/data/Functional_Gene_Interaction_Prediction/PathwayCommons/train_pairs.csv", sep='\t', header=True, index=False)
a=1

