import networkx as nx
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
import os
from utils import *
from collections import Counter

def read_gmt(file_path):
    pathways = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            pathway_name = parts[0]
            genes = parts[2:]  # Skip the first two columns which are pathway name and description
            pathways[pathway_name] = genes
    return pathways

# gmt_file_path = 'Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/msigdb.v2023.2.Hs.entrez.gmt'
# pathways = read_gmt(gmt_file_path)
# all_genes = []
# for i in list(pathways.values()):
#     all_genes += i
# gene_ids = list(set(all_genes)) #42416

# # genetype_dict = load_json("/home/yuwei/projects_backup/Gene_Ontology_PT/R_Annotation/entrezid_genetype_full.json")

# #'protein coding': 19486, 'pseudo': 13207, 'ncRNA': 8475, 'snoRNA': 708, 'other': 399, 'snRNA': 64, 'rRNA': 24, 'tRNA': 22, 'unknown': 6, 'scRNA': 3

# map = pd.DataFrame()
# map['EntrezID'] = gene_ids

# entrez_refseq = load_json("ID_mapping_dicts/entrez_to_other/entrez_refseq.json")
# refseq = []
# for e in tqdm(gene_ids):
#     if str(e) in entrez_refseq.keys():
#         refseq.append(entrez_refseq[str(e)])
#     else:
#         refseq.append('-1')
#         print(e)
# map['RefSeq_ID'] = refseq
# map.to_csv("Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/map.csv", sep='\t', header=True, index=False)


# genes_type = {str(i): genetype_dict[str(i)] for i in all_genes if str(i) in genetype_dict.keys()}

train_positive_pairs = pd.read_csv("Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/train_positive_gene_pairs.txt", sep='\t', header=None)
train_negative_pairs = pd.read_csv("Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/train_negative_gene_pairs.txt", sep='\t', header=None)
test_positive_pairs = pd.read_csv("Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/test_positive_gene_pairs.txt", sep='\t', header=None)
test_negative_pairs = pd.read_csv("Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/test_negative_gene_pairs.txt", sep='\t', header=None)

def sample_out(positive, negative, n):
    filtered_positive_df = positive[positive.groupby(0).cumcount() < n]
    filtered_negative_df = negative[negative.groupby(0).cumcount() < n]
    filtered_positive_df.columns = ['Gene1', 'Gene2', 'Co_present']
    filtered_negative_df.columns = ['Gene1', 'Gene2', 'Co_present']
    combined_df= pd.concat([filtered_positive_df, filtered_negative_df])
    return combined_df


entrez = list(pd.read_csv("Benchmark_Evaluation/input_embeddings/GOA_entrez.txt", sep='\t', header=None, dtype=str)[0])
train_positive_pairs = train_positive_pairs[train_positive_pairs[0].astype(str).isin(entrez)]
train_negative_pairs = train_negative_pairs[train_negative_pairs[0].astype(str).isin(entrez)]
test_positive_pairs = test_positive_pairs[test_positive_pairs[0].astype(str).isin(entrez)]
test_negative_pairs = test_negative_pairs[test_negative_pairs[0].astype(str).isin(entrez)]

train_positive_pairs = train_positive_pairs[train_positive_pairs[1].astype(str).isin(entrez)]
train_negative_pairs = train_negative_pairs[train_negative_pairs[1].astype(str).isin(entrez)]
test_positive_pairs = test_positive_pairs[test_positive_pairs[1].astype(str).isin(entrez)]
test_negative_pairs = test_negative_pairs[test_negative_pairs[1].astype(str).isin(entrez)]

entrez = list(pd.read_csv("Benchmark_Evaluation/input_embeddings/gene2vec_id.csv", sep='\t', header=None, dtype=str)[0])
train_positive_pairs = train_positive_pairs[train_positive_pairs[0].astype(str).isin(entrez)]
train_negative_pairs = train_negative_pairs[train_negative_pairs[0].astype(str).isin(entrez)]
test_positive_pairs = test_positive_pairs[test_positive_pairs[0].astype(str).isin(entrez)]
test_negative_pairs = test_negative_pairs[test_negative_pairs[0].astype(str).isin(entrez)]

train_positive_pairs = train_positive_pairs[train_positive_pairs[1].astype(str).isin(entrez)]
train_negative_pairs = train_negative_pairs[train_negative_pairs[1].astype(str).isin(entrez)]
test_positive_pairs = test_positive_pairs[test_positive_pairs[1].astype(str).isin(entrez)]
test_negative_pairs = test_negative_pairs[test_negative_pairs[1].astype(str).isin(entrez)]

train_df = sample_out(positive=train_positive_pairs, negative=train_negative_pairs, n=18)
train_df.to_csv("Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/train_pairs.csv", sep='\t', header=True, index=False)
print(len(train_df))
test_df = sample_out(positive=test_positive_pairs, negative=test_negative_pairs, n=18)
test_df.to_csv("Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/test_pairs.csv", sep='\t', header=True, index=False)
print(len(test_df))
a=1