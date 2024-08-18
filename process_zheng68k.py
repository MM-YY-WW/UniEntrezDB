import json
from utils import *
from collections import Counter
import h5py
import torch
import pandas as pd

# Map gene names in zheng68k dataset to entrez ids and save as a csv
# gene_names = pd.read_csv("Benchmark_Evaluation/data/Single_Cell_Type_Annotation/Zheng68k/gene_name.csv", sep ='\t', header=None)
# name_id_dict = load_json("ID_mapping_dicts/name_entrez_id.json")
# gene_ids = []
# for gn in gene_names[0]:
#     if gn in name_id_dict.keys():
#         gene_ids.append(str(name_id_dict[gn]))
#     else:
#         gene_ids.append("-1")
# pd.Series(gene_ids).to_csv("Benchmark_Evaluation/data/Single_Cell_Type_Annotation/Zheng68k/gene_entrez.csv", sep='\t', header=False, index=False)
# mapped = len(gene_names[0])-gene_ids.count('-1')
# print(f'there are {mapped} out of {len(gene_names[0])} mapped successfully from name to EntrezID')


# convert h5df to .pt

# hdf5_file = "Benchmark_Evaluation/input_embeddings/ontoprotein.hdf5"
# with h5py.File(hdf5_file, 'r') as f:
#     # Assuming the dataset is stored under the key 'data' and 'ids'
#     data = [list(v) for v in tqdm(list(f.values()))]
#     ids = list(f.keys())

# # Convert to PyTorch tensor with type torch.float32
# tensor_data = torch.tensor(data, dtype=torch.float32)

# # Save the tensor to a .pt file
# torch.save(tensor_data, "Benchmark_Evaluation/input_embeddings/ontoprotein.pt")

# # Convert the IDs to a DataFrame and save to a CSV file
# ids_df = pd.DataFrame(ids, columns=["ID"])
# ids_df.to_csv("Benchmark_Evaluation/input_embeddings/ids.csv", sep='\t', header=False, index=False)


emb1 = torch.load("Benchmark_Evaluation/input_embeddings/DNA_bert/dna_bert2_embedding_unget.pt")
emb2 = torch.load("Benchmark_Evaluation/input_embeddings/DNA_bert/Embedding_dna_bert2.pt")

id2 = pd.read_csv("Benchmark_Evaluation/input_embeddings/DNA_bert/EntrezID_dna_bert2.csv", sep='\t', header=None, dtype=str)
id1 = pd.read_csv("Benchmark_Evaluation/input_embeddings/DNA_bert/EntrezID_unget.csv", sep='\t', header=None, dtype=str)
a=1