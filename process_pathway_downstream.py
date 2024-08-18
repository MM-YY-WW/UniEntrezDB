import itertools
import random
from pathlib import Path
from tqdm import tqdm
from utils import *
import math
import tqdm_pathos
from multiprocessing import Pool
import multiprocessing as mp
import random

def load_pathways(msigdb_file):
    pathway_dict = {}
    with open(msigdb_file, 'r') as file:
        for line in file:
            parts = line.strip().split("\t")
            if len(parts) > 2:
                pathway_name = parts[0]
                genes = parts[2:]
                pathway_dict[pathway_name] = genes
    return pathway_dict

def convert_to_ids(genes, name_id_dict):
    return [name_id_dict.get(gene, "-1") for gene in genes if gene in name_id_dict]

def create_gene_pathway_dict(pathway_dict):
    """ Converts gene names in pathways to Entrez IDs and maps each ID to the pathways it appears in. """
    ids = list(pd.read_csv("Benchmark_Evaluation/data/all_map.csv", usecols=['EntrezID'], dtype=str, sep='\t')['EntrezID'])
    gene_pathway_dict = {}
    for pathway, genes in tqdm(pathway_dict.items()):
        for gene in genes:
            if gene in ids:
                if gene in gene_pathway_dict:
                    gene_pathway_dict[gene].add(pathway)
                else:
                    gene_pathway_dict[gene] = {pathway}
    return gene_pathway_dict

def check_pair(pair, gene_pathway_dict):
    gene1, gene2 = pair
    pathways1 = gene_pathway_dict[gene1]
    pathways2 = gene_pathway_dict[gene2]

    if set(pathways1) & set(pathways2):
        return (gene1, gene2, 1)
    else:
        return (gene1, gene2, 0)
    
def process_pairs(gene_pathway_dict):
    all_genes = list(gene_pathway_dict.keys())
    all_pairs = tqdm(list(itertools.combinations(all_genes, 2)))

    with Pool(processes=70) as pool:
        results = pool.starmap(check_pair, [(pair, gene_pathway_dict) for pair in all_pairs])

    positive_pairs = [result for result in results if result[2] == 1]
    negative_pairs = [result for result in results if result[2] == 0]

    return positive_pairs, negative_pairs


def save_pairs(filename, pairs):
    with open(filename, 'w') as file:
        for gene1, gene2, label in tqdm(pairs):
            file.write(f"{gene1}\t{gene2}\t{label}\n")

# Load pathways
def generate_pair(msigdb_file, positive_file_path, negative_file_path):
    pathway_dict = load_pathways(msigdb_file)

    gene_pathway_dict = create_gene_pathway_dict(pathway_dict)
    new_gene_pathway_dict= {str(key): ",".join(value) for key, value in gene_pathway_dict.items()}
    save_json(new_gene_pathway_dict, "Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/gene_pathway_dict", "w")
    # Generate and process pairs
    
    random.seed = 42
    gene = list(gene_pathway_dict.keys())
    train = random.sample(gene, k=int(len(gene)*0.8))
    test = list(set(gene) - set(train))
    train_pathway_dict = {trainkey: gene_pathway_dict[trainkey] for trainkey in train}
    test_pathway_dict = {testkey: gene_pathway_dict[testkey] for testkey in test}
    
    # positive_pairs, negative_pairs = process_pairs(gene_pathway_dict)
    train_positive_pairs, train_negative_pairs = process_pairs(train_pathway_dict)
    save_pairs("Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/train_positive_gene_pairs.txt", train_positive_pairs)
    save_pairs("Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/train_negative_gene_pairs.txt", train_negative_pairs)
    test_positive_pairs, test_negative_pairs = process_pairs(test_pathway_dict)
    save_pairs("Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/test_positive_gene_pairs.txt", test_positive_pairs)
    save_pairs("Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/test_negative_gene_pairs.txt", test_negative_pairs)
    # Save pairs to files
    # save_pairs(positive_file_path, positive_pairs)
    # save_pairs(negative_file_path, negative_pairs)
    
generate_pair(msigdb_file='Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/msigdb.v2023.2.Hs.entrez.gmt',
              positive_file_path='Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/positive_gene_pairs.txt',
              negative_file_path='Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/negative_gene_pairs.txt')

