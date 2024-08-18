from Bio import Entrez, SeqIO
import pandas as pd
from tqdm import tqdm
from utils import *
from tqdm import tqdm
from multiprocessing import Pool
import time 
# refseq_entrez_dict= load_json("ID_mapping_dicts/other_to_entrez/refseq_entrez.json")
# entrez_refseq_dict= {value:key for key, value in tqdm(refseq_entrez_dict.items())}
# print(len(refseq_entrez_dict), len(entrez_refseq_dict))
# save_json(entrez_refseq_dict, "ID_mapping_dicts/entrez_to_other/entrez_refseq.json", "w")
# entrez_refseq = load_json("ID_mapping_dicts/entrez_to_other/entrez_refseq.json")
# entrez = pd.read_csv("Benchmark_Evaluation/data/Single_Cell_Type_Annotation/Zheng68k/gene_entrez.csv", sep='\t', header=None)
# refseq = []
# for e in entrez[0]:
#     if str(e) in entrez_refseq.keys():
#         refseq.append(entrez_refseq[str(e)])
#     else:
#         refseq.append('-1')
#         print(e)
# pd.Series(refseq).to_csv("Benchmark_Evaluation/data/Single_Cell_Type_Annotation/Zheng68k/refseq_id.csv", sep='\t', header=False, index=False)
# a=1
# Entrez.api_key = "9876c755e7d56d672c5b2ff013ea929a1108" #zheng68k
# Entrez.api_key ="009e62f966e0672614e96867d694f488e108" #GGI
Entrez.api_key ="5d7381e537080a511422cd9079f63bafcc08" #msigdb

def fetch_sequences(refseq_id):

    handle = Entrez.efetch(db="nucleotide", id=refseq_id, rettype="gb", retmode="text")
    records = SeqIO.parse(handle, "genbank")
    
    dna_sequence = None
    rna_sequence = None
    protein_sequence = None
    if refseq_id == '-1':
        return refseq_id, dna_sequence, rna_sequence, protein_sequence
    
    for record in records:
        dna_sequence = record.seq
        for feature in record.features:
            if feature.type == "mRNA":
                rna_sequence = feature.location.extract(record).seq
            elif feature.type == "CDS":
                if "translation" in feature.qualifiers:
                    protein_sequence = feature.qualifiers['translation'][0]
        break  # Assuming we only need the first record
    
    handle.close()
    print(f"RefSeq ID: {refseq_id} DNA Sequence: {dna_sequence[:10]} RNA Sequence: {rna_sequence[:10]} Protein Sequence: {protein_sequence[:10]}")
    return refseq_id, dna_sequence, rna_sequence, protein_sequence

def fetch_sequences_with_delay(refseq_id):
    try:
        result = fetch_sequences(refseq_id)
        print(result[0], result[1][:10], result[2][:10], result[3][:10])
        return result
    except Exception as e:
        print(e)
        return refseq_id, None, None, None

def process_refseq_ids(refseq_ids):
    with Pool(processes=5) as pool:  # Adjust the number of processes as needed
        results = list(pool.imap(fetch_sequences_with_delay, refseq_ids))
    return results

# Load RefSeq IDs
all_id = pd.read_csv("Benchmark_Evaluation/data/Single_Cell_Type_Annotation/Zheng68k/map.csv", sep='\t', names=['RefSeq_ID']) #zheng68k
# all_id = pd.read_csv("Benchmark_Evaluation/data/Functional_Gene_Interaction_Prediction/PathwayCommons/map.csv", sep='\t')  #GGI
all_id = pd.read_csv("Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/map.csv", sep='\t')  #msigdb

# Process RefSeq IDs
# results = process_refseq_ids(all_id['RefSeq_ID'])
results = process_refseq_ids(["NP_809729","YP_005029954"])
# pd.Series(results).to_csv("Benchmark_Evaluation/data/Single_Cell_Type_Annotation/Zheng68k/temp_results.csv", sep='\t') #zheng68k
# pd.Series(results).to_csv("Benchmark_Evaluation/data/Functional_Gene_Interaction_Prediction/PathwayCommons/temp_results.csv", sep='\t') #GGI
pd.Series(results).to_csv("Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/temp_results.csv", sep='\t') #msigdb
# Collect results
dna = []
protein = []
rna = []
for refseq_id, dna_seq, rna_seq, protein_seq in results:
    dna.append(dna_seq)
    protein.append(protein_seq)
    rna.append(rna_seq)

# Save results

all_id['dna'] = dna
all_id['rna'] = rna
all_id['protein'] = protein
# all_id.to_csv("Benchmark_Evaluation/data/Single_Cell_Type_Annotation/zheng68k/map.csv", sep='\t', header=True, index=False) #zheng68k
# all_id.to_csv("Benchmark_Evaluation/data/Functional_Gene_Interaction_Prediction/PathwayCommons/map.csv", sep='\t', header=True, index=False)  #GGI
# all_id.to_csv("Benchmark_Evaluation/data/Pathway_Co_present_Prediction/MsigDB/map.csv", sep='\t', header=True, index=False)  #msigdb
