import os
import pandas as pd
import csv
import os
import json
from tqdm import tqdm
import argparse
import numpy as np
import multiprocessing as mp
import tqdm_pathos
from multiprocessing import Pool
from utils import *


def get_gene_govec(goannos, args):
    embedding = np.load(args.go_embedding_path)
    id = np.load(args.go_id_path)
    indices = [list(id).index(anno) if anno in list(id) else -1 for anno in goannos]
    extracted_elements = embedding[indices]
    #qualifier = [[qua] for idx, qua in zip(indices, qualifier) if idx != -1]
    mean_pooled = np.mean(extracted_elements, axis=0)
    return mean_pooled


def generate_gene_govec_mp(args, noIEA):
    # gene_EntrezID, gene_GOAnno, gene_Qualifier = get_goall(args)
    gene_EntrezID, gene_GOAnno = get_gofromcsv(noIEA)
    result_folder_path = os.path.join(args.result_folder, f'noIEA:{noIEA}_noroot')
    if not os.path.exists(result_folder_path):
        os.makedirs(result_folder_path)
    print('\nStart Generating embedding, the progress bar does not stuck, may take around 30 minutes')
    results = tqdm_pathos.starmap(get_gene_govec, [(anno, args) for anno in gene_GOAnno], n_cpus=int(args.ratio_cpu*mp.cpu_count()))
    #results = tqdm_pathos.starmap(get_gene_govec, [(anno, qua, args) for anno, qua in zip(gene_GOAnno[:2], gene_Qualifier[:2])], n_cpus=1)
    result_to_save = np.array([[int(id)] + list(res) for id, res in zip(gene_EntrezID, results)])
    file_path = os.path.join(result_folder_path,f'{args.method}_{len(gene_EntrezID)}.txt')
    np.savetxt(file_path,result_to_save, fmt='%s', delimiter='\t')
    
    print(f'\nFinished, file saves to {file_path}')

# def get_gofromcsv(noIEA):
#     all_anno_csv = pd.read_csv("Go_annotation/extra_go/processed_uni_all_database.csv", sep='\t', usecols=["EntrezID", "Evidence_Code", "GO_ID"])
#     if noIEA:
#         entrez_set = set(all_anno_csv[all_anno_csv['Evidence_Code']!= 'IEA']['EntrezID'])
#     else:
#         entrez_set = set(all_anno_csv['EntrezID'])
#     entrez_list = []
#     annotation_list = []
#     for es in tqdm(entrez_set):
#         if ';' not in es:
#             entrez_list.append(es)
#             annotation_list.append(list(all_anno_csv[all_anno_csv['EntrezID']!= es]['GO_ID']))
#     df = pd.DataFrame()
#     df['EntrezID'] = entrez_list
#     combine_anno = [','.join(i) for i in annotation_list]
#     if noIEA:
#         result_path = "Go_annotation/extra_go/Embeddings/noIEA_annotation.csv"
#         df['noIEA_Annotations'] = combine_anno
#     else:
#         result_path = "Go_annotation/extra_go/Embeddings/All_annotation.csv"
#         df['All_Annotations'] = combine_anno
#     df.to_csv(result_path, sep='\t', header=True, index=False)
#     return entrez_list, annotation_list 

all_anno_csv = pd.read_csv("Go_annotation/extra_go/downstream_subfile.csv", sep='\t', usecols=["EntrezID", "Evidence_Code", "GO_ID"], dtype=str)

def process_entrez(es, noIEA):
    global all_anno_csv
    if noIEA:
        go_ids = list(all_anno_csv[(all_anno_csv['EntrezID'] == str(es)) & (all_anno_csv['Evidence_Code'] != 'IEA')]['GO_ID'])
    else:
        go_ids = list(all_anno_csv[all_anno_csv['EntrezID'] == str(es)]['GO_ID'])
    go_ids = list(set(go_ids) - {'GO:008150', 'GO:0003674', "GO:0005575"})
    if go_ids:
        print(list(set(go_ids)))
        return es, list(set(go_ids))
    print(es)
    return None

def get_gofromcsv(noIEA):
    # all_anno_csv = pd.read_csv("Go_annotation/extra_go/processed_uni_all_database.csv", sep='\t', usecols=["EntrezID", "Evidence_Code", "GO_ID"])
    entrez_set = list(pd.read_csv("/home/yuwei/data/projects_backup/UniEntrezDB/Benchmark_Evaluation/data/all_entrez.txt", sep='\t', header=None)[0])
    # all_anno_csv = all_anno_csv[all_anno_csv['EntrezID'].isin(entrez_set)]
    # all_anno_csv.to_csv("Go_annotation/extra_go/downstream_subfile.csv", sep='\t', header=True, index=False)
    with Pool(processes=mp.cpu_count()) as pool:
        results = list(tqdm(pool.starmap(process_entrez, [(es,  noIEA) for es in entrez_set]), total=len(entrez_set)))
    
    # Filter out None results
    results = [res for res in results if res is not None]

    entrez_list, annotation_list = zip(*results) if results else ([], [])
    
    df = pd.DataFrame()
    df['EntrezID'] = entrez_list
    combine_anno = [','.join(i) for i in annotation_list]
    if noIEA:
        result_path = "Go_annotation/extra_go/Embeddings/noIEA_noroot_annotation.csv"
        df['noIEA_Annotations'] = combine_anno
    else:
        result_path = "Go_annotation/extra_go/Embeddings/All_noroot_annotation.csv"
        df['All_Annotations'] = combine_anno
    df.to_csv(result_path, sep='\t', header=True, index=False)
    return list(entrez_list), list(annotation_list)

def get_goall(args):
    """_summary_
        Obtain GO annotation for each gene id 
        return three list:
        GeneID, GO annotation and qualifier is positive or negative
    Args:
        args (_type_): _description_
    """
    evidence_levels = [l.strip() for l in args.evidence_level.split(',')]
    object_types = [t.strip() for t in args.object_type.split(',')]
    qualifier = [q.strip() for q in args.qualifier.split(',')]
    evidence_level_symbol = {'-1': ['TAS', 'ND', 'IC', 'NAS'],
                              '0': ['IEA'],
                              '1': ['ISO', 'IGC', 'ISM', 'ISS', 'RCA', 'ISA'],
                              '2': ['IGI', 'HDA', 'EXP', 'IMP', 'IKR', 'IPI', 'HGI', 'HEP', 'IEP', 'IDA', 'HTP', 'IBA', 'HMP']
                            }
    object_type_dict = {
            'protein': ["protein","protein_complex"],
        
                'DNA': ["sense_intronic_ncRNA_gene","miRNA_gene","sense_overlap_ncRNA_gene","RNase_P_RNA_gene","telomerase_RNA_gene",
                        "lincRNA_gene","pseudogene","tRNA_gene","SRP_RNA_gene","protein_coding_gene","RNase_MRP_RNA_gene","gene_segment",
                        "gene","ncRNA_gene","rRNA_gene","snRNA_gene","scRNA_gene","transposable_element_gene","snoRNA_gene","lncRNA_gene"],
            
                'RNA': ["piRNA","snoRNA","tRNA","snRNA","hammerhead_ribozyme","RNA","scRNA","RNase_P_RNA","lnc_RNA","telomerase_RNA",
                        "ncRNA","RNase_MRP_RNA","guide_RNA","ribozyme","rRNA","miRNA","antisense_lncRNA","mRNA","antisense_RNA",
                        "tmRNA","SRP_RNA"],
            
             'others': ["biological region","gene_product"]
                        }         
    qualifier_level_dict = {
        'pos': ["enables","acts_upstream_of_or_within_positive_effect","acts_upstream_of_or_within","located_in","colocalizes_with",
                "part_of","acts_upstream_of_positive_effect","involved_in","acts_upstream_of","is_active_in",
                "contributes_to","acts_upstream_of_negative_effect","acts_upstream_of_or_within_negative_effect",],
        
        'neg': ["NOT|located_in","NOT|is_active_in","NOT|involved_in","NOT|acts_upstream_of_or_within_negative_effect",
                "NOT|acts_upstream_of_or_within_positive_effect","NOT|contributes_to","NOT|part_of","NOT|enables"
                "NOT|acts_upstream_of","NOT|colocalizes_with","NOT|acts_upstream_of_or_within", ]
    }
    target_evidence_codes = []
    for level in evidence_levels:
        target_evidence_codes += evidence_level_symbol[level]
    
    target_object_types = []
    for type in object_types:
        target_object_types += object_type_dict[type]
                
    target_qualifier = []
    for type in qualifier:
        target_qualifier += qualifier_level_dict[type]      
             
    anno_dict = load_json(args.DBID_EntrezID_fulldict_path)
    #gene_entrezID : [goid, goid, goid...]
    result_gene_EntrezID = []
    result_gene_GOAnno = []
    result_gene_Qualifier = []
    anno_count=0
    gene_count=0
    print(f"\n===Start to obtain chosen annotation: \n evidence_level{evidence_levels} \n object_types: {object_types} \n qualifier: {qualifier} \n")
    for geneid, type_dict in tqdm(anno_dict.items()):
        gene_anno_list = []
        gene_anno_qualifier_list = []
        for type, evidence_dict in type_dict.items():
            if type in target_object_types:
                for evidence, qualifier_dict in evidence_dict.items():
                    if evidence in target_evidence_codes:
                        for qualifier, goid_list in qualifier_dict.items():
                            if qualifier in target_qualifier:
                                gene_anno_list += goid_list
                                anno_count += len(goid_list)
                                if qualifier in qualifier_level_dict['neg']:
                                    gene_anno_qualifier_list += [-1]*len(goid_list)
                                else:
                                    gene_anno_qualifier_list += [1]*len(goid_list)
                                assert len(gene_anno_list) == len(gene_anno_qualifier_list)
        if gene_anno_list != []:
            gene_count +=1
            result_gene_EntrezID.append(geneid)
            result_gene_GOAnno.append(gene_anno_list)
            result_gene_Qualifier.append(gene_anno_qualifier_list)
                    
    print(f'\nThere are {anno_count} annotations used for {gene_count} genes embedding generation')
    obsolete_to_new = load_json(args.obsolete_to_new_dict)
    print("Replace Obsolete GO terms....")
    for m in range(len(result_gene_GOAnno)):
        GOAnno = result_gene_GOAnno[m]
        for n in range(len(GOAnno)):
            if GOAnno[n] in obsolete_to_new.keys():
                GOAnno[n] = obsolete_to_new[GOAnno[n]]
        GOAnno = list(set(GOAnno))
        result_gene_GOAnno[m] = GOAnno

    if args.augment_go_anno !='':
        print("Augment GO Annotations...")
        augment_go_anno_dict=load_json(args.augment_go_anno)
        for i in range(len(result_gene_GOAnno)):
            GOAnno = result_gene_GOAnno[i]
            new_go = []
            for j in range(len(GOAnno)):
                new_go += augment_go_anno_dict[GOAnno[j]]
            result_gene_GOAnno[i] = list(set(GOAnno+new_go))
            

        
    # parser.add_argument("--obsolete_to_new_dict", type=str, default='Go_annotation/result/temp/obsolete_to_new.json', help='path of dictionary to replace the obsolete go id')
    # parser.add_argument("--augment_go_anno", type=str, default="", help='the path to augment the go annotation to all related nodes, (Go_annotation/result/temp/augment_go_dict.json)')
    # do not use result_gene_qualifier, need to set gene qualifier filter to pos only because the negative ones is proven not related to this gene or gene product
    return result_gene_EntrezID, result_gene_GOAnno, result_gene_Qualifier

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--go_id_path', type=str, default='Node_embedding_mxbai/go_id.npy', 
                        help="The path to the numpy file that saves the go-term id")
    parser.add_argument('--result_folder', type=str, default='Go_annotation/extra_go/Embeddings', 
                        help="The path to the numpy file that saves the go-term id")
    parser.add_argument('--go_embedding_path', type=str, default='Node_embedding_mxbai/47727.npy', 
                        help="The path to the numpy file that saves all the go-term embedding")
    parser.add_argument('--method', type=str, default='mean',
                        help="currently only has mean method, TODO more model")
    parser.add_argument('--DBID_EntrezID_fulldict_path', type=str, default='Go_annotation/result/temp/full_annotation_dict.json',
                        help="path to the annotation dictionary")
    # parser.add_argument('--mp', type=int, default=0,
    #                     help="Use multiprocessing or not")
    parser.add_argument('--ratio_cpu', type=float, default=1,
                        help='how many percent of cpu cores will be used to perform the multiprocessing')
    # parser.add_argument('--interval', type=int, default=100,
    #                     help="save when processing how many samples")
    parser.add_argument('--evidence_level', type=str, default='-1, 0, 1, 2',
                        help='consider which evidence based on evidence level, input the level and separate by comma\
                        default = "-1, 0, 1, 2"\
                        -1: the confidence in this level is uncertain, evidence code contains: TAS, ND, IC, NAS\
                        0: the confidence in this level is not confident, evidence code contains: IEA\
                        1: the confidence in this level is kind of confident, evidence code contains: ISO, IGC, ISM, ISS, RCA, ISA\
                        2: the confidence in this level is very confident, evidence code contains: IGI, HDA, EXP, IMP, IKR, IPI, HGI, HEP, IEP, IDA, HTP, IBA, HMP\ '
                        )
    parser.add_argument('--object_type', type=str, default='protein, DNA, RNA, others',
                        help='consider which type of object of the gene, input a string of the type to be selected separate by comma\
                            default="protein, DNA, RNA, others"\
                        "protein": ["protein","protein_complex"],\
                            "DNA": ["sense_intronic_ncRNA_gene","miRNA_gene","sense_overlap_ncRNA_gene","RNase_P_RNA_gene","telomerase_RNA_gene",\
                                    "lincRNA_gene","pseudogene","tRNA_gene","SRP_RNA_gene","protein_coding_gene","RNase_MRP_RNA_gene","gene_segment",\
                                    "gene","ncRNA_gene","rRNA_gene","snRNA_gene","scRNA_gene","transposable_element_gene","snoRNA_gene","lncRNA_gene"],\
                            "RNA": ["piRNA","snoRNA","tRNA","snRNA","hammerhead_ribozyme","RNA","scRNA","RNase_P_RNA","lnc_RNA","telomerase_RNA",\
                                    "ncRNA","RNase_MRP_RNA","guide_RNA","ribozyme","rRNA","miRNA","antisense_lncRNA","mRNA","antisense_RNA",\
                                    "tmRNA","SRP_RNA"],\
                        "others" : ["biological region","gene_product"]\
                        ')
    parser.add_argument('--qualifier', type=str, default='pos',
                        help='which kind of qualifier is chosen, separate by comma\
                        default="pos, neg"\
                        "pos":[ "enables","acts_upstream_of_or_within_positive_effect","acts_upstream_of_or_within","located_in","colocalizes_with",\
                                "part_of","acts_upstream_of_positive_effect","involved_in","acts_upstream_of","is_active_in",\
                                "contributes_to","acts_upstream_of_negative_effect","acts_upstream_of_or_within_negative_effect"]\
                        "neg":[ "NOT|located_in","NOT|is_active_in","NOT|involved_in", "NOT|acts_upstream_of_or_within_negative_effect",\
                                "NOT|acts_upstream_of_or_within_positive_effect","NOT|contributes_to","NOT|part_of","NOT|enables"\
                                "NOT|acts_upstream_of","NOT|colocalizes_with","NOT|acts_upstream_of_or_within", ]\
                        ')
    parser.add_argument("--obsolete_to_new_dict", type=str, default='Go_annotation/result/temp/obsolete_to_new.json', help='path of dictionary to replace the obsolete go id')
    parser.add_argument("--augment_go_anno", type=str, default="", help='the path to augment the go annotation to all related nodes, (Go_annotation/result/temp/augment_go_dict.json)')
    args = parser.parse_args()

    if not os.path.exists(args.result_folder):
        os.makedirs(args.result_folder)

    generate_gene_govec_mp(args=args, noIEA = True)
    generate_gene_govec_mp(args=args, noIEA = False)
    
if __name__ =="__main__":
    main()
    
    
# python get_gene_embedding_mp.py \
# --go_id_path 'Node_embedding_mxbai/go_id.npy' \
# --result_folder 'Gene_emb_simple_gaf' \
# --go_embedding_path 'Node_embedding_mxbai/47727.npy' \
# --method mean \
# --DBID_EntrezID_fulldict_path 'Go_annotation/result/temp/full_annotation_dict.json' \
# --ratio_cpu 0.25 \
# --evidence_level '1, 2' \
# --object_type 'protein, DNA, RNA, others' \
# --qualifier 'pos, neg'

# python get_gene_embedding_mp.py \
# --go_id_path 'Node_embedding_mxbai/go_id.npy' \
# --result_folder 'Gene_emb_simple_gaf' \
# --go_embedding_path 'Node_embedding_mxbai/47727.npy' \
# --method mean \
# --DBID_EntrezID_fulldict_path 'Go_annotation/result/temp/full_annotation_dict.json' \
# --ratio_cpu 0.7 \
# --evidence_level '-1, 0, 1, 2' \
# --object_type 'protein, DNA, RNA, others' \
# --qualifier 'pos'