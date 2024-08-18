from utils import *
import pandas as pd
import multiprocessing as mp

processed_columns = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'Aspect', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By']
columns17 = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'Qualifier', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'With_or_From', 'Aspect', 'DB_Object_Name', 'DB_Object_Synonym', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By', 'Annotation_Extension', 'Gene_Product_Form_ID']

gene_association = pd.read_csv("Go_annotation/extra_go/SGD/annotation_files/gene_association.sgd.gaf", sep='\t', header=None)
gene_association.columns = columns17
sgd = pd.read_csv("Go_annotation/extra_go/SGD/annotation_files/sgd.gaf", sep='\t', header=None)
sgd.columns = columns17




combine = pd.concat([gene_association[processed_columns], sgd[processed_columns]])
combine = combine.drop_duplicates(keep='last')
combine.to_csv("Go_annotation/extra_go/SGD/combine.gaf", sep='\t', header=True, index=False)



print(len(combine))
combine = pd.read_csv("Go_annotation/extra_go/SGD/combine.gaf", sep='\t')
print(set(combine['DB']))

db_to_id = {'UniProtKB':'ID_mapping/finish/UniProtKB_entrez.json', 
            'CGD': 'ID_mapping/finish/CGD_entrez.json',
            'EuPathDB': 'ID_mapping/finish/EupathDB_entrez.json', 'VEuPathDB': "ID_mapping/finish/EupathDB_entrez.json", 
            'dictyBase': "ID_mapping/finish/dictybase_entrez.json", 
            "FB": "ID_mapping/finish/FB_entrez.json", 
            "RefSeq": "ID_mapping/finish/refseq_entrez.json",
            "PR": "ID_mapping/finish/not_implement.json",
            "MGI": "ID_mapping/finish/mgi_entrez.json",
            "JaponicusDB": "ID_mapping/finish/japonicusdb_entrez.json",
            "PomBase": "ID_mapping/finish/pombase_entrez.json",
            "PseudoCAP": "ID_mapping/finish/pseudoCAP_entrez.json",
            "SGD": "ID_mapping/finish/sgd_entrez.json"}

mapping_dict = {db: load_json(db_to_id[db]) for db in set(combine['DB'])}

def process_row(index):
    print(index)
    global mapping_dict
    obj = combine.iloc[index]
    if obj['DB_Object_ID'] in mapping_dict[obj['DB']]:
        return str(mapping_dict[obj['DB']][obj['DB_Object_ID']])
    else:
        return '-1'

args = range(len(combine))
print('start mapping')
with mp.Pool(mp.cpu_count()) as pool:
    entrez_list = list(pool.imap(process_row, args))
    
# entrez_list = []  
# for i in tqdm(range(len(combine))):
#     obj = combine.iloc[i]
#     if obj['DB_Object_ID'] in mapping_dict[obj['DB']].keys():
#         entrez_list.append(str(mapping_dict[obj['DB']][obj['DB_Object_ID']]))
#     else:
#         entrez_list.append('-1')
unientrez = combine.copy()
unientrez['EntrezID'] = entrez_list
unientrez = unientrez[unientrez['EntrezID']!= '-1']
print('\n')
print(len(unientrez))
if not os.path.exists("Go_annotation/extra_go/SGD/unientrez_files"):
    os.makedirs("Go_annotation/extra_go/SGD/unientrez_files")
unientrez.to_csv("Go_annotation/extra_go/SGD/unientrez_files/uni_all.csv", sep='\t', header=True, index=False)
