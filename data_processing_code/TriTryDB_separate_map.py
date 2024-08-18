from utils import *
import pandas as pd
import multiprocessing as mp
from collections import Counter 

processed_columns = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'Aspect', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By']
columns16 = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'Qualifier', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'With_or_From', 'Aspect', 'DB_Object_Name', 'DB_Object_Synonym', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By', 'Annotation_Extension', 'Gene_Product_Form_ID']




# combine = pd.read_csv("Go_annotation/extra_go/AmoebaDB/combined.gaf", comment='!', names=columns16, sep='\t')
# combine = combine[processed_columns]
# combine = combine.drop_duplicates(keep='last')
# combine.to_csv("Go_annotation/extra_go/AmoebaDB/combine.gaf", sep='\t', header=True, index=False)



combine = pd.read_csv("Go_annotation/extra_go/TriTryDB/combine.gaf", sep='\t')
print(len(combine))
print(Counter(combine['DB']))

db_to_id = {'UniProtKB':'ID_mapping/finish/UniProtKB_entrez.json', 
            'CGD': 'ID_mapping/finish/CGD_entrez.json',
            'EuPathDB': 'ID_mapping/finish/EupathDB_entrez.json', 'VEuPathDB': "ID_mapping/finish/EupathDB_entrez.json", 
            'dictyBase': "ID_mapping/finish/dictybase_entrez.json", 
            "FB": "ID_mapping/finish/FB_entrez.json", 
            "RefSeq": "ID_mapping/finish/refseq_entrez.json",
            "PR": "ID_mapping/finish/not_implement.json", "GeneDB": "ID_mapping/finish/not_implement.json",
            "MGI": "ID_mapping/finish/mgi_entrez.json",
            "JaponicusDB": "ID_mapping/finish/japonicusdb_entrez.json",
            "PomBase": "ID_mapping/finish/pombase_entrez.json",
            "Xenbase": "ID_mapping/finish/Xenbase_entrez.json",
            "ZFIN": "ID_mapping/finish/ZFIN_entrez.json",
            "TriTrypDB": "ID_mapping/finish/TriTryDB_entrez.json"}


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
print(len(combine),len(unientrez), len(unientrez['Evidence_Code'] == 'IEA'))
if not os.path.exists("Go_annotation/extra_go/TriTryDB/unientrez_files"):
    os.makedirs("Go_annotation/extra_go/TriTryDB/unientrez_files")
unientrez.to_csv("Go_annotation/extra_go/TriTryDB/unientrez_files/uni_all.csv", sep='\t', header=True, index=False)
