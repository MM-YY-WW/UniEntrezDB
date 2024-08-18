from utils import *
import pandas as pd
import multiprocessing as mp


processed_columns = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'Aspect', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By']


#===============Combine

# combine = pd.DataFrame()
# for filename in tqdm(os.listdir("Go_annotation/extra_go/CrytoDB/processed_annotation_files")):
#     anno_file = pd.read_csv(os.path.join("Go_annotation/extra_go/CrytoDB/processed_annotation_files", filename), sep='\t')
#     combine = pd.concat([combine, anno_file], axis = 0)
# combine.to_csv("Go_annotation/extra_go/CrytoDB/combined.gaf", sep='\t', header=True, index=False)


#=================UniEntrez
def process_row(index):
    print(index)
    global mapping_dict
    obj = combine.iloc[index]
    if obj['DB_Object_ID'] in mapping_dict[obj['DB']]:
        return str(mapping_dict[obj['DB']][obj['DB_Object_ID']])
    else:
        return '-1'

combine = pd.read_csv("Go_annotation/extra_go/CrytoDB/combined.gaf", sep='\t')
print(set(combine['DB']))
db_to_id = {'UniProtKB':'ID_mapping/finish/UniProtKB_entrez.json', 'CGD': 'ID_mapping/finish/CGD_entrez.json','EuPathDB': 'ID_mapping/finish/EupathDB_entrez.json', 'VEuPathDB': "ID_mapping/finish/EupathDB_entrez.json"}
mapping_dict = {db: load_json(db_to_id[db]) for db in set(combine['DB'])}

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
unientrez.to_csv("Go_annotation/extra_go/CrytoDB/unientrez_files/uni_all.csv", sep='\t', header=True, index=False)
