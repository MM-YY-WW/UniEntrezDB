from utils import *
import multiprocessing as mp

# combine dictionary from eupathdb
# json_list = ['ID_mapping/finish/AmoebaDB_entrez.json', "ID_mapping/finish/CrytoDB_entrez.json",'ID_mapping/finish/FungiDB_entrez.json','ID_mapping/finish/GiardiaDB_entrez.json','ID_mapping/finish/ToxoDB_entrez.json','ID_mapping/finish/plasmoDB_entrez.json',"ID_mapping/finish/TriTryDB_entrez.json"]

# all_json  = {}
# for jl in tqdm(json_list):
#     json_file = load_json(jl)
#     all_json = {**all_json, **json_file}
#     print(len(all_json))
    
# save_json(all_json, "ID_mapping/finish/EupathDB_entrez.json", "w")
# a=1
# find fungi DB have wrong rows and corrected 
# filename_wrong = []
# for filename in tqdm(os.listdir("Go_annotation/extra_go/FungiDB/processed_annotation_files")):
#     file = pd.read_csv(f"Go_annotation/extra_go/FungiDB/processed_annotation_files/{filename}", sep='\t', comment='!', usecols=['Aspect'])
#     if len(Counter(file['Aspect'])) != 3:
#         corrected_file = pd.read_csv(f"Go_annotation/extra_go/FungiDB/processed_annotation_files/{filename}", sep='\t', comment='!')
#         corrected_file.to_csv(f"Go_annotation/extra_go/FungiDB/processed_annotation_files/{filename}", sep='\t', header=True, index=False)
#         filename_wrong.append(filename)

# print(filename_wrong)
# print(len(filename_wrong))
# processed_columns = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'Aspect', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By']

#combine = pd.read_csv("Go_annotation/extra_go/FungiDB/combined.gaf", sep='\t')
combine = pd.read_csv("Go_annotation/extra_go/FungiDB/unientrez_files/uni_all.csv", sep='\t')
combine = combine.dropna(subset=['DB'])
#corrected_file = corrected_file[corrected_file['Aspect'].isin(['C', 'F', 'P'])]
combine = combine[combine['Aspect'].isin(['C', 'F', 'P'])]
combine = combine[combine['Evidence_Code'].isin(['TAS', 'ND', 'IC', 'NAS', 'IEA','ISO', 'IGC', 'ISM', 'ISS', 'RCA', 'ISA','IGI', 'HDA', 'EXP', 'IMP', 'IKR', 'IPI', 'HGI', 'HEP', 'IEP', 'IDA', 'HTP', 'IBA', 'HMP'])]
combine = combine[combine['GO_ID'].str.contains('GO:')]
combine = combine[combine['Taxon'].str.contains('taxon')]

# combine.to_csv("Go_annotation/extra_go/FungiDB/combined.gaf", sep='\t', header=True, index=False)
# combine.to_csv("Go_annotation/extra_go/FungiDB/combined.gaf", sep='\t', header=True, index=False)
print(len(combine))
print(set(combine['DB']))
print(set(combine['Aspect']))
combine.to_csv("Go_annotation/extra_go/FungiDB/unientrez_files/uni_all.csv", sep='\t', header=True, index=False)

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
            "Xenbase": "ID_mapping/finish/Xenbase_entrez.json",
            "ZFIN": "ID_mapping/finish/ZFIN_entrez.json"}

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
print(len(combine),len(unientrez))
if not os.path.exists("Go_annotation/extra_go/FungiDB/unientrez_files"):
    os.makedirs("Go_annotation/extra_go/FungiDB/unientrez_files")
unientrez.to_csv("Go_annotation/extra_go/FungiDB/unientrez_files/uni_all.csv", sep='\t', header=True, index=False)
