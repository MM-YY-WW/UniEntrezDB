from utils import *
import pandas as pd
import multiprocessing as mp


#combine_all = pd.DataFrame()
# columns = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'Aspect', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By']
# columns16 = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'Qualifier', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'With_or_From', 'Aspect', 'DB_Object_Name', 'DB_Object_Synonym', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By', 'Annotation_Extension', 'Gene_Product_Form_ID']

# entrez_map = load_json("ID_mapping/finish/AmoebaDB_entrez.json")
# annotation_folder = 'Go_annotation/extra_go/AmoebaDB/annotaion_files'
# length_track = {}
# evidence_level_symbol = ['TAS', 'ND', 'IC', 'NAS', 'IEA','ISO', 'IGC', 'ISM', 'ISS', 'RCA', 'ISA','IGI', 'HDA', 'EXP', 'IMP', 'IKR', 'IPI', 'HGI', 'HEP', 'IEP', 'IDA', 'HTP', 'IBA', 'HMP']
# for filename in tqdm(os.listdir(annotation_folder)):
#     if '.gaf' in filename:
#         anno_file = pd.read_csv(os.path.join("Go_annotation/extra_go/AmoebaDB/annotaion_files", filename), sep='\t', header=None, skiprows=1)
#         anno_file = anno_file.dropna(axis=1, how='any')
#         length = len(anno_file.iloc[0])
#         if length == 12:
#             anno_file.columns = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'Aspect', 'DB_Object_Name', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By']
#         if length == 13:
#             anno_file.columns = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'With_or_From', 'Aspect', 'DB_Object_Name', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By']
#         length_track[filename] = length
#         anno_file.to_csv(os.path.join("Go_annotation/extra_go/AmoebaDB/nonan_annotation_files", f'nonan_{filename}'), header=True, index=False, sep='\t')

# print(Counter(length_track.values()))
# save_json(length_track, "Go_annotation/extra_go/AmoebaDB/length_track.json", "w")

# for filename in tqdm(os.listdir("Go_annotation/extra_go/AmoebaDB/nonan_annotation_files")):
#     if '.gaf' in filename:
#         anno_file = pd.read_csv(os.path.join("Go_annotation/extra_go/AmoebaDB/nonan_annotation_files", filename), sep='\t', nrows=1)
#         to_combine = anno_file[['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'Aspect', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By']]
#         to_combine['filename'] = filename
#         print(filename)
#         print(to_combine)
#         combine_all = pd.concat([combine_all, to_combine], axis=0, ignore_index=True)
# combine_all.to_csv("Go_annotation/extra_go/AmoebaDB/combined_goa_sample.gaf", index=False, sep='\t')

# DB	DB_Object_ID	DB_Object_Symbol	GO_ID	DB:Reference	Evidence_Code	Aspect	DB_Object_Name	DB_Object_Type	Taxon	Date	Assigned_By
# VEuPathDB	NF0000010	NF0000010	GO:0000160	GO_REF:0000002	IEA	P	unspecified product	protein_coding	taxon:9000000014	20231012	interpro2go

# /home/yuwei/projects_backup/Gene_Ontology_PT/Go_annotation/extra_go/AmoebaDB/nonan_annotation_files/nonan_AmoebaDB-66_NfowleriATCC30863_GO.gaf

        # db =0
        # db_object_id =1
        # db_object_symbol=2
        # go_id = -1
        # reference = -1
        # evidence= -1
        # aspect = -1
        # taxon = -1
        # for i in range(anno_file.iloc[1]):
        #     if 'GO:' in anno_file.iloc[1][i]:
        #         go_id = i
        #     if "REF" in anno_file.iloc[1][i]:
        #         reference = i
        #     if anno_file.iloc[1][i] in evidence_level_symbol:
        #         evidence = i
        #     if anno_file.iloc[1][i] in ['F', 'P', 'C']:
        #         aspect = i
        #     if 'taxon' in anno_file.iloc[1][i]:
        #         taxon = i 
        # # anno_file.columns = columns
        # combine_all = pd.concat([combine_all, anno_file], axis=0)
    
#EuPathDB	EHI8A_106390	EHI8A_106390		GO:0005515	GO_REF:0000002	IEA		F	protein kinase, putative		transcript	taxon:885319	20180720	EuPathDB	


	

processed_columns = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'Aspect', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By']
columns16 = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'Qualifier', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'With_or_From', 'Aspect', 'DB_Object_Name', 'DB_Object_Synonym', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By', 'Annotation_Extension', 'Gene_Product_Form_ID']




# combine = pd.read_csv("Go_annotation/extra_go/AmoebaDB/combined.gaf", comment='!', names=columns16, sep='\t')
# combine = combine[processed_columns]
# combine = combine.drop_duplicates(keep='last')
# combine.to_csv("Go_annotation/extra_go/AmoebaDB/combine.gaf", sep='\t', header=True, index=False)



combine = pd.read_csv("Go_annotation/extra_go/AmoebaDB/combine.gaf", sep='\t')
print(len(combine))
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
print(len(combine),len(unientrez), len(unientrez['Evidence_Code'] == 'IEA'))
if not os.path.exists("Go_annotation/extra_go/AmoebaDB/unientrez_files"):
    os.makedirs("Go_annotation/extra_go/AmoebaDB/unientrez_files")
unientrez.to_csv("Go_annotation/extra_go/AmoebaDB/unientrez_files/uni_all.csv", sep='\t', header=True, index=False)
