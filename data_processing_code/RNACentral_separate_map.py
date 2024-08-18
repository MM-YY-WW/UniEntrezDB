from utils import *
import pandas as pd
import multiprocessing as mp

processed_columns = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'Aspect', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By']
columns17 = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'Qualifier', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'With_or_From', 'Aspect', 'DB_Object_Name', 'DB_Object_Synonym', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By', 'Annotation_Extension', 'Gene_Product_Form_ID']

# combine = pd.DataFrame()
# for filenames in tqdm(os.listdir("Go_annotation/extra_go/RNACentral/annotation_files")):
#     anno_files = pd.read_csv(os.path.join("Go_annotation/extra_go/RNACentral/annotation_files", filenames), sep='\t', comment='!', header=None, names =columns17)
#     combine = pd.concat([combine, anno_files[processed_columns]], axis=0)
# combine = combine.drop_duplicates(keep='last')
# combine.to_csv("Go_annotation/extra_go/RNACentral/combine.gaf", sep='\t', header=True, index=False)



#combine = pd.read_csv("Go_annotation/extra_go/RNACentral/combine.gaf", sep='\t', usecols=['DB'])
#print(set(combine['DB']))
#print(len(combine))
#combine = pd.read_csv("Go_annotation/extra_go/RNACentral/combine.gaf", sep='\t')
#combine = combine[combine['DB'] !='PDB']

db_to_id = {'UniProtKB':'ID_mapping/finish/UniProtKB_entrez.json', 
            'CGD': 'ID_mapping/finish/CGD_entrez.json',
            'EuPathDB': 'ID_mapping/finish/EupathDB_entrez.json', 'VEuPathDB': "ID_mapping/finish/EupathDB_entrez.json", 
            'dictyBase': "ID_mapping/finish/dictybase_entrez.json", 
            "FB": "ID_mapping/finish/FB_entrez.json", 
            "RefSeq": "ID_mapping/finish/refseq_entrez.json",
            "PR": "ID_mapping/finish/not_implement.json","ComplexPortal": "ID_mapping/finish/not_implement.json", "PDB": "ID_mapping/finish/not_implement.json",
            "MGI": "ID_mapping/finish/mgi_entrez.json",
            "JaponicusDB": "ID_mapping/finish/japonicusdb_entrez.json",
            "RNAcentral": "ID_mapping/finish/RNACentral_entraz.json"}

#mapping_dict = {db: load_json(db_to_id[db]) for db in {'UniProtKB', 'PDB', 'RNAcentral', 'ComplexPortal'}}

def process_row(index):
    print(index)
    global mapping_dict
    obj = combine.iloc[index]
    if obj['DB_Object_ID'] in mapping_dict[obj['DB']]:
        return str(mapping_dict[obj['DB']][obj['DB_Object_ID']])
    else:
        return '-1'

# for filename in tqdm(os.listdir("Go_annotation/extra_go/RNACentral/processed_annotation_files")):
#     print(filename)
#     if filename == "goa_uniprot_gcrp.gaf":
#         combine  = pd.read_csv(os.path.join("Go_annotation/extra_go/RNACentral/processed_annotation_files", filename), sep='\t', comment="!", usecols=['DB', 'DB_Object_ID'])
#         args = range(len(combine))
#         with mp.Pool(mp.cpu_count()) as pool:
#             entrez_list = list(pool.imap(process_row, args))
        
#         # entrez_list = []  
#         # for i in tqdm(range(len(combine))):
#         #     obj = combine.iloc[i]
#         #     if obj['DB_Object_ID'] in mapping_dict[obj['DB']].keys():
#         #         entrez_list.append(str(mapping_dict[obj['DB']][obj['DB_Object_ID']]))
#         #     else:
#         #         entrez_list.append('-1')
#         unientrez = combine.copy()
#         pd.Series(entrez_list).to_csv("Go_annotation/extra_go/RNACentral/goa_uniprot_gcrp_names", sep='\t', header=True, index=False)
#         unientrez['EntrezID'] = entrez_list
#         #unientrez = unientrez[unientrez['EntrezID']!= '-1']
#         #unientrez = unientrez.dropna(subset=['DB'])
#         #corrected_file = corrected_file[corrected_file['Aspect'].isin(['C', 'F', 'P'])]
#         # unientrez = unientrez[unientrez['Aspect'].isin(['C', 'F', 'P'])]
#         # unientrez = unientrez[unientrez['Evidence_Code'].isin(['TAS', 'ND', 'IC', 'NAS', 'IEA','ISO', 'IGC', 'ISM', 'ISS', 'RCA', 'ISA','IGI', 'HDA', 'EXP', 'IMP', 'IKR', 'IPI', 'HGI', 'HEP', 'IEP', 'IDA', 'HTP', 'IBA', 'HMP'])]
#         # unientrez = unientrez[unientrez['GO_ID'].str.contains('GO:')]
#         # unientrez = unientrez[unientrez['Taxon'].str.contains('taxon')]
#         print('\n')
#         #print(len(combine),len(unientrez), len(unientrez['Evidence_Code'] == 'IEA'))
#         print(len(combine),len(unientrez))
#         if not os.path.exists("Go_annotation/extra_go/RNACentral/unientrez_files"):
#             os.makedirs("Go_annotation/extra_go/RNACentral/unientrez_files")
#         unientrez.to_csv(f"Go_annotation/extra_go/RNACentral/unientrez_files/{len(unientrez)}_{len(combine)}_{filename[:-4]}.csv", sep='\t', header=True, index=False)

combine = pd.DataFrame()
for ufilename in tqdm(os.listdir("Go_annotation/extra_go/RNACentral/unientrez_files")):
#     check_col = pd.read_csv(os.path.join("Go_annotation/extra_go/uniprot/unientrez_files", ufilename), sep='\t', comment='!', nrows=2, usecols=processed_columns)
#     all_col[ufilename] = len(check_col.columns)
# print(Counter(all_col.values()))
    if ufilename == "40320744_377980902_goa_uniprot_gcrp.csv":
        entrez_ids = list(pd.read_csv("Go_annotation/extra_go/RNACentral/goa_uniprot_gcrp_names.csv", sep='\t',comment='!')['0'])
        print(len(entrez_ids))
        fulluniprot = pd.read_csv("Go_annotation/extra_go/RNACentral/processed_annotation_files/goa_uniprot_gcrp.gaf", sep='\t', comment='!')
        print(len(fulluniprot), len(entrez_ids)==len(fulluniprot))
        assert len(entrez_ids) == len(fulluniprot)
        print('finish read')
        fulluniprot['EntrezID'] = entrez_ids
        fulluniprot = fulluniprot[fulluniprot['EntrezID']!= '-1']
        fulluniprot = fulluniprot.dropna(subset=['DB'])
        print('start combine')
        combine = pd.concat([combine, fulluniprot], axis = 0)
    elif str(ufilename.split("_")[1]) != '0':
        unifile = pd.read_csv(os.path.join("Go_annotation/extra_go/RNACentral/unientrez_files", ufilename), sep='\t', comment='!', usecols=processed_columns)
        combine = pd.concat([combine, unifile], axis=0)
combine = combine[combine['EntrezID']!= '-1']
combine = combine.dropna(subset=['DB'])
combine = combine[combine['Aspect'].isin(['C', 'F', 'P'])]
combine = combine[combine['Evidence_Code'].isin(['TAS', 'ND', 'IC', 'NAS', 'IEA','ISO', 'IGC', 'ISM', 'ISS', 'RCA', 'ISA','IGI', 'HDA', 'EXP', 'IMP', 'IKR', 'IPI', 'HGI', 'HEP', 'IEP', 'IDA', 'HTP', 'IBA', 'HMP'])]
combine = combine[combine['GO_ID'].str.contains('GO:')]
combine = combine[combine['Taxon'].str.contains('taxon')]
combine.to_csv("Go_annotation/extra_go/RNACentral/unientrez_files/uni_all.csv", sep='\t', header=True, index=False)