import pandas as pd
# from utils import *
from collections import Counter
import os
processed_columns = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'Aspect', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By']
processed_columns_entrez = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'Aspect', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By', 'EntrezID']

# dict_keys=['pombase', 'SGD', 'TAIR', 'CGD', 'MGI', 'Xenbase', 'ZFIN', 'PseudoCAP', 'dictybase', 
#            'refseq', 'PlasmoDB', 'JaponicusDB', 'AmoebaDB', 'CrytoDB', 'FungiDB', 'GiardiaDB', 'ToxoDB', 'FB']
all_database = pd.read_csv("Go_annotation/extra_go/uni_all_database.csv", sep='\t', )
statistic_dict = load_json("Go_annotation/extra_go/statistics/statistics_separate.json")
nIEA= {k: statistic_dict[k]['noIEA']for k in statistic_dict.keys()}
files = []
extra_go_folder = "Go_annotation/extra_go"
statistic_dict  = {}
for folder in tqdm(os.listdir(extra_go_folder)):
    uni_file_path = f"{extra_go_folder}/{folder}/unientrez_files/uni_all.csv"
    if os.path.exists(uni_file_path):
        print(uni_file_path)
        # uni_file = pd.read_csv(uni_file_path, sep='\t', usecols=['DB', 'Date','Evidence_Code','Aspect','Taxon','DB_Object_Type'])
        uni_file = pd.read_csv(uni_file_path, sep='\t')
        files.append(uni_file)
        statistic_dict[folder] = {'total_length' :len(uni_file),
                        'DB': Counter(uni_file['DB']),
                        #'Date': Counter(uni_file['Date']),
                        'Evidence_Code': Counter(uni_file['Evidence_Code']),
                        'Aspect': Counter(uni_file['Aspect']),
                        'Taxon': Counter(uni_file['Taxon']),
                        #'noIEA_Date': Counter(uni_file[uni_file['Evidence_Code'] != 'IEA']['Date']),
                        'noIEA_Aspect': Counter(uni_file[uni_file['Evidence_Code'] != 'IEA']['Aspect']),
                        'noIEA_Taxon': Counter(uni_file[uni_file['Evidence_Code'] != 'IEA']['Taxon']),
                        'noIEA' : len(uni_file[uni_file['Evidence_Code'] != 'IEA']),
                        "DB_Object_Type" : Counter(uni_file['DB_Object_Type']),
                        "noIEA_DB_Object_Type" : Counter(uni_file[uni_file['Evidence_Code'] != 'IEA']['DB_Object_Type'])
                        }
    save_json(statistic_dict, "Go_annotation/extra_go/statistics/statistics_separate.json", "w")
   
combine = pd.concat(files)
combine = combine.drop_duplicates()
process_combine = combine[processed_columns_entrez]
process_combine = process_combine[process_combine['Aspect'].isin(['C', 'F', 'P'])]
process_combine = process_combine[process_combine['Evidence_Code'].isin(['TAS', 'ND', 'IC', 'NAS', 'IEA','ISO', 'IGC', 'ISM', 'ISS', 'RCA', 'ISA','IGI', 'HDA', 'EXP', 'IMP', 'IKR', 'IPI', 'HGI', 'HEP', 'IEP', 'IDA', 'HTP', 'IBA', 'HMP'])]
process_combine = process_combine[process_combine['GO_ID'].str.contains('GO:')]
process_combine = process_combine[process_combine['Taxon'].str.contains('taxon')]
unifided_entrez = []
for i in tqdm(process_combine['EntrezID']):
    if not ";" in str(i):
        unifided_entrez.append(str(int(float(i))))
    else:
        unifided_entrez.append(i)
process_combine['EntrezID'] = unifided_entrez



combine.to_csv("Go_annotation/extra_go/uni_all_database.csv", sep='\t', header=True, index=False)
    

a=1