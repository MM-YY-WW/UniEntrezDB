from utils import *
import pandas as pd
import multiprocessing as mp

#===============Unify Columns

processed_columns = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'Aspect', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By']

# columns16 = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'Qualifier', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'With_or_From', 'Aspect', 'DB_Object_Name', 'DB_Object_Synonym', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By', 'Annotation_Extension', 'Gene_Product_Form_ID']

# cgd_from_cgd = pd.read_csv("Go_annotation/extra_go/CGD/annotation_files/cgd_from_cgd.gaf", sep='\t', header=None, skiprows=49)
# cgd_from_cgd.columns = ['DB','DB_Object_ID','DB_Object_Symbol', 'Qualifier', 'GO_ID', 'DB:Reference','Evidence_Code','With_or_From','Aspect','DB_Object_Name', 'DB_Object_Synonym', 'DB_Object_Type', 'Taxon','Date','Assigned_By','Annotation_Extension', 'Gene_Product_Form_ID']
# processed_cgd_from_cgd = cgd_from_cgd[processed_columns]
# processed_cgd_from_cgd.to_csv("Go_annotation/extra_go/CGD/processed_annotation_files/processed_cgd_from_cgd.gaf", header=True, index=False, sep='\t')

# cgd = pd.read_csv("Go_annotation/extra_go/CGD/annotation_files/cgd.gaf", sep='\t', header=None, skiprows=56)
# cgd.columns = ['DB','DB_Object_ID','DB_Object_Symbol', 'Qualifier', 'GO_ID', 'DB:Reference','Evidence_Code','With_or_From','Aspect','DB_Object_Name', 'DB_Object_Synonym', 'DB_Object_Type', 'Taxon','Date','Assigned_By','Annotation_Extension', 'Gene_Product_Form_ID']
# processed_cgd = cgd[processed_columns]
# processed_cgd.to_csv("Go_annotation/extra_go/CGD/processed_annotation_files/processed_cgd.gaf", header=True, index=False, sep='\t')


# gene_association = pd.read_csv("Go_annotation/extra_go/CGD/annotation_files/gene_association.cgd", sep='\t', header=None, skiprows=21)
# gene_association.columns = ['DB','DB_Object_ID','DB_Object_Symbol', 'Qualifier', 'GO_ID', 'DB:Reference','Evidence_Code','With_or_From','Aspect','DB_Object_Name', 'DB_Object_Synonym', 'DB_Object_Type', 'Taxon','Date','Assigned_By','Annotation_Extension', 'Gene_Product_Form_ID']
# processed_gene_association = gene_association[processed_columns]
# processed_gene_association.to_csv("Go_annotation/extra_go/CGD/processed_annotation_files/processed_gene_association.cgd", header=True, index=False, sep='\t')

# GOslim_gene_association = pd.read_csv("Go_annotation/extra_go/CGD/annotation_files/GOslim_gene_association.cgd", sep='\t', header=None, skiprows=11)
# GOslim_gene_association.columns = ['DB','DB_Object_ID','DB_Object_Symbol', 'Qualifier', 'GO_ID', 'DB:Reference','Evidence_Code','With_or_From','Aspect','DB_Object_Name', 'DB_Object_Synonym', 'DB_Object_Type', 'Taxon','Date','Assigned_By']
# processed_GOslim_gene_association = GOslim_gene_association[processed_columns]
# processed_GOslim_gene_association.to_csv("Go_annotation/extra_go/CGD/processed_annotation_files/processed_GOslim_gene_association.cgd", header=True, index=False, sep='\t')


#===============Combine

# combine = pd.DataFrame()
# for filename in tqdm(os.listdir("Go_annotation/extra_go/CGD/processed_annotation_files")):
#     anno_file = pd.read_csv(os.path.join("Go_annotation/extra_go/CGD/processed_annotation_files", filename), sep='\t')
#     combine = pd.concat([combine, anno_file], axis = 0)
#combine.to_csv("Go_annotation/extra_go/CGD/combined.gaf", sep='\t', header=True, index=False)


#=================UniEntrez
def process_row(index):
    print(index)
    global mapping_dict
    obj = combine.iloc[index]
    if obj['DB_Object_ID'] in mapping_dict[obj['DB']]:
        return str(mapping_dict[obj['DB']][obj['DB_Object_ID']])
    else:
        return '-1'

combine = pd.read_csv("Go_annotation/extra_go/CGD/combined.gaf", sep='\t')
print(set(combine['DB']))
db_to_id = {'UniProtKB':'ID_mapping/finish/UniProtKB_entrez.json', 'CGD': 'ID_mapping/finish/CGD_entrez.json'}
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
unientrez.to_csv("Go_annotation/extra_go/CGD/unientrez_files/uni_all.csv", sep='\t', header=True, index=False)

