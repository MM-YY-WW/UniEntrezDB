from utils import *
import pandas as pd


processed_columns = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'Aspect', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By']
file = pd.read_csv("Go_annotation/extra_go/refseq/annotation_files/GCF_027475565.1-RS_2023_04_gene_ontology.gaf", sep='\t',skiprows=8)
combine = file[['!#DB', 'GeneID', 'Symbol', 'GO_ID', 'Reference', "Evidence_Code", "Aspect", "Type", "Taxon", "Date", "Assigned_By"]]
combine.columns = processed_columns
combine['EntrezID'] = list(combine['DB_Object_ID'])
combine.to_csv("Go_annotation/extra_go/refseq/unientrez_files/uni_all.csv", sep='\t', header=True, index=False)
print(len(combine))
