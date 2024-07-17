# UnientrezDB: Large-scale Gene Ontology Annotation Dataset and Evaluation Benchmarks with Unified Entrez Gene Identifiers

This repository contains the official implementation for the paper "UnientrezDB: Large-scale Gene Ontology Annotation Dataset and Evaluation Benchmarks with Unified Entrez Gene Identifiers." Our work focuses on providing a comprehensive dataset and benchmarks for evaluating gene ontology annotations using a unified system of Entrez Gene Identifiers.

## UniEntrez GOA:

Pretrain dataset: Maps GO annotations of gene and gene products to Gene Entrez ID. 

**Pretrained Dataset**: [Download Pretrained Dataset](https://drive.google.com/file/d/1DsXufybeSgEXrx8szkF0kuhASmAVOaU-/view?usp=sharing)

Contains the following 12 columns:

- DB
- DB Object ID
- DB Object Symbol
- GO ID
- DB:Reference(|DB:Reference)
- Evidence Code
- Aspect
- DB Object Type
- Taxon(|taxon)
- Date
- Assigned By
- EntrezID

![Figure 1](https://github.com/MM-YY-WW/UniEntrezDB/blob/main/Figures/goa.png)


## UniEntrez Evaluation Benchmarks

Evaluation Benchmarks datasets:

**Downstream Task Datasets**: [Download Datasets](https://drive.google.com/file/d/1fSRXO26jr1XcFn7GKqRoN_CZUbuEY8Cj/view?usp=sharing)

- Pathway Co-present Prediction (Gene-level)
- Functional Gene Interaction Prediction (Gene-level)
- Protein-Protein Interaction (Protein-level)
- Single-Cell Type Annotation (Cell-level)

![Figure 2](https://github.com/MM-YY-WW/UniEntrezDB/blob/main/Figures/evaluation_benchmark.png)


## ID Mapping Between Public Biological Databases

**ID Mapping Dictionaries**: [Download ID Mapping dictionaries](https://drive.google.com/file/d/1La80B3hUibbe94FghkTIx80DRzPfwYix/view?usp=sharing)

![Figure 3](https://github.com/MM-YY-WW/UniEntrezDB/blob/main/Figures/id_mapping_to_Database.png)


## Reproduce our work
To run the downstream tasks, follow the steps at [overall_results.ipynb](overeall_results.ipynb)


### Embeddings Download

To download the embeddings generated in this study, use the following link.

- **Embeddings(GOA, Gene2Vec, DNABert, OntoProtein)**: [Download Embeddings](https://drive.google.com/file/d/1OcAnUT6CJEDsQk2hPlPE2tpf-hL9nDA4/view?usp=sharing)

### Our Results

![Figure 4](https://github.com/MM-YY-WW/UniEntrezDB/blob/main/Figures/Results.png)










