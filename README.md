# UnientrezDB: Large-scale Gene Ontology Annotation Dataset and Evaluation Benchmarks with Unified Entrez Gene Identifiers

This repository contains the official implementation for the paper "UnientrezDB: Large-scale Gene Ontology Annotation Dataset and Evaluation Benchmarks with Unified Entrez Gene Identifiers." Our work focuses on providing a comprehensive dataset and benchmarks for evaluating gene ontology annotations using a unified system of Entrez Gene Identifiers.

## Overall Average Results

#### [DISEASES](https://diseases.jensenlab.org/About)

| Methods | AP\% | APOP\% | AUCROC\% | 
| :------ | ---------: | -------------: | -------------: | 
| DNABert-2 | 5.10 $\pm$ 3.98 | 44.72 $\pm$ 9.93 | 52.11 $\pm$ 0.82 |
| Gene2Vec | 5.52 $\pm$ 4.00 | 55.18 $\pm$ 9.96 | 52.09 $\pm$ 0.85 |
| OntoProtein | 4.08 $\pm$ 3.49 | 11.95 $\pm$ 3.67 | 52.24 $\pm$ 1.00 |
| GOA_Emb | 5.32 $\pm$ 4.10 | 48.97 $\pm$ 11.39 | 51.97 $\pm$ 0.94 |
| GOA_Emb + Gene2Vec | 6.83 $\pm$ 3.46 | 81.72 $\pm$ 23.00 | 54.75 $\pm$ 1.29 |
| GOA_Emb + DNABert-2 | 4.16 $\pm$ 3.47 | 14.96 $\pm$ 11.94 | 52.20$\pm$ 0.98 |


#### [DisGeNET](https://www.disgenet.org/)

| Methods | AP\% | APOP\% | AUCROC\% | 
| :------ | ---------: | -------------: | -------------: | 
| DNABert-2 | 3.33 $\pm$ 1.53 | 64.01 $\pm$ 13.46 | 53.54 $\pm$ 1.21 |
| Gene2Vec | 2.59 $\pm$ 1.66 | 42.68 $\pm$ 5.12 | 52.15 $\pm$ 0.72 |
| OntoProtein | 2.01 $\pm$ 1.44 | 9.905 $\pm$ 1.61 | 51.86 $\pm$ 0.73 |
| GOA_Emb | 2.57 $\pm$ 1.65 | 40.81 $\pm$ 6.46 | 52.06 $\pm$ 0.78 |
| GOA_Emb + Gene2Vec | 3.12 $\pm$ 1.55 | 59.72 $\pm$ 10.59 | 53.41 $\pm$ 0.79 |
| GOA_Emb + DNABert-2 | 2.01 $\pm$ 1.44 | 9.91 $\pm$ 1.60 | 51.86 $\pm$ 0.73 |


### DNABert-2:

#### [DISEASES](https://diseases.jensenlab.org/About)

| Network | AP\% | APOP\% | AUCROC\% | 
| :------ | ---------: | -------------: | -------------: | 
| BioGRID | 3.18 | 35.43 | 51.19 |
| BioPlex | 6.57 | 52.60 | 51.25 |
| ComPPIHumanInt | 3.31 | 43.74 | 52.34 |
| ConsensusPathDB | 3.49 | 54.73 | 53.39 |
| FunCoup | 3.11 | 34.97 | 51.50 |
| HIPPIE | 3.21 | 40.71 | 52.23 |
| HumanNet | 3.66 | 55.98 | 52.26 |
| HuMAP | 4.17 | 54.49 | 52.31 |
| HuRI | 6.27 | 34.61 | 51.47 |
| OmniPath | 3.58 | 49.02 | 52.41 |
| PCNet | 3.28 | 43.08 | 52.05 |
| ProteomeHD | 18.91 | 23.85 | 51.48 |
| SIGNOR | 5.08 | 42.76 | 51.46 |
| STRING | 3.59 | 60.06 | 54.20 |

#### [DisGeNET](https://www.disgenet.org/)

| Network | AP\% | APOP\% | AUCROC\% | 
| :------ | ---------: | -------------: | -------------: | 
| BioGRID | 2.64 | 67.97 | 53.70 |
| BioPlex | 3.39 | 46.30 | 53.79 |
| ComPPIHumanInt | 2.88 | 73.77 | 54.30 |
| ConsensusPathDB | 2.63 | 65.54 | 53.80 |
| FunCoup | 2.77 | 69.70 | 54.15 |
| HIPPIE | 2.82 | 80.19 | 54.27 |
| HumanNet | 2.80 | 78.67 | 53.86 |
| HuMAP | 2.77 | 65.57 | 53.72 |
| HuRI | 3.41 | 51.86 | 51.80 |
| OmniPath | 3.06 | 78.76 | 54.68 |
| PCNet | 2.48 | 64.11 | 53.29 |
| ProteomeHD | 8.74 | 32.09 | 49.88 |
| SIGNOR | 3.42 | 50.97 | 54.20 |
| STRING | 2.76 | 70.70 | 54.17 |




### Gene2Vec:

#### [DISEASES](https://diseases.jensenlab.org/About)

| Network | AP\% | APOP\% | AUCROC\% | 
| :------ | ---------: | -------------: | -------------: | 
| BioGRID | 3.93|60.55|52.82 |
| BioPlex | 6.47|49.25|51.56 |
| ComPPIHumanInt | 3.78|59.31|52.58 |
| ConsensusPathDB | 3.87|58.21|52.52 |
| FunCoup | 3.89|63.59|52.07 |
| HIPPIE | 3.63|53.73|52.23 |
| HumanNet | 3.78|60.85|53.10 |
| HuMAP | 4.32|59.06|52.21 |
| HuRI | 6.38|38.57|50.58 |
| OmniPath | 4.2|64.95|53.41 |
| PCNet | 3.91|60.76|52.72 |
| ProteomeHD | 19.49|30.16|50.43 |
| SIGNOR | 6.15 | 65.30 | 51.58 |
| STRING | 3.54|48.26|51.49 |

#### [DisGeNET](https://www.disgenet.org/)

| Network | AP\% | APOP\% | AUCROC\% | 
| :------ | ---------: | -------------: | -------------: | 
| BioGRID | 1.81|45.33|52.67 |
| BioPlex | 3.13|39.96|51.62 |
| ComPPIHumanInt | 1.86|43.55|52.01 |
| ConsensusPathDB | 1.94|49.13|52.67 |
| FunCoup | 1.97|48.16|52.09 |
| HIPPIE | 1.68|38.18|51.82 |
| HumanNet | 1.67|37.88|51.90 |
| HuMAP | 1.93|42.75|52.52 |
| HuRI | 3.17|42.43|50.37 |
| OmniPath | 1.88|44.14|51.99 |
| PCNet | 1.71|43.7|53.21 |
| ProteomeHD | 8.15|29.16|53.47 |
| SIGNOR | 3.59 | 49.75 | 51.89 |
| STRING | 1.76|43.41|51.90 |



### OntoProtein:

#### [DISEASES](https://diseases.jensenlab.org/About)

| Network | AP\% | APOP\% | AUCROC\% | 
| :------ | ---------: | -------------: | -------------: | 
| BioGRID | 2.63|12.15|52.46 |
| BioPlex | 4.94|14.12|52.49 |
| ComPPIHumanInt | 2.65|13.93|52.93 |
| ConsensusPathDB | 2.62|13.75|52.81 |
| FunCoup | 2.58|13.5|52.75 |
| HIPPIE | 2.67|12.97|52.72 |
| HumanNet | 2.65|13.97|52.93 |
| HuMAP | 3.05|12.15|52.32 |
| HuRI | 5.11|7.08|51.15 |
| OmniPath | 2.76|11.37|52.29 |
| PCNet | 2.60 |11.64|52.290 |
| ProteomeHD | 16.29|1.01|49.00 |
| SIGNOR | 4.00 | 16.85 | 52.59 |
| STRING | 2.61|12.87|52.68 |

#### [DisGeNET](https://www.disgenet.org/)

| Network | AP\% | APOP\% | AUCROC\% | 
| :------ | ---------: | -------------: | -------------: | 
| BioGRID | 1.37|10.59|52.19 |
| BioPlex | 2.46|10.26|52.72 |
| ComPPIHumanInt | 1.45|10.72|52.36 |
| ConsensusPathDB | 1.36|7.29|51.39 |
| FunCoup | 1.45|11.95|52.36 |
| HIPPIE | 1.37|10.96|52.31 |
| HumanNet | 1.36|12.09|52.27 |
| HuMAP | 1.52|10.13|51.81 |
| HuRI | 2.23|7.81|51.02 |
| OmniPath | 1.43|11.05|51.84 |
| PCNet | 1.35|8.69|51.72 |
| ProteomeHD | 7.01|6.83|49.73 |
| SIGNOR | 2.44 | 9.36 | 52.08 |
| STRING | 1.37|10.95|52.26 |



### GOA_Emb:

#### [DISEASES](https://diseases.jensenlab.org/About)

| Network | AP\% | APOP\% | AUCROC\% | 
| :------ | ---------: | -------------: | -------------: | 
| BioGRID | 3.18|35.43|51.19 |
| BioPlex | 6.57|52.6|51.25 |
| ComPPIHumanInt | 3.31|43.74|52.34 |
| ConsensusPathDB |3.49|54.73|53.39 |
| FunCoup | 3.11|34.97|51.50 |
| HIPPIE | 3.21|40.71|52.23 |
| HumanNet | 3.78|60.85|53.10 |
| HuMAP | 4.17|54.49|52.31 |
| HuRI | 6.38|38.57|50.58 |
| OmniPath | 4.20|64.95|53.41 |
| PCNet | 3.91|60.76|52.72 |
| ProteomeHD | 19.49|30.16|50.43 |
| SIGNOR | 6.15 | 65.30 | 51.58 |
| STRING | 3.54|48.26|51.49 |

#### [DisGeNET](https://www.disgenet.org/)

| Network | AP\% | APOP\% | AUCROC\% | 
| :------ | ---------: | -------------: | -------------: | 
| BioGRID | 1.88|46.09|51.46 |
| BioPlex | 2.90|28.40|52.09 |
| ComPPIHumanInt | 1.8|34.72|51.88 |
| ConsensusPathDB | 1.89|49.36|52.37 |
| FunCoup | 1.8|35.68|51.83 |
| HIPPIE | 1.90|43.88|51.97 |
| HumanNet | 1.67|37.88|51.90 |
| HuMAP | 1.93|42.75|52.52 |
| HuRI | 3.17|42.43|50.37 |
| OmniPath | 1.88|44.14|51.99 |
| PCNet | 1.71|43.70|53.21 |
| ProteomeHD | 8.15|29.16|53.47 |
| SIGNOR | 3.59 | 49.75 | 51.89 |
| STRING | 1.76|43.41|51.90 |


### GOA_Emb + Gene2Vec 

#### [DISEASES](https://diseases.jensenlab.org/About)

| Network | AP\% | APOP\% | AUCROC\% | 
| :------ | ---------: | -------------: | -------------: | 
| BioGRID |6.15|103.22|56.26 |
| BioPlex | 6.67|50.37|53.02 |
| ComPPIHumanInt | 5.24|90.44|55.27 |
| ConsensusPathDB | 5.63|99.76|56.13 |
| FunCoup |5.28|96.98|55.94 |
| HIPPIE | 6.23|100.55|55.82 |
| HumanNet |5.51|93.34|54.92 |
| HuMAP | 5.29|67.07|53.53 |
| HuRI | 7.82|64.83|53.23 |
| OmniPath | 4.94|85.69|54.96 |
| PCNet | 6.13|102.7|55.18 |
| ProteomeHD | 19.03|21.35|52.24 |
| SIGNOR | 6.35 | 73.22 | 53.88 |
| STRING | 5.33|94.55|56.17 |

#### [DisGeNET](https://www.disgenet.org/)

| Network | AP\% | APOP\% | AUCROC\% | 
| :------ | ---------: | -------------: | -------------: | 
| BioGRID | 2.30|60.76|53.49 |
| BioPlex | 3.50|50.85|53.63 |
| ComPPIHumanInt | 2.49|67.72|53.61 |
| ConsensusPathDB | 2.47|68.76|54.10 |
| FunCoup | 2.58|60.81|53.11 |
| HIPPIE | 2.41|66.58|54.05 |
| HumanNet | 2.48|69.17|53.30 |
| HuMAP |2.59|63.02|54.01|
| HuRI | 3.35|45.55|51.73 |
| OmniPath | 2.59|68.18|54.43 |
| PCNet | 2.46|65.74|53.84 |
| ProteomeHD | 8.53|31.93|51.60 |
| SIGNOR | 3.38 | 51.04 | 53.26 |
| STRING | 2.51|65.93|53.60 |
### GOA_Emb + DNABert-2


#### [DISEASES](https://diseases.jensenlab.org/About)

| Network | AP\% | APOP\% | AUCROC\% | 
| :------ | ---------: | -------------: | -------------: | 
| BioGRID | 2.63|12.15|52.46 |
| BioPlex |4.94|14.12|52.49 |
| ComPPIHumanInt | 2.65|13.93|52.93|
| ConsensusPathDB |2.62|13.75|52.81 |
| FunCoup | 2.58|13.5|52.75 |
| HIPPIE | 2.67|12.97|52.72 |
| HumanNet | 3.66|55.98|52.26 |
| HuMAP | 3.05|12.15|52.32 |
| HuRI | 5.11|7.08|51.15 |
| OmniPath | 2.76|11.37|52.29 |
| PCNet |2.6|11.64|52.29 |
| ProteomeHD | 16.29|1.01|49.0 |
| SIGNOR | 5.08 | 42.76 | 51.46 |
| STRING | 2.61|12.87|52.68 |

#### [DisGeNET](https://www.disgenet.org/)

| Network | AP\% | APOP\% | AUCROC\% | 
| :------ | ---------: | -------------: | -------------: | 
| BioGRID | 1.37|10.59|52.19 |
| BioPlex |2.46|10.26|52.72 |
| ComPPIHumanInt | 1.45|10.72|52.36 |
| ConsensusPathDB | 1.36|7.29|51.39 |
| FunCoup | 1.45|11.95|52.36 |
| HIPPIE | 1.37|10.96|52.31 |
| HumanNet |1.36|12.09|52.27 |
| HuMAP | 1.52|10.13|51.81 |
| HuRI | 2.23|7.81|51.02 |
| OmniPath |1.43|11.05|51.84 |
| PCNet |1.35|8.69|51.72 |
| ProteomeHD | 7.01|6.83|49.73 |
| SIGNOR | 2.44|9.36|52.08 |
| STRING | 1.37|10.95|52.26 |

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


### Embeddings Download

To download the embeddings generated in this study, use the following link.

- **Embeddings(GOA, Gene2Vec, DNABert, OntoProtein)**: [Download Embeddings](https://drive.google.com/file/d/1xAVhiQtGgyTgqmmI2FOO4VpVYRlbgK30/view?usp=sharing)

To run the downstream tasks, follow the steps at [overall_results.ipynb](overeall_results.ipynb)

### Our Results

![Figure 4](https://github.com/MM-YY-WW/UniEntrezDB/blob/main/Figures/Results.png)










