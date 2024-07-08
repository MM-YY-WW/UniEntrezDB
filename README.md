# UnientrezDB: Large-scale Gene Ontology Annotation Dataset and Evaluation Benchmarks with Unified Entrez Gene Identifiers

This repository contains the official implementation for the paper "UnientrezDB: Large-scale Gene Ontology Annotation Dataset and Evaluation Benchmarks with Unified Entrez Gene Identifiers." Our work focuses on providing a comprehensive dataset and benchmarks for evaluating gene ontology annotations using a unified system of Entrez Gene Identifiers.

## Figures from the Paper

Below are several figures from the paper that illustrate key concepts and results:

![Figure 1](https://github.com/MM-YY-WW/UniEntrezDB/blob/main/Figures/goa.png)
*Figure 1: Description of what the figure represents.*

![Figure 2](https://github.com/MM-YY-WW/UniEntrezDB/blob/main/Figures/evaluation_benchmark.png)
*Figure 2: Description of what the figure represents.*

![Figure 3](path/to/figure3.png)
*Figure 3: Description of what the figure represents.*

## Dataset Download

You can download the datasets used in our study through the following links.

- **Pretrained Dataset**: [Download Pretrained Dataset](https://drive.google.com/file/d/1DsXufybeSgEXrx8szkF0kuhASmAVOaU-/view?usp=sharing)
- **Downstream Task Datasets**: [Download Datasets](https://drive.google.com/file/d/1fSRXO26jr1XcFn7GKqRoN_CZUbuEY8Cj/view?usp=sharing)
- **ID Mapping Dicts**: [Download ID Mapping dictionaries](https://drive.google.com/file/d/1La80B3hUibbe94FghkTIx80DRzPfwYix/view?usp=sharing)

## Embeddings Download

To download the embeddings generated in this study, use the following links. Each link corresponds to embeddings tailored for different aspects of the study:

- **Embedding 1**: [Download Embedding 1](https://example.com/embedding1)
- **Embedding 2**: [Download Embedding 2](https://example.com/embedding2)
- **Embedding 3**: [Download Embedding 3](https://example.com/embedding3)

## Running Downstream Tasks

To run the downstream tasks, follow these steps. Below are also sample results and code snippets that you can use to reproduce the findings of our study.

### Requirements

Before running the scripts, make sure to install all required packages:

```bash
pip install -r requirements.txt

```

### Use Gene Embedding as Input 

The format of gene embedding is .pt file (default dimension 1024) with the corresponding Gene Entrez ID in .csv format


### Script

To test embedding performance of the downstream tasks:

```
bash all_in_one.sh \
-Co_Present True \
-GGI True \
-PPI True \
-Cell_Type True \
-Gene_Embedding_Path '' \
-Gene_ID_Path '' \
-Result_Folder_Path '' \
```

## 