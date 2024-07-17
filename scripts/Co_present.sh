#!/bin/bash

choose_emb=$1

# Define the common parameters
epoch=100
seed=0
learning_rate=1e-3
batch_size=233
dataset="msigdb"
num_folds=5
classifier_ckpt_path=''
master_port=12368

if [ "$choose_emb" == "Gene2Vec" ]; then
  gene_entrezID_path="Embeddings/Gene2Vec/gene2vec_id.csv"
  gene_emb_path="Embeddings/Gene2Vec/gene2vec.pt"
  model_name="Co_present_Gene2Vec"
  concat_emb=""

elif [ "$choose_emb" == "DNABert-2" ]; then
  gene_entrezID_path="Embeddings/DNABert/dnabert_allid.csv"
  gene_emb_path="Embeddings/DNABert/dnabert_all.pt"
  model_name="Co_present_DNABert2"
  concat_emb=""

elif [ "$choose_emb" == "OntoProtein" ]; then
  gene_entrezID_path="Embeddings/ontoprotein/ids.csv"
  gene_emb_path="Embeddings/ontoprotein/ontoprotein.pt"
  model_name="Co_present_OntoProtein"
  concat_emb=""

elif [ "$choose_emb" == "GOA_Emb" ]; then
  gene_entrezID_path="Embeddings/GOA/gene_entrezID.csv"
  gene_emb_path="Embeddings/GOA/mean_20001.pt"
  model_name="Co_present_GOA_Emb"
  concat_emb=""

elif [ "$choose_emb" == "GOA_Emb+Gene2Vec" ]; then
  gene_entrezID_path="Embeddings/GOA/gene_entrezID.csv"
  gene_emb_path="Embeddings/GOA/mean_20001.pt"
  model_name="Co_present_GOA_Emb_Gene2Vec"
  concat_emb="gene2vec"

elif [ "$choose_emb" == "GOA_Emb+DNABert-2" ]; then
  gene_entrezID_path="Embeddings/GOA/gene_entrezID.csv"
  gene_emb_path="Embeddings/GOA/mean_20001.pt"
  model_name="Co_present_GOA_Emb_DNABert2"
  concat_emb="dnabert"
else
  echo "Invalid embedding choice: $choose_emb"
  exit 1
fi

# Training and testing loop
for current_fold in {0..4}; do
  # Training command
  CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7,8 python -m torch.distributed.launch --nnodes=1 --node_rank=0 --nproc_per_node=8 --master_port=$master_port ddp_copresent.py \
  --epoch=$epoch \
  --seed=$seed \
  --learning_rate=$learning_rate \
  --batch_size=$batch_size \
  --dataset=$dataset \
  --num_folds=$num_folds \
  --classifier_ckpt_path=$classifier_ckpt_path \
  --gene_entrezID_path=$gene_entrezID_path \
  --gene_emb_path=$gene_emb_path \
  --model_name="$model_name" \
  --current_fold=$current_fold \
  --concat_emb=$concat_emb

  # Testing command
  CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7,8 python -m torch.distributed.launch --nnodes=1 --node_rank=0 --nproc_per_node=8 --master_port=$master_port ddp_copresent_test.py \
  --classifier_ckpt_path="ckpts/${model_name}_${current_fold}_best.pth" \
  --epoch=1 \
  --seed=$seed \
  --learning_rate=$learning_rate \
  --batch_size=$batch_size \
  --dataset=$dataset \
  --num_folds=$num_folds \
  --gene_entrezID_path=$gene_entrezID_path \
  --gene_emb_path=$gene_emb_path \
  --model_name="${model_name}_test" \
  --current_fold=$current_fold \
  --concat_emb=$concat_emb 

done

# Additional test commands outside the loop
for current_fold in {0..4}; do
# for current_fold in {0..0}; do
  CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7,9 python -m torch.distributed.launch --nnodes=1 --node_rank=0 --nproc_per_node=8 --master_port=$master_port ddp_copresent_test.py \
  --classifier_ckpt_path="ckpts/${model_name}_${current_fold}_best.pth" \
  --epoch=1 \
  --seed=$seed \
  --learning_rate=$learning_rate \
  --batch_size=$batch_size \
  --dataset=$dataset \
  --num_folds=$num_folds \
  --gene_entrezID_path=$gene_entrezID_path \
  --gene_emb_path=$gene_emb_path \
  --model_name="${model_name}_test" \
  --current_fold=$current_fold \
  --concat_emb=$concat_emb 

done