#!/bin/bash

choose_emb=$1

# Define the common parameters
epoch=100
seed=2021
learning_rate=1e-3
batch_size=111
dataset="zheng68k"
num_folds=5
master_port=12368
data_path="Benchmark_Evaluation_data/Zheng68k/Zheng68K.h5ad"
grad_acc=1

if [ "$choose_emb" == "Gene2Vec" ]; then
  model_name="SingleCell_Gene2Vec"
  text_emb_path="Embeddings/Gene2Vec/Gene2Vec_15779_16906.txt"
  second_text_emb_path="Embeddings/GOA/yesIEA_yesRoot_8467_16906.txt"
  use_first_emb=1

elif [ "$choose_emb" == "DNABert-2" ]; then
  model_name="SingleCell_DNABert2"
  text_emb_path="Embeddings/DNABert/10185_16906.txt"
  second_text_emb_path="Embeddings/GOA/yesIEA_yesRoot_8467_16906.txt"
  use_first_emb=1

elif [ "$choose_emb" == "OntoProtein" ]; then
  model_name="SingleCell_OntoProtein"
  text_emb_path="Embeddings/ontoprotein/16873_16906.txt"
  second_text_emb_path="Embeddings/GOA/yesIEA_yesRoot_8467_16906.txt"
  use_first_emb=1

elif [ "$choose_emb" == "GOA_Emb" ]; then
  model_name="SingleCell_GOA_Emb"
  text_emb_path="Embeddings/GOA/yesIEA_yesRoot_8467_16906.txt"
  second_text_emb_path="Embeddings/Gene2Vec/Gene2Vec_15779_16906.txt"
  use_first_emb=1

elif [ "$choose_emb" == "GOA_Emb+Gene2Vec" ]; then
  model_name="SingleCell_GOA_Emb_Gene2Vec"
  text_emb_path="Embeddings/Gene2Vec/Gene2Vec_15779_16906.txt"
  second_text_emb_path="Embeddings/GOA/yesIEA_yesRoot_8467_16906.txt"
  use_first_emb=0


elif [ "$choose_emb" == "GOA_Emb+DNABert-2" ]; then
  model_name="SingleCell_GOA_Emb_DNABert2"
  text_emb_path="Embeddings/DNABert/10185_16906.txt"
  second_text_emb_path="Embeddings/GOA/yesIEA_yesRoot_8467_16906.txt"
  use_first_emb=0



else
  echo "Invalid embedding choice: $choose_emb"
  exit 1
fi

# Training and testing loop
for current_fold in {0..0}; do
# for current_fold in {0..4}; do
  # Training command
  CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7,8 python -m torch.distributed.launch --nnodes=1 --node_rank=0 --nproc_per_node=8 --master_port=$master_port ddp_singlecell.py \
 --epoch $epoch \
 --batch_size $batch_size \
 --seed $seed \
 --grad_acc $grad_acc \
 --learning_rate $learning_rate \
 --data_path $data_path \
 --model_name $model_name \
 --text_emb_path $text_emb_path \
 --second_text_emb_path $second_text_emb_path \
 --use_first_emb $use_first_emb \
 --current_fold $current_fold

  # Testing command
  CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7,8 python -m torch.distributed.launch --nnodes=1 --node_rank=0 --nproc_per_node=8 --master_port=$master_port ddp_singlecell_test.py \
  --model_path="ckpts/${model_name}_${current_fold}_best.pth" \
  --epoch 1 \
  --batch_size $batch_size \
  --seed $seed \
  --grad_acc $grad_acc \
  --learning_rate $learning_rate \
  --data_path $data_path \
  --model_name="${model_name}_test" \
  --text_emb_path $text_emb_path \
  --second_text_emb_path $second_text_emb_path \
  --use_first_emb $use_first_emb \
  --current_fold $current_fold

done
