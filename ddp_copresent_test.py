import argparse
import os
import random
import torch
from tqdm.auto import tqdm
import numpy as np
from tqdm import tqdm
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, f1_score 
from utils import *
from GGIPNN_pytorch import GGIPNN
import torch.distributed as dist
from torch.utils.data import DataLoader, Dataset, DistributedSampler
from torch.nn.parallel import DistributedDataParallel as DDP

class GenePairDataset_embinput(Dataset):
    def __init__(self, dataframe, gene_emb_dict, gene2vec_emb_dict):
        self.dataframe = dataframe
        self.gene_emb_dict = gene_emb_dict
        self.gene2vec_emb_dict = gene2vec_emb_dict
        
    def __len__(self):
        return len(self.dataframe)

    def __getitem__(self, idx):
        row = self.dataframe.iloc[idx]
        if self.gene2vec_emb_dict == {}:
            gene1, gene2, label = self.gene_emb_dict[str(row['Gene1'])], self.gene_emb_dict[str(row['Gene2'])], row['Co_present']
        else:
            gene1 = torch.cat([self.gene_emb_dict[str(row['Gene1'])],self.gene2vec_emb_dict[str(row['Gene1'])]])
            gene2 = torch.cat([self.gene_emb_dict[str(row['Gene2'])],self.gene2vec_emb_dict[str(row['Gene2'])]])
            label = row['Co_present']
        return torch.FloatTensor(gene1), torch.FloatTensor(gene2), label
        

def create_fold_dataloaders(train_val_df, test_df, batch_size, num_folds, rank, world_size, gene_emb_dict, gene2vec_emb_dict):
    # Create dataset for the test set
    test_dataset = GenePairDataset_embinput(test_df, gene_emb_dict=gene_emb_dict, gene2vec_emb_dict=gene2vec_emb_dict)
    test_sampler = DistributedSampler(test_dataset, num_replicas=world_size, rank=rank, shuffle=False)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, sampler=test_sampler)

    skf = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=42)
    fold_dataloaders = []
    for train_idx, val_idx in skf.split(train_val_df, train_val_df['Co_present']):
        train_df = train_val_df.iloc[train_idx]
        val_df = train_val_df.iloc[val_idx]

        train_dataset = GenePairDataset_embinput(train_df, gene_emb_dict=gene_emb_dict, gene2vec_emb_dict=gene2vec_emb_dict)
        val_dataset = GenePairDataset_embinput(val_df, gene_emb_dict=gene_emb_dict, gene2vec_emb_dict=gene2vec_emb_dict)

        train_sampler = DistributedSampler(train_dataset, num_replicas=world_size, rank=rank, shuffle=True)
        val_sampler = DistributedSampler(val_dataset, num_replicas=world_size, rank=rank, shuffle=True)

        train_loader = DataLoader(train_dataset, batch_size=batch_size, sampler=train_sampler)
        val_loader = DataLoader(val_dataset, batch_size=batch_size, sampler=val_sampler)

        fold_dataloaders.append((train_loader, val_loader))
    return fold_dataloaders, test_loader


def train(model, dataloader, optimizer, criterion, device, epoch, scheduler):
    model.train()
    dist.barrier()
    dataloader.sampler.set_epoch(epoch)
    total_loss = 0
    all_preds = []
    all_labels = []

    for batch in dataloader:
        embeddings_gene1 = batch[0].to(device)
        embeddings_gene2 = batch[1].to(device)
        labels = batch[2].to(device)
        optimizer.zero_grad()
        embeddings = torch.cat([embeddings_gene1, embeddings_gene2], dim=-1)
        outputs = model(embeddings)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
        preds = outputs.argmax(dim=-1)
        all_preds.extend(preds.detach().cpu().numpy())
        all_labels.extend(labels.detach().cpu().numpy())

    accuracy = accuracy_score(all_labels, all_preds)
    f1 = f1_score(all_labels, all_preds, average='macro')
    scheduler.step()
    
    return total_loss / len(dataloader), accuracy, f1

def evaluate(model, dataloader, criterion, device, epoch):
    model.eval()
    dist.barrier()
    total_loss = 0
    all_preds = []
    all_labels = []

    with torch.no_grad():
        for batch in dataloader:
            embeddings_gene1 = batch[0].to(device)
            embeddings_gene2 = batch[1].to(device)
            labels = batch[2].to(device)
            embeddings = torch.cat([embeddings_gene1, embeddings_gene2], dim=-1)
            outputs = model(embeddings)
            loss = criterion(outputs, labels)
            total_loss += loss.item()
            preds = outputs.argmax(dim=-1)
            all_preds.extend(preds.cpu().numpy())
            all_labels.extend(labels.cpu().numpy())

    accuracy = accuracy_score(all_labels, all_preds)
    f1 = f1_score(all_labels, all_preds, average='macro')
    return total_loss / len(dataloader), accuracy, f1

def seed_all(seed_value, cuda_deterministic=False):
    """
    set all random seeds
    """
    random.seed(seed_value)
    os.environ['PYTHONHASHSEED'] = str(seed_value)
    np.random.seed(seed_value)
    torch.manual_seed(seed_value)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed_value)
        torch.cuda.manual_seed_all(seed_value)
    # Speed-reproducibility tradeoff https://pytorch.org/docs/stable/notes/randomness.html
    if cuda_deterministic:  # slower, more reproducible
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    else:  # faster, less reproducible
        torch.backends.cudnn.deterministic = False
        torch.backends.cudnn.benchmark = True
        
def save_best_ckpt(epoch, model, optimizer, scheduler, losses, model_name, ckpt_folder):
    """
    save checkpoint
    """
    torch.save(
        {
            'epoch': epoch,
            'model_state_dict': model.module.state_dict(),
            #'to_out_state_dict': 
            'optimizer_state_dict': optimizer.state_dict(),
            'scheduler_state_dict': scheduler.state_dict(),
            'losses': losses,
        },
        f'{ckpt_folder}{model_name}_best.pth'
    )
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epoch", type=int, default=100, help='Number of epochs.')
    parser.add_argument("--seed", type=int, default=2021, help='Random seed.')
    parser.add_argument("--batch_size", type=int, default=233, help='Number of batch size.')
    parser.add_argument("--learning_rate", type=float, default=1e-4, help='Learning rate.')
    parser.add_argument("--device", type=int, default=1, help='which gpu to use')
    parser.add_argument("--dataset", type=str, default='gene2vec', help="choose the dataset used to finetune from [mye, ms, pancreas, zheng68k]")
    parser.add_argument("--num_folds", type=int, default=5, help="k for k fold cross validation")
    parser.add_argument("--current_fold", type=int, default=0, help="perform experiment of which fold")
    parser.add_argument("--local_rank", type=int, default=-1,
                        help="local rank for ddp")
    parser.add_argument("--model_name", type=str, default='', help="the model name of current experiment")
    parser.add_argument("--train_file_path", type=str, default= 'Benchmark_Evaluation_data/MsigDB/train_pairs.csv', 
                        help="the path to save or load positive gene pairs")
    parser.add_argument("--test_file_path", type=str, default= 'Benchmark_Evaluation_data/MsigDB/test_pairs.csv', 
                        help="the path to save or load negative gene pairs")
    parser.add_argument('--gene_emb_path', type=str, default='go_bert/Gene2Vec_baseline/gene2vec_dim_200_iter_9.txt', 
                        help='this file should be a txt file with the first column be the gene entrez id and the rest columns be the gene embedding, dtype float32')
    parser.add_argument("--gene_entrezID_path", type=str, default='', help="corresponding geneEntrezID")
    parser.add_argument("--classifier_ckpt_path", type=str, default='', 
                        help="path to load classifier checkpoint")
    parser.add_argument("--concat_emb", type=str, default="", help="if pass this argument, concat the gobert embedding with gene2vec embedding")
    args = parser.parse_args()
    rank = int(os.environ["RANK"])
    local_rank = args.local_rank
    is_master = local_rank == 0
    
    dist.init_process_group(backend='nccl')
    torch.cuda.set_device(local_rank)
    print(local_rank)
    device = torch.device("cuda", local_rank)
    world_size = torch.distributed.get_world_size()
    seed_all(args.seed + torch.distributed.get_rank())

    gene_emb = torch.load(args.gene_emb_path)
    gene_entrezID = list(pd.read_csv(args.gene_entrezID_path, header=None)[0])
    gene_emb_normalized = gene_emb
    gene_emb_dict = {str(int(id)):ge for id, ge in zip(gene_entrezID, gene_emb_normalized)}
    emb_size = len(gene_emb[0])
    gene2vec_emb_dict = {}
    if args.concat_emb == 'gene2vec':
        gene2vec_emb = torch.load('Embeddings/Gene2Vec/gene2vec.pt')
        gene2vec_entrezID = list(pd.read_csv('Embeddings/Gene2Vec/gene2vec_id.csv', header=None)[0])
        gene2vec_emb_dict = {str(int(id2)):ge2 for id2, ge2 in zip(gene2vec_entrezID, gene2vec_emb)}
        emb_size += len(gene2vec_emb[0])
    elif args.concat_emb == 'dnabert':
        gene2vec_emb = torch.load('Embeddings/DNABert/dnabert_all.pt')
        gene2vec_entrezID = list(pd.read_csv('Embeddings/DNABert/dnabert_allid.csv', header=None)[0])
        gene2vec_emb_dict = {str(int(id2)):ge2 for id2, ge2 in zip(gene2vec_entrezID, gene2vec_emb)}
        emb_size += len(gene2vec_emb[0])


    # Define a simple classifier
    model = GGIPNN(sequence_length=2, num_classes=2, embedding_size=emb_size, hidden_dimension=100)
    if args.classifier_ckpt_path != '':
        path = args.classifier_ckpt_path
        ckpt = torch.load(path)
        model.load_state_dict(ckpt['model_state_dict'])
    model.to(device)    
    model = DDP(model, device_ids=[local_rank], output_device=local_rank, find_unused_parameters=True)

    train_val_df = pd.read_csv(args.train_file_path, sep='\t')
    test_df = pd.read_csv(args.test_file_path, sep='\t')
    dataloaders, test_loader = create_fold_dataloaders(train_val_df=train_val_df,test_df=test_df,  batch_size=args.batch_size, num_folds=args.num_folds, rank=rank, world_size=world_size, gene_emb_dict = gene_emb_dict, gene2vec_emb_dict=gene2vec_emb_dict)
    criterion = torch.nn.CrossEntropyLoss()

    dist.barrier()        
    for epoch in tqdm(range(1, args.epoch + 1)):
        test_loss, test_acc, test_f1 = evaluate(model, test_loader, criterion, device, epoch)
        if is_master:
            print(f'Fold {args.current_fold}, Epoch {epoch}/{args.epoch}')

            #record the results to CSV file
            all_results = pd.read_csv("results_overview.csv", sep='\t')
            if args.model_name[11:-5] not in list(all_results['Model_name']):
                new_row = {col: 0 for col in all_results.columns}
                new_row['Model_name'] = args.model_name[11:-5]
                all_results = pd.concat([all_results, pd.DataFrame([new_row])], ignore_index=True)
            all_results.loc[all_results['Model_name'] == args.model_name[11:-5], f'Copresent_Fold{args.current_fold}_Acc'] = round(test_acc, 4)
            all_results.loc[all_results['Model_name'] == args.model_name[11:-5], f'Copresent_Fold{args.current_fold}_F1'] = round(test_f1, 4)
            all_results.to_csv("results_overview.csv", sep='\t', header=True, index=False)
                
            print(f'Test Loss: {test_loss:.4f}, Test Accuracy: {test_acc:.4f}, Val F1: {test_f1:.4f}')

if __name__ == "__main__":
    main()