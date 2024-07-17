# -*- coding: utf-8 -*-
import os
import argparse
import random
import random
from functools import reduce
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, ShuffleSplit, StratifiedShuffleSplit, StratifiedKFold
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix, precision_recall_fscore_support, classification_report
import torch
from torch import nn
from torch.optim import Adam, SGD, AdamW
from torch.utils.data import DataLoader, Dataset
from torch.utils.data.distributed import DistributedSampler
from torch.nn.parallel import DistributedDataParallel as DDP
import torch.distributed as dist
from tqdm import tqdm
# from performer_pytorch import PerformerLM
import scanpy as sc
import anndata as ad
from utils import *
import pickle as pkl
parser = argparse.ArgumentParser()
parser.add_argument("--local_rank", type=int, default=-1, help='Local process rank.')
parser.add_argument("--bin_num", type=int, default=5, help='Number of bins.')
parser.add_argument("--gene_num", type=int, default=16906, help='Number of genes.')
parser.add_argument("--epoch", type=int, default=100, help='Number of epochs.')
parser.add_argument("--seed", type=int, default=2021, help='Random seed.')
parser.add_argument("--batch_size", type=int, default=111, help='Number of batch size.')
parser.add_argument("--learning_rate", type=float, default=1e-4, help='Learning rate.')
parser.add_argument("--grad_acc", type=int, default=60, help='Number of gradient accumulation.')
parser.add_argument("--valid_every", type=int, default=1, help='Number of training epochs between twice validation.')
parser.add_argument("--pos_embed", type=bool, default=True, help='Using Gene2vec encoding or not.')
parser.add_argument("--data_path", type=str, default='./data/Zheng68K.h5ad', help='Path of data for finetune.')
parser.add_argument("--model_path", type=str, default='', help='Path of pretrained model.')
parser.add_argument("--ckpt_dir", type=str, default='./ckpts/', help='Directory of checkpoint to save.')
parser.add_argument("--model_name", type=str, default='finetune', help='Finetuned model name.')
parser.add_argument("--current_fold", type=int, default=0)
parser.add_argument("--text_emb_path", type=str, default="")
parser.add_argument("--text_method", type=str, default="text", help="origin, mask_notext, text")
parser.add_argument("--all_gene_num", type=int, default=0)
parser.add_argument("--second_text_emb_path", type=str, default="")
parser.add_argument("--use_first_emb", type=int, default=0, help="if set to 1, then use the first_emb")
args = parser.parse_args()
rank = int(os.environ["RANK"])
local_rank = args.local_rank
is_master = local_rank == 0
if is_master:
    log_directory = f"./log"
    if not os.path.exists(log_directory):
        os.makedirs(log_directory)
    log_file_path = os.path.join(log_directory, f"{args.model_name}_{args.current_fold}.log")
    logging.basicConfig(filename=log_file_path, level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', filemode='a')
    logger= logging.getLogger()
    
    
SEED = args.seed
EPOCHS = args.epoch
BATCH_SIZE = args.batch_size
GRADIENT_ACCUMULATION = args.grad_acc
LEARNING_RATE = args.learning_rate
SEQ_LEN = args.gene_num + 1
VALIDATE_EVERY = args.valid_every

PATIENCE = 100
UNASSIGN_THRES = 0.0

CLASS = args.bin_num + 2
POS_EMBED_USING = args.pos_embed

model_name = args.model_name
ckpt_dir = args.ckpt_dir

dist.init_process_group(backend='nccl')
torch.cuda.set_device(local_rank)
device = torch.device("cuda", local_rank)
world_size = torch.distributed.get_world_size()

seed_all(SEED + torch.distributed.get_rank())


class SCDataset(Dataset):
    def __init__(self, data, label):
        super().__init__()
        self.data = data
        self.label = label

    def __getitem__(self, index):
        rand_start = random.randint(0, self.data.shape[0]-1)
        full_seq = self.data[rand_start].toarray()[0]
        full_seq[full_seq > (CLASS - 2)] = CLASS - 2
        full_seq = torch.from_numpy(full_seq).long()
        full_seq = torch.cat((full_seq, torch.tensor([0]))).to(device)
        seq_label = self.label[rand_start]
        return full_seq, seq_label

    def __len__(self):
        return self.data.shape[0]
    
    
def create_fold_dataloaders(data, batch_size, num_folds, rank, world_size):
    sss = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=SEED)
    for index_train, index_val in sss.split(data, label):
        data_train, label_train = data[index_train], label[index_train]
        data_test, label_test = data[index_val], label[index_val]
        
        # Create dataset for the test set
        test_dataset = SCDataset(data_test, label_test)
        test_sampler = DistributedSampler(test_dataset, num_replicas=world_size, rank=rank, shuffle=False)
        test_loader = DataLoader(test_dataset, batch_size=batch_size, sampler=test_sampler)

        skf = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=SEED)
        fold_dataloaders = []
        for train_idx, val_idx in skf.split(data_train, label_train):
            fdata_train, flabel_train = data_train[train_idx], label_train[train_idx]
            fdata_val, flabel_val = data_train[val_idx], label_train[val_idx]

            train_dataset = SCDataset(fdata_train, flabel_train)
            val_dataset = SCDataset(fdata_val, flabel_val)

            train_sampler = DistributedSampler(train_dataset, num_replicas=world_size, rank=rank, shuffle=True)
            val_sampler = DistributedSampler(val_dataset, num_replicas=world_size, rank=rank, shuffle=True)

            train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, sampler=train_sampler)
            val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE, sampler=val_sampler)

            fold_dataloaders.append((train_loader, val_loader))
    return fold_dataloaders, test_loader, test_sampler, train_sampler, val_sampler   

class Identity_concat(torch.nn.Module):
    def __init__(self, dropout = 0., h_dim = 100, out_dim = 10, text_dim = 1024, args=None):
        super(Identity_concat, self).__init__()
        self.conv1 = nn.Conv2d(1, 1, (1, text_dim+1))
        self.act = nn.ReLU()
        self.fc1 = nn.Linear(in_features=args.gene_num, out_features=512, bias=True)
        self.act1 = nn.ReLU()
        self.dropout1 = nn.Dropout(dropout)
        self.fc2 = nn.Linear(in_features=512, out_features=h_dim, bias=True)
        self.act2 = nn.ReLU()
        self.dropout2 = nn.Dropout(dropout)
        self.fc3 = nn.Linear(in_features=h_dim, out_features=out_dim, bias=True)   
        
    def forward(self, x, text_embedding):
        #'concat': cat([N, M, K], [N, M, T]) = [N, M, K+T]
        #x [N, M, K]
        non_zero_mask = text_embedding.sum(dim=1) != 0
        x=x[:,non_zero_mask]
        text_embedding= text_embedding[non_zero_mask]
        text_embedding = text_embedding.repeat(x.shape[0], 1,1) #[N, M, T]
        x=x.unsqueeze(-1)
        x = torch.cat((x, text_embedding), dim=2) #[N, M, K+T]
        x = x[:,None,:,:]
        x = self.conv1(x)
        x = self.act(x)
        x = x.view(x.shape[0],-1)
        x = self.fc1(x)
        x = self.act1(x)
        x = self.dropout1(x)
        x = self.fc2(x)
        x = self.act2(x)
        x = self.dropout2(x)
        x = self.fc3(x)
        return x


data = sc.read_h5ad(args.data_path)
label_dict, label = np.unique(np.array(data.obs['celltype']), return_inverse=True)  # Convert strings categorical to integrate categorical, and label_dict[label] can be restored
#store the label dict and label for prediction
with open('label_dict', 'wb') as fp:
    pkl.dump(label_dict, fp)
with open('label', 'wb') as fp:
    pkl.dump(label, fp)
class_num = np.unique(label, return_counts=True)[1].tolist()
class_weight = torch.tensor([(1 - (x / sum(class_num))) ** 2 for x in class_num])
label = torch.from_numpy(label)
data = data.X

acc = []
f1 = []
f1w = []
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
pred_list = pd.Series(['un'] * data.shape[0])

dataloaders, test_loader, test_sampler, train_sampler, val_sampler = create_fold_dataloaders(data=data, batch_size=args.batch_size, num_folds=5, rank=rank, world_size=world_size)

train_loader, val_loader = dataloaders[args.current_fold]

text_embedding = None

if args.text_method =="text":
    text_embedding = np.loadtxt(args.text_emb_path) 
    #text_embedding = np.vstack((text_embedding,np.zeros((1, text_embedding.shape[1])))).astype(np.float32)     
    text_embedding = torch.tensor(text_embedding, requires_grad=False)
    text_embedding = torch.nn.functional.normalize(text_embedding, p=2, dim=1)
    if args.second_text_emb_path:
        second_text_embedding = np.loadtxt(args.second_text_emb_path)
        second_text_embedding = torch.tensor(second_text_embedding, requires_grad=False)
        second_text_embedding = torch.nn.functional.normalize(second_text_embedding, p=2, dim=1)
        if args.use_first_emb == 1:
            concatenated_embeddings = np.zeros((text_embedding.shape[0], text_embedding.shape[1]), dtype=np.float32)
            for i in range(text_embedding.shape[0]):
                if text_embedding[i].sum(dim=0) != 0 and second_text_embedding[i].sum(dim=0) != 0:
                    concatenated_embeddings[i] = text_embedding[i]
        else:    
            concatenated_embeddings = np.zeros((text_embedding.shape[0], text_embedding.shape[1] + second_text_embedding.shape[1]), dtype=np.float32)
            for i in range(text_embedding.shape[0]):
                if text_embedding[i].sum(dim=0) != 0 and second_text_embedding[i].sum(dim=0) != 0:
                    concatenated_embeddings[i] = torch.cat([text_embedding[i], second_text_embedding[i]])
        args.gene_num = sum([1 for i in concatenated_embeddings if i.sum() != 0])
        text_embedding = torch.tensor(concatenated_embeddings, requires_grad=False)
        
    zero_row = torch.zeros(1, text_embedding.shape[1], dtype=torch.float32)
    text_embedding = torch.cat((text_embedding, zero_row), dim=0)
    model = Identity_concat(dropout=0., h_dim=128, out_dim=label_dict.shape[0], text_dim = text_embedding.shape[1],args=args)

    
if args.model_path != '':
    path = args.model_path
    ckpt = torch.load(path)
    model.load_state_dict(ckpt['model_state_dict'])
model = model.to(device)
model = DDP(model, device_ids=[local_rank], output_device=local_rank)


optimizer = Adam(model.parameters(), lr=LEARNING_RATE)
scheduler = CosineAnnealingWarmupRestarts(
    optimizer,
    first_cycle_steps=15,
    cycle_mult=2,
    max_lr=LEARNING_RATE,
    min_lr=1e-6,
    warmup_steps=5,
    gamma=0.9
)
loss_fn = nn.CrossEntropyLoss(weight=None).to(local_rank)

dist.barrier()
trigger_times = 0
max_acc = 0.0
for i in tqdm(range(1, EPOCHS+1)):
    train_loader.sampler.set_epoch(i)
    model.train()
    dist.barrier()
    running_loss = 0.0
    cum_acc = 0.0
    for index, (data, labels) in enumerate(train_loader):
        index += 1
        data, labels = data.to(device), labels.to(device)
        if index % GRADIENT_ACCUMULATION != 0:
            with model.no_sync():
                logits = model(data, text_embedding)
                loss = loss_fn(logits, labels)
                loss.backward()
        if index % GRADIENT_ACCUMULATION == 0:
            logits = model(data, text_embedding)
            loss = loss_fn(logits, labels)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), int(1e6))
            optimizer.step()
            optimizer.zero_grad()
        running_loss += loss.item()
        softmax = nn.Softmax(dim=-1)
        final = softmax(logits)
        final = final.argmax(dim=-1)
        pred_num = labels.size(0)
        correct_num = torch.eq(final, labels).sum(dim=-1)
        cum_acc += torch.true_divide(correct_num, pred_num).mean().item()
    epoch_loss = running_loss / index
    epoch_acc = 100 * cum_acc / index
    epoch_loss = get_reduced(epoch_loss, local_rank, 0, world_size)
    epoch_acc = get_reduced(epoch_acc, local_rank, 0, world_size)
    if is_master:
        print(f'    == Fold {args.current_fold} Epoch: {i} | Training Loss: {epoch_loss:.6f} | Accuracy: {epoch_acc:6.4f}%  ==')
        logging.info(f'    ==  Epoch: {i} | Training Loss: {epoch_loss:.6f} | Accuracy: {epoch_acc:6.4f}%  ==')
        print(optimizer.param_groups[0]['lr'])
    dist.barrier()
    scheduler.step()

    if i % VALIDATE_EVERY == 0:
        model.eval()
        dist.barrier()
        running_loss = 0.0
        predictions = []
        truths = []
        with torch.no_grad():
            for index, (data_v, labels_v) in enumerate(val_loader):
                index += 1
                data_v, labels_v = data_v.to(device), labels_v.to(device)
                logits = model(data_v, text_embedding)
                loss = loss_fn(logits, labels_v)
                running_loss += loss.item()
                softmax = nn.Softmax(dim=-1)
                final_prob = softmax(logits)
                final = final_prob.argmax(dim=-1)
                final[np.amax(np.array(final_prob.cpu()), axis=-1) < UNASSIGN_THRES] = -1
                predictions.append(final)
                truths.append(labels_v)
            del data_v, labels_v, logits, final_prob, final
            # gather
            predictions = distributed_concat(torch.cat(predictions, dim=0), len(val_sampler.dataset), world_size)
            truths = distributed_concat(torch.cat(truths, dim=0), len(val_sampler.dataset), world_size)
            no_drop = predictions != -1
            predictions = np.array((predictions[no_drop]).cpu())
            truths = np.array((truths[no_drop]).cpu())
            cur_acc = accuracy_score(truths, predictions)
            f1 = f1_score(truths, predictions, average='macro')
            val_loss = running_loss / index
            val_loss = get_reduced(val_loss, local_rank, 0, world_size)
            if is_master:
                print(f'    == Fold {args.current_fold} Epoch: {i} | validation acc {cur_acc:.6f}  Validation Loss: {val_loss:.6f} | F1 Score: {f1:.6f}  ==')
                print(confusion_matrix(truths, predictions))
                #print(classification_report(truths, predictions, target_names=label_dict.tolist(), digits=4))
                logging.info(f'    ==  Epoch: {i} | validation acc {cur_acc:.6f} Validation Loss: {val_loss:.6f} | F1 Score: {f1:.6f}  ==')
                logging.info(confusion_matrix(truths, predictions))
                #logging.info(classification_report(truths, predictions, target_names=label_dict.tolist(), digits=4))
            if cur_acc > max_acc:
                max_acc = cur_acc
                trigger_times = 0
                save_best_ckpt(i, model, optimizer, scheduler, val_loss, f'{model_name}_{str(args.current_fold)}', ckpt_dir)
            else:
                trigger_times += 1
                if trigger_times > PATIENCE:
                    break
    del predictions, truths