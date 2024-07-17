#from Bio import Entrez
import numpy as np
import json
#import requests
import scanpy as sc
import torch
import numpy as np
import pandas as pd
import math
from sklearn.metrics import roc_auc_score, accuracy_score, f1_score, precision_score, recall_score
import random 
import os
import matplotlib.pyplot as plt
from datetime import datetime
from tqdm import tqdm
import logging
#import umap
from PIL import Image, ImageDraw
from torch.optim.lr_scheduler import _LRScheduler


def get_entrez_id(gene_name):
    # Always provide your email address to NCBI to let them know who is using their services.
    Entrez.email = "your_email@example.com"
    if 'ENSG' in gene_name:
        try:
            # Search in NCBI Gene database for Ensembl ID
            search_result = Entrez.esearch(db="gene", term=f"{gene_name}[Ensembl]")
            record = Entrez.read(search_result)
            search_result.close()

            # Fetch the first result's Entrez ID
            entrez_ids = record['IdList']
            if entrez_ids:
                return entrez_ids[0]
            else:
                return -1
        except Exception as e:
            return f"An error occurred: {e}"
    else:
        try:
            # Search in NCBI Gene database
            search_result = Entrez.esearch(db="gene", term=f"{gene_name}[Gene Name]")
            record = Entrez.read(search_result)
            search_result.close()

            # Fetch the first result's Entrez ID
            entrez_ids = record['IdList']
            if entrez_ids:
                return entrez_ids[0]
            else:
                return -1
        except Exception as e:
            return f"An error occurred: {e}"

# # Example usage
# gene_name = "MT-ATP6"
# entrez_id = get_entrez_id(gene_name)
# print(f"The Entrez ID for {gene_name} is: {entrez_id}")
class CosineAnnealingWarmupRestarts(_LRScheduler):
    """
        optimizer (Optimizer): Wrapped optimizer.
        first_cycle_steps (int): First cycle step size.
        cycle_mult(float): Cycle steps magnification. Default: -1.
        max_lr(float): First cycle's max learning rate. Default: 0.1.
        min_lr(float): Min learning rate. Default: 0.001.
        warmup_steps(int): Linear warmup step size. Default: 0.
        gamma(float): Decrease rate of max learning rate by cycle. Default: 1.
        last_epoch (int): The index of last epoch. Default: -1.
    """
    
    def __init__(self,
                 optimizer : torch.optim.Optimizer,
                 first_cycle_steps : int,
                 cycle_mult : float = 1.,
                 max_lr : float = 0.1,
                 min_lr : float = 0.001,
                 warmup_steps : int = 0,
                 gamma : float = 1.,
                 last_epoch : int = -1
        ):
        assert warmup_steps < first_cycle_steps
        
        self.first_cycle_steps = first_cycle_steps # first cycle step size
        self.cycle_mult = cycle_mult # cycle steps magnification
        self.base_max_lr = max_lr # first max learning rate
        self.max_lr = max_lr # max learning rate in the current cycle
        self.min_lr = min_lr # min learning rate
        self.warmup_steps = warmup_steps # warmup step size
        self.gamma = gamma # decrease rate of max learning rate by cycle
        
        self.cur_cycle_steps = first_cycle_steps # first cycle step size
        self.cycle = 0 # cycle count
        self.step_in_cycle = last_epoch # step size of the current cycle
        
        super(CosineAnnealingWarmupRestarts, self).__init__(optimizer, last_epoch)
        
        # set learning rate min_lr
        self.init_lr()

    def init_lr(self):
        self.base_lrs = []
        for param_group in self.optimizer.param_groups:
            param_group['lr'] = self.min_lr
            self.base_lrs.append(self.min_lr)
    
    def get_lr(self):
        if self.step_in_cycle == -1:
            return self.base_lrs
        elif self.step_in_cycle < self.warmup_steps:
            return [(self.max_lr - base_lr)*self.step_in_cycle / self.warmup_steps + base_lr for base_lr in self.base_lrs]
        else:
            return [base_lr + (self.max_lr - base_lr) \
                    * (1 + math.cos(math.pi * (self.step_in_cycle-self.warmup_steps) \
                                    / (self.cur_cycle_steps - self.warmup_steps))) / 2
                    for base_lr in self.base_lrs]

    def step(self, epoch=None):
        if epoch is None:
            epoch = self.last_epoch + 1
            self.step_in_cycle = self.step_in_cycle + 1
            if self.step_in_cycle >= self.cur_cycle_steps:
                self.cycle += 1
                self.step_in_cycle = self.step_in_cycle - self.cur_cycle_steps
                self.cur_cycle_steps = int((self.cur_cycle_steps - self.warmup_steps) * self.cycle_mult) + self.warmup_steps
        else:
            if epoch >= self.first_cycle_steps:
                if self.cycle_mult == 1.:
                    self.step_in_cycle = epoch % self.first_cycle_steps
                    self.cycle = epoch // self.first_cycle_steps
                else:
                    n = int(math.log((epoch / self.first_cycle_steps * (self.cycle_mult - 1) + 1), self.cycle_mult))
                    self.cycle = n
                    self.step_in_cycle = epoch - int(self.first_cycle_steps * (self.cycle_mult ** n - 1) / (self.cycle_mult - 1))
                    self.cur_cycle_steps = self.first_cycle_steps * self.cycle_mult ** (n)
            else:
                self.cur_cycle_steps = self.first_cycle_steps
                self.step_in_cycle = epoch
                
        self.max_lr = self.base_max_lr * (self.gamma**self.cycle)
        self.last_epoch = math.floor(epoch)
        for param_group, lr in zip(self.optimizer.param_groups, self.get_lr()):
            param_group['lr'] = lr

def cosine_similarity(vector1, vector2):
    """Calculate the cosine similarity between two vectors."""
    dot_product = np.dot(vector1, vector2)
    norm_vector1 = np.linalg.norm(vector1)
    norm_vector2 = np.linalg.norm(vector2)
    similarity = dot_product / (norm_vector1 * norm_vector2)
    return similarity

def load_json(filepath):
    with open(filepath, 'r') as file:
        result = json.load(file)

    return result

def save_json(file_to_save, filepath, mode):
    with open(filepath, mode) as file:
        json.dump(file_to_save, file, indent=4)


def calculate_metric(y, y_pred):
    accuracy = accuracy_score(y, y_pred)
    f1 = f1_score(y, y_pred, average='weighted', zero_division=0)
    precision = precision_score(y, y_pred, average='weighted', zero_division=0)
    recall = recall_score(y, y_pred, average='weighted', zero_division=0)
    return accuracy, f1, precision, recall

def set_seed(seed):
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.cuda.manual_seed_all(seed)
    if torch.cuda.is_available():
        torch.backends.cudnn.deterministric = True
        torch.backends.cudnn.benchmark = False

def umap_visualization(embedding, labels, path=''):
    reducer = umap.UMAP()
    umap_embeddings = reducer.fit_transform(embedding)
    plt.figure(figsize=(10, 6))
    plt.scatter(umap_embeddings[:, 0], umap_embeddings[:, 1], c=labels, cmap='viridis', s=5)
    plt.colorbar()
    plt.title('UMAP Visualization')
    plt.xlabel('UMAP Dimension 1')
    plt.ylabel('UMAP Dimension 2')
    #path='./gif/temp/' + datetime.now().strftime('%Y_%m_%d_%H_%M_%S') + '.png'
    plt.savefig(path)
    return path

def result_gif(train_frame_paths, test_frame_paths, gif_path):
    if len(train_frame_paths) != len(test_frame_paths):
        raise ValueError("train_frames and test_frames must be of the same length")
    combined_frames = []
    for train_path, test_path in zip(train_frame_paths, test_frame_paths):
        with Image.open(train_path) as train_frame, Image.open(test_path) as test_frame:
            width, height = train_frame.size
            combined_image = Image.new('RGB', (2 * width, height), (255, 255, 255))  
            combined_image.paste(train_frame, (0, 0))
            combined_image.paste(test_frame, (width, 0))
            combined_frames.append(combined_image.copy())
    combined_frames[0].save(gif_path, save_all=True, append_images=combined_frames[1:], duration=500, loop=0)

def plot_results(train_acc_list, test_acc_list, train_precision_list, test_precision_list, train_f1_list, test_f1_list, train_recall_list, test_recall_list, args, fold):
    # Create a directory to store the plots
    if not os.path.exists("./plot_results"):
        os.makedirs("./plot_results")

    # Generate a unique timestamp for the plot filename
    timestamp = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')

    # Set up the figure and subplots
    fig, axs = plt.subplots(2, 4, figsize=(16, 8))
    fig.suptitle(f"lr = {args.learning_rate}, epoch = {args.epoch}, dataset = {args.dataset}, seed = {args.seed}, batch_size = {args.batch_size}, model_name = {args.model_name}, hidden_k_size = {args.hidden_k_size}, missing_summary = {args.missing_summary}, scbert_out_model = {args.scbert_out_model}")
    logging.info(f"lr = {args.learning_rate}, epoch = {args.epoch}, dataset = {args.dataset}, seed = {args.seed}, batch_size = {args.batch_size}, model_name = {args.model_name}, hidden_k_size = {args.hidden_k_size}, missing_summary = {args.missing_summary}, scbert_out_model = {args.scbert_out_model}")
    # Plot training accuracy
    axs[0, 0].plot(train_acc_list, label='Train Accuracy')
    axs[0, 0].set_title('Train Accuracy')

    # Plot testing accuracy
    axs[1, 0].plot(test_acc_list, label='Test Accuracy')
    axs[1, 0].set_title('Test Accuracy')

    # Plot training precision
    axs[0, 1].plot(train_precision_list, label='Train Precision')
    axs[0, 1].set_title('Train Precision')

    # Plot testing precision
    axs[1, 1].plot(test_precision_list, label='Test Precision')
    axs[1, 1].set_title('Test Precision')

    # Plot training F1 score
    axs[0, 2].plot(train_f1_list, label='Train F1 Score')
    axs[0, 2].set_title('Train F1 Score')

    # Plot testing F1 score
    axs[1, 2].plot(test_f1_list, label='Test F1 Score')
    axs[1, 2].set_title('Test F1 Score')

    # Plot training recall
    axs[0, 3].plot(train_recall_list, label='Train Recall')
    axs[0, 3].set_title('Train Recall')

    # Plot testing recall
    axs[1, 3].plot(test_recall_list, label='Test Recall')
    axs[1, 3].set_title('Test Recall')

    # Add legends
    axs[0, 0].legend()
    axs[1, 0].legend()
    axs[0, 1].legend()
    axs[1, 1].legend()
    axs[0, 2].legend()
    axs[1, 2].legend()
    axs[0, 3].legend()
    axs[1, 3].legend()
    if args.model_name == 'scBERT':
        plot_filename = f"./plot_results/{timestamp}_{args.dataset}_{args.model_name}_{args.scbert_out_model}_Fold{fold}.png"
    else:
        plot_filename = f"./plot_results/{timestamp}_{args.dataset}_{args.model_name}_Fold{fold}.png"
    logging.info(plot_filename)
    plt.savefig(plot_filename)
    plt.show()

def get_go_by_idlist(gene_id_list, args):
    """_summary_
        Obtain GO annotation for each gene id 
        return three list:
        GeneID, GO annotation and qualifier is positive or negative
    Args:
        args (_type_): _description_
    """
    evidence_levels = [l.strip() for l in args.evidence_level.split(',')]
    object_types = [t.strip() for t in args.object_type.split(',')]
    qualifier = [q.strip() for q in args.qualifier.split(',')]
    evidence_level_symbol = {'-1': ['TAS', 'ND', 'IC', 'NAS'],
                              '0': ['IEA'],
                              '1': ['ISO', 'IGC', 'ISM', 'ISS', 'RCA', 'ISA'],
                              '2': ['IGI', 'HDA', 'EXP', 'IMP', 'IKR', 'IPI', 'HGI', 'HEP', 'IEP', 'IDA', 'HTP', 'IBA', 'HMP']
                            }
    object_type_dict = {
            'protein': ["protein","protein_complex"],
        
                'DNA': ["sense_intronic_ncRNA_gene","miRNA_gene","sense_overlap_ncRNA_gene","RNase_P_RNA_gene","telomerase_RNA_gene",
                        "lincRNA_gene","pseudogene","tRNA_gene","SRP_RNA_gene","protein_coding_gene","RNase_MRP_RNA_gene","gene_segment",
                        "gene","ncRNA_gene","rRNA_gene","snRNA_gene","scRNA_gene","transposable_element_gene","snoRNA_gene","lncRNA_gene"],
            
                'RNA': ["piRNA","snoRNA","tRNA","snRNA","hammerhead_ribozyme","RNA","scRNA","RNase_P_RNA","lnc_RNA","telomerase_RNA",
                        "ncRNA","RNase_MRP_RNA","guide_RNA","ribozyme","rRNA","miRNA","antisense_lncRNA","mRNA","antisense_RNA",
                        "tmRNA","SRP_RNA"],
            
             'others': ["biological region","gene_product"]
                        }         
    qualifier_level_dict = {
        'pos': ["enables","acts_upstream_of_or_within_positive_effect","acts_upstream_of_or_within","located_in","colocalizes_with",
                "part_of","acts_upstream_of_positive_effect","involved_in","acts_upstream_of","is_active_in",
                "NOT|acts_upstream_of_or_within_negative_effect","contributes_to",],
        
        'neg': ["NOT|located_in","NOT|is_active_in","NOT|involved_in","acts_upstream_of_negative_effect",
                "NOT|acts_upstream_of_or_within_positive_effect","NOT|contributes_to","NOT|part_of","NOT|enables"
                "NOT|acts_upstream_of","NOT|colocalizes_with","NOT|acts_upstream_of_or_within", "acts_upstream_of_or_within_negative_effect",]
    }
    target_evidence_codes = []
    for level in evidence_levels:
        target_evidence_codes += evidence_level_symbol[level]
    
    target_object_types = []
    for type in object_types:
        target_object_types += object_type_dict[type]
                
    target_qualifier = []
    for type in qualifier:
        target_qualifier += qualifier_level_dict[type]      
             
    anno_dict = load_json(args.DBID_EntrezID_fulldict_path)
    #gene_entrezID : [goid, goid, goid...]
    result_gene_EntrezID = []
    result_gene_GOAnno = []
    result_gene_Qualifier = []
    anno_count=0
    gene_count=0
    print(f"\n===Start to obtain chosen annotation: \n evidence_level{evidence_levels} \n object_types: {object_types} \n qualifier: {qualifier} \n")
    for geneid in tqdm(gene_id_list):
        gene_anno_list = []
        gene_anno_qualifier_list = []
        for type, evidence_dict in anno_dict[str(geneid)].items():
            if type in target_object_types:
                for evidence, qualifier_dict in evidence_dict.items():
                    if evidence in target_evidence_codes:
                        for qualifier, goid_list in qualifier_dict.items():
                            if qualifier in target_qualifier:
                                gene_anno_list += goid_list
                                anno_count += len(goid_list)
                                gene_anno_qualifier_list += [qualifier]*len(goid_list)
                                assert len(gene_anno_list) == len(gene_anno_qualifier_list)
        if gene_anno_list != []:
            gene_count +=1
            result_gene_EntrezID.append(str(geneid))
            result_gene_GOAnno.append(gene_anno_list)
            result_gene_Qualifier.append(gene_anno_qualifier_list)
                    
    print(f'\nThere are {anno_count} annotations used for {gene_count} genes embedding generation')
    
    return result_gene_EntrezID, result_gene_GOAnno, result_gene_Qualifier