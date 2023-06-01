import os
import argparse
import tqdm 
import json 
import time 
import pickle 

import numpy as np
import pandas as pd

import torch
from torch.nn import CrossEntropyLoss

from tokenizers import Tokenizer
from models.progen.modeling_progen import ProGenForCausalLM


########################################################################
# model

def create_model(ckpt, fp16):
    if fp16:
        return ProGenForCausalLM.from_pretrained(ckpt, revision='float16', torch_dtype=torch.float16, low_cpu_mem_usage=True)
    else:
        return ProGenForCausalLM.from_pretrained(ckpt)


def create_tokenizer_custom(file):
    with open(file, 'r') as f:
        return Tokenizer.from_str(f.read())

########################################################################
# fitness

def calc_fitness(model, prots, tokenizer, device='cuda:0', model_context_len=1024, fp16=False, reduction='mean'):
    loss_list = []
    loss_fn = CrossEntropyLoss()
    with torch.no_grad():
        with torch.cuda.amp.autocast(enabled=fp16):
            for prot in tqdm.tqdm(prots):
                loss_val = 0
                sequence_chunks=[]
                if len(prot) <= model_context_len:
                    sequence_chunks = [prot]
                else:
                    len_target_seq = len(prot)
                    num_windows = 1 + int(len_target_seq / model_context_len)
                    start=0
                    for window_index in range(1, num_windows+1):
                        sequence_chunks.append(prot[start:start+model_context_len])
                        start += model_context_len
                
                for chunk in sequence_chunks:
                    is_reversed = False
                    for p in [chunk, chunk[::-1]]:
                        ids = torch.tensor(tokenizer.encode(p).ids).to(device)
                        if len(ids) == 1:
                            continue 
                        input_ids = ids[:-1]
                        targets = ids[1:]
                        logits=model(input_ids).logits
                        # remove terminals
                        bos_token, eos_token = 3, 4
                        if targets[-1] in [bos_token, eos_token]:
                            logits = logits[:-1, ...]
                            targets = targets[:-1]
                        # For the rare case with a protein of length 1025 or 1026, skips the 1 length eos or bos piece on the end  
                        if len(targets) == 0:
                            continue 
                        assert (targets == bos_token).sum() == 0
                        assert (targets == eos_token).sum() == 0

                        # remove unused logits
                        first_token, last_token = 5, 29
                        logits = logits[:, first_token:(last_token+1)]
                        targets = targets - first_token

                        assert logits.shape[1] == (last_token - first_token + 1)
                        is_reversed = True 
                        loss = loss_fn(target=targets.view(-1), input=logits.view(-1,logits.size(-1)))
                        loss_val += -loss.item()
                loss_val /= 2.0 #normalizing for mirroring
                
                if reduction=='mean':
                    loss_val /= len(prot) #average by seq length

                loss_list += [loss_val]
    return np.array(loss_list) 

def get_mutated_sequence(focus_seq, mutant, start_idx=1, AA_vocab="ACDEFGHIKLMNPQRSTVWY"):
    """
    Helper function that mutates an input sequence (focus_seq) via an input mutation triplet (substitutions only).
    Mutation triplet are typically based on 1-indexing: start_idx is used for switching to 0-indexing.
    """
    mutated_seq = list(focus_seq)
    for mutation in mutant.split(":"):
        try:
            from_AA, position, to_AA = mutation[0], int(mutation[1:-1]), mutation[-1]
        except:
            print("Issue with mutant: "+str(mutation))
        relative_position = position - start_idx
        assert (from_AA==focus_seq[relative_position]), "Invalid from_AA or mutant position: "+str(mutation)+" from_AA: "+str(from_AA) + " relative pos: "+str(relative_position) + " focus_seq: "+str(focus_seq)
        assert (to_AA in AA_vocab) , "Mutant to_AA is invalid: "+str(mutation)
        mutated_seq[relative_position] = to_AA
    return "1"+"".join(mutated_seq)+"2"

def main():
    """
    Main script to score sets of mutated protein sequences (substitutions or indels) with Progen2.
    """
    starttime = time.time()
    # (0) constants
    models_151M = [ 'progen2-small' ]
    models_754M = [ 'progen2-medium', 'progen2-oas', 'progen2-base' ]
    models_2B = [ 'progen2-large', 'progen2-BFD90' ]
    models_6B = [ 'progen2-xlarge' ]
    models = models_151M + models_754M + models_2B + models_6B

    parser = argparse.ArgumentParser(description='Tranception scoring')
    parser.add_argument('--checkpoint', default="/n/groups/marks/projects/marks_lab_and_oatml/protein_transformer/baseline_models/progen2/progen2-small", type=str, help='Name of or path to Progen2 model')
    parser.add_argument('--dataset_reference_file', default='/home/pn73/Tranception/proteingym/ProteinGym_reference_file_substitutions.csv', type=str, help='Path of DMS folder')
    parser.add_argument('--dataset_folder', default='/n/groups/marks/projects/marks_lab_and_oatml/protein_transformer/Tranception_open_source/DMS_files/ProteinGym_substitutions', type=str, help='Path of DMS folder')
    parser.add_argument('--target_sequence_index', type=int, default=None, help='Path of DMS folder')
    parser.add_argument('--target_sequence_index_range_start', type=int, default=None)
    parser.add_argument('--target_sequence_index_range_end', type=int, default=None)
    parser.add_argument('--output_scores_folder', default=None, type=str, help='Name of folder to write model scores to')
    parser.add_argument('--indel_mode', action='store_true', help='Whether to score sequences with insertions and deletions')
    parser.add_argument("--score_wt_only", action="store_true", help="Pass this to score only the wildtype sequences in the mapping file")
    parser.add_argument('--mutated_sequence_column', default='mutant', type=str, help='Name of column with mutated protein sequence')
    parser.add_argument('--fp16', action='store_true', help='Whether to score sequences with half precision')
    parser.add_argument("--tokenizer_file", default="tokenizer.json", type=str)
    parser.add_argument("--label_column", default="label", type=str, help="Name of column with label")
    args = parser.parse_args()

    if not os.path.isdir(args.output_scores_folder):
        os.makedirs(args.output_scores_folder, exist_ok=True)

    model = create_model(ckpt=args.checkpoint, fp16=args.fp16).cuda()
    config = json.load(open(args.checkpoint+os.sep+'config.json',"r"))
    print("Maximum context length: {}".format(config['n_positions']))
    tokenizer = create_tokenizer_custom(file=args.tokenizer_file)

    target_seq_mapfile = pd.read_csv(args.dataset_reference_file)
    list_target_sequences = target_seq_mapfile["target_sequence_id"].tolist()
    if args.target_sequence_index != None:
        target_sequence_ids = [list_target_sequences[args.target_sequence_index]]
        target_sequence_file_names = [target_seq_mapfile["target_sequence_filename"][target_seq_mapfile["target_sequence_id"]==target_sequence_ids[0]].values[0]]
        target_seqs = [target_seq_mapfile["target_seq"][target_seq_mapfile["target_sequence_id"]==target_sequence_ids[0]].values[0].upper()]
        if "mutated_sequence_column" in target_seq_mapfile.columns:
            mutated_sequence_columns = [target_seq_mapfile["mutated_sequence_column"][args.target_sequence_index]]
        else:
            mutated_sequence_columns = [args.mutated_sequence_column]
    elif args.target_sequence_index_range_start != None and args.target_sequence_index_range_end != None: 
        target_sequence_ids = list_target_sequences[args.target_sequence_index_range_start:args.target_sequence_index_range_end]
        target_sequence_file_names = [target_seq_mapfile["target_sequence_filename"][target_seq_mapfile["target_sequence_id"]==target_sequence_id].values[0] for target_sequence_id in target_sequence_ids]
        target_seqs = [target_seq_mapfile["target_seq"][target_seq_mapfile["target_sequence_id"]==target_sequence_id].values[0].upper() for target_sequence_id in target_sequence_ids]
        if "mutated_sequence_column" in target_seq_mapfile.columns:
            mutated_sequence_columns = [target_seq_mapfile["mutated_sequence_column"][i] for i in range(args.target_sequence_index_range_start, args.target_sequence_index_range_end)]
        else:
            mutated_sequence_columns = [args.mutated_sequence_column for i in range(args.target_sequence_index_range_start, args.target_sequence_index_range_end)]

    print("Computing scores for: {} with Progen2: {}".format(target_sequence_ids, args.Progen2_model_name_or_path))
    for i in range(len(target_sequence_ids)):
        print(target_sequence_ids[i])
        if not args.score_wt_only:
            variant_data = pd.read_csv(args.dataset_folder + os.sep + target_sequence_file_names[i], low_memory=False)
            if mutated_sequence_columns[i] not in variant_data.columns:
                print(f"{mutated_sequence_columns[i]} not present in dataframe")
                exit()
            else:
                mutant_column = mutated_sequence_columns[i]
            variant_data["mutated_sequence_processed"] = variant_data[mutant_column].apply(lambda x: "1" + x + "2")
            model_scores = calc_fitness(model=model, prots=np.array(variant_data['mutated_sequence_processed']), model_context_len=int(config['n_positions']), tokenizer=tokenizer, fp16=args.fp16)
            focus_seq_score = calc_fitness(model=model, prots=np.array("1" + target_seqs[i] + "2"), model_context_len=int(config['n_positions']), tokenizer=tokenizer, fp16=args.fp16)
            variant_data['Progen2_score']=model_scores
            variant_data['Progen2_score_wt_ratio']= model_scores - focus_seq_score[0]
            scoring_filename = args.output_scores_folder+os.sep+target_sequence_ids[i]+'.csv'
            variant_data["model_score"] = variant_data["Progen2_score_wt_ratio"]
            if args.label_column in variant_data.columns:
                variant_data = variant_data.rename(columns={args.label_column:"label", mutant_column:"mutated_sequence"})
                variant_data[["mutated_sequence",'Progen2_score','Progen2_score_wt_ratio',"label"]].to_csv(scoring_filename, index=False)
            else:
                variant_data = variant_data.rename(columns={mutant_column:"mutated_sequence"})
                variant_data[["mutated_sequence",'Progen2_score','Progen2_score_wt_ratio']].to_csv(scoring_filename, index=False)
        else:
            model_scores = calc_fitness(model=model, prots=np.array(["1" + target_seqs[i] + "2"]), model_context_len=int(config["n_positions"]), tokenizer=tokenizer, fp16=args.fp16)
            result_dict = {"mutated_sequence":[target_seqs[i]], "Progen2_score":model_scores, "Progen2_score_wt_ratio":0.0}
            result_df = pd.DataFrame(result_dict)
            scoring_filename = args.output_scores_folder+os.sep+target_sequence_ids[i]+'.csv'
            result_df.to_csv(scoring_filename, index=False)
        print(f"Time taken to run:{time.time() - starttime}")

if __name__ == '__main__':
    main()
