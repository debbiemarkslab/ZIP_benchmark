import pandas as pd 
from tqdm import tqdm 
import os 

def create_mapfile(variant_folder, output_file, alignment_folder=None, protein_sequence_column="protein_sequence", mutant_column="mutated_sequence"):
    mapping_df = {"target_sequence_id":[],"target_sequence_filename":[],"target_seq":[],"MSA_filename":[],"weight_file_name":[],"MSA_start":[],"MSA_end":[],"mutated_sequence_column":[]} 
    print("Generating mapping file") 
    for file in tqdm(os.listdir(variant_folder)):
        if not file.endswith(".csv"):
            continue 
        name = os.path.splitext(file)[0]
        df = pd.read_csv(variant_folder + os.sep + file)
        mapping_df["target_sequence_id"].append(name)
        mapping_df["target_sequence_filename"].append(file)
        mapping_df["target_seq"].append(df[protein_sequence_column].values[0])
        mapping_df["mutated_sequence_column"].append(mutant_column)
        if alignment_folder is not None:
            if os.path.exists(f"{alignment_folder}/{name}.a2m"):
                mapping_df["MSA_filename"].append(f"{name}.a2m")
                mapping_df["weight_file_name"].append(None)
                mapping_df["MSA_start"].append(1)
                mapping_df["MSA_end"].append(len(df[protein_sequence_column].values[0]))
            else:
                print(f"Warning: alignment file {name}.a2m not found in {alignment_folder}. Skipping alignment information for this file.")
                mapping_df["MSA_filename"].append(None)
                mapping_df["weight_file_name"].append(None)
                mapping_df["MSA_start"].append(None)
                mapping_df["MSA_end"].append(None)
        else:
            mapping_df["MSA_filename"].append(None)
            mapping_df["weight_file_name"].append(None)
            mapping_df["MSA_start"].append(None)
            mapping_df["MSA_end"].append(None)
    output_folder = os.path.dirname(output_file)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    mapping_df = pd.DataFrame(mapping_df)
    mapping_df.to_csv(output_file,index=False)
