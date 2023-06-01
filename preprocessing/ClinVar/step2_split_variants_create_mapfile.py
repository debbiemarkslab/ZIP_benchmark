import argparse 
import pandas as pd 
import os 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split ClinVar data and create mapfile")
    parser.add_argument("--input_file", type=str, help="Path to processed ClinVar file")
    parser.add_argument("--alignment_folder", type=str, help="Path to alignment folder")
    parser.add_argument("--output_folder", type=str, help="Path to output folder")
    parser.add_argument("--mapping_file", type=str, help="Path to save mapping file to")
    args = parser.parse_args()

    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    df = pd.read_csv(args.input_file)

    # making a unique refseq id column (combining the few cases where there are two refseq ids for one protein)
    def make_refseq_unique(group):
        assert len(group["protein_sequence"].unique()) == 1
        if len(group["refseq_protein_id"].unique()) != 1:
            group["refseq_unique_id"] = "-".join(group["refseq_protein_id"].unique().tolist())
        else:
            group["refseq_unique_id"] = group["refseq_protein_id"]
        return group
    df_unique_refseq = df.groupby("protein_sequence", group_keys=True).apply(make_refseq_unique)

    counter = 0
    unique_refseq_to_protein = []

    for refseq_protein_id in df_unique_refseq["refseq_unique_id"].unique():
        subset = df_unique_refseq[df_unique_refseq["refseq_unique_id"] == refseq_protein_id]
        protein_sequences = subset["protein_sequence"].unique()
        assert len(protein_sequences) == 1, "There should only be one unique protein_sequence per refseq id"
        
        counter += len(subset)
        subset = subset.rename(columns={"protein_sequence_mutated":"mutated_sequence"})
        subset.to_csv(os.path.join(args.output_folder, f"{refseq_protein_id}.csv"), index=False)
        # Save RefSeq - protein mapping
        seq = protein_sequences[0]
        unique_refseq_to_protein.append({"protein_sequence": seq, "refseq_unique_id": refseq_protein_id})

    print(f"{counter} variants processed across {df_unique_refseq['refseq_unique_id'].nunique()} unique refseq ids")
    df_refseq_to_protein = pd.DataFrame(unique_refseq_to_protein)

    # Making mapping file 
    # TODO: Adjust this to use the mapfile function in utils 
    # Add other fields in ProteinGym format:
    df_tranception_mapping = pd.DataFrame()
    df_tranception_mapping["target_sequence_id"] = df_refseq_to_protein["refseq_unique_id"]
    df_tranception_mapping["target_sequence_filename"] = df_refseq_to_protein["refseq_unique_id"] + ".csv"
    df_tranception_mapping["target_seq"] = df_refseq_to_protein["protein_sequence"]
    df_tranception_mapping["mutated_sequence_column"] = "mutated_sequence"
    df_tranception_mapping["MSA_filename"] = df_refseq_to_protein["refseq_unique_id"] + ".a2m"
    for id_val in df_refseq_to_protein["refseq_unique_id"]:
        filenames = os.listdir(args.alignment_folder) 
        assert id_val + ".a2m" in filenames 
    df_tranception_mapping["weight_file_name"] = None 
    df_tranception_mapping["MSA_start"] = 1 
    df_tranception_mapping["MSA_end"] = df_tranception_mapping["target_seq"].apply(len)
    df_tranception_mapping = df_tranception_mapping.sort_values(by="target_sequence_id")

    if not os.path.exists(os.path.dirname(args.mapping_file)):
        os.makedirs(os.path.dirname(args.mapping_file))
    # Save to CSV
    df_tranception_mapping.to_csv(args.mapping_file, index=False)