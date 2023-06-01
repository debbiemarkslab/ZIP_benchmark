import pandas as pd
import numpy as np

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
    return "".join(mutated_seq)

def DMS_file_cleanup(DMS_filename, target_seq, start_idx=1, end_idx=None, DMS_mutant_column='mutant', DMS_phenotype_name='score', DMS_directionality=1, AA_vocab = "ACDEFGHIKLMNPQRSTVWY"):
    """
    Function to process the raw indel DMS assay data (eg., removing invalid mutants, aggregate silent mutations).
    """
    DMS_data = pd.read_csv(DMS_filename, low_memory=False)
    end_idx = start_idx + len(target_seq) - 1 if end_idx is None else end_idx
    DMS_data['mutated_sequence'] = DMS_data[DMS_mutant_column]
    
    DMS_data=DMS_data[DMS_data['mutated_sequence'].notnull()].copy()
    
    DMS_data[DMS_phenotype_name]=pd.to_numeric(DMS_data[DMS_phenotype_name],errors='coerce')
    DMS_data=DMS_data[np.isfinite(DMS_data[DMS_phenotype_name])]
    DMS_data.dropna(subset = [DMS_phenotype_name], inplace=True)
    DMS_data['label'] = DMS_data[DMS_phenotype_name] * DMS_directionality
    DMS_data=DMS_data[['mutated_sequence','label']]
    DMS_data=DMS_data.groupby('mutated_sequence').mean().reset_index()
    DMS_data=DMS_data[['mutated_sequence','label']]
    return DMS_data

