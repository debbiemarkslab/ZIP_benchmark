import argparse 
import pandas as pd 


"""
This script takes in the raw DDD file and uses the given index file to add columns 
for gene, HGVSp, refseq id and protein consequence. 
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add columns for gene, HGVSp, refseq id and protein consequence to the raw DDD file.")
    parser.add_argument("--raw_DDD_file", help="The raw DDD file.")
    parser.add_argument("--raw_ASD_file", help="The raw ASD file.")
    parser.add_argument("--raw_ASD_replicate_file", help="The raw file for ASD replicate data")
    parser.add_argument("--index_file", help="The index file.")
    parser.add_argument("--output_file", help="The output file.")
    args = parser.parse_args()

    index_df = pd.read_csv(args.index_file, sep="\t")

    ddd_df = pd.read_csv(args.raw_DDD_file, sep="\t")
    ddd_df = ddd_df.drop(["study","consequence"], axis="columns")
    ddd_df = ddd_df.merge(index_df, on=["chrom","pos","ref","alt"], how="left")
    ddd_df = ddd_df[["chrom","pos","ref","alt","gene","HGVSp","protein","Consequence", "study"]]

    asd_df = pd.read_csv(args.raw_ASD_file, sep="\t")
    asd_df = asd_df.rename(columns={"Chrom":"chrom","Position":"pos","Ref":"ref","Alt":"alt"})
    asd_df = asd_df[["chrom","pos","ref","alt"]]
    asd_df = asd_df.merge(index_df, on=["chrom","pos","ref","alt"], how="left")
    asd_df = asd_df[["chrom","pos","ref","alt","gene","HGVSp","protein","Consequence", "study"]]

    asd_rep_df = pd.read_csv(args.raw_ASD_replicate_file, sep="\t")
    asd_rep_df = asd_rep_df.rename(columns={"Chrom":"chrom","Position":"pos","Ref":"ref","Alt":"alt"})
    asd_rep_df = asd_rep_df.merge(index_df, on=["chrom","pos","ref","alt"], how="left")
    asd_rep_df = asd_rep_df[["chrom","pos","ref","alt","gene","HGVSp","protein","Consequence", "study"]]

    # merge the two dataframes
    df = pd.concat([ddd_df, asd_df, asd_rep_df], axis=0)
    df.to_csv(args.output_file, sep="\t", index=False)