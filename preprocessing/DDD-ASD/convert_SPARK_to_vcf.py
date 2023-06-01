import pandas as pd 
import os 
import argparse 

if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description="Convert SPARK data from tsv to vcf format")
    parser.add_argument("--input_tsv_file", type=str, help="Path to input tsv file")
    parser.add_argument("--output_file", type=str, help="Path to output vcf file")
    args = parser.parse_args()

    df = pd.read_csv(args.input_tsv_file, sep="\t")
    df = df.sort_values(by=["Chrom", "Position"])
    with open(args.output_file, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##INFO=<ID=ID,Number=1,Type=String,Description='Input ID'>\n")
        f.write("##contig=<ID=1>\n")
        f.write("##contig=<ID=2>\n")
        f.write("##contig=<ID=3>\n")
        f.write("##contig=<ID=4>\n")
        f.write("##contig=<ID=5>\n")
        f.write("##contig=<ID=6>\n")
        f.write("##contig=<ID=7>\n")
        f.write("##contig=<ID=8>\n")
        f.write("##contig=<ID=9>\n")
        f.write("##contig=<ID=10>\n")
        f.write("##contig=<ID=11>\n")
        f.write("##contig=<ID=12>\n")
        f.write("##contig=<ID=13>\n")
        f.write("##contig=<ID=14>\n")
        f.write("##contig=<ID=15>\n")
        f.write("##contig=<ID=16>\n")
        f.write("##contig=<ID=17>\n")
        f.write("##contig=<ID=18>\n")
        f.write("##contig=<ID=19>\n")
        f.write("##contig=<ID=20>\n")
        f.write("##contig=<ID=21>\n")
        f.write("##contig=<ID=22>\n")
        f.write("##contig=<ID=X>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i,row in df.iterrows():
            varid = row["VarID"].replace(":","_")
            # Replacing last ocurrence of : with / to match vcf format
            varid = varid[::-1].replace("_","/",1)[::-1]
            f.write(f"{row['Chrom']}\t{row['Position']}\t{varid}\t{row['Ref']}\t{row['Alt']}\t.\t.\tID={row['Pheno']}\n")