import argparse 
import os 
import pandas as pd 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert DDD data from text to vcf format") 
    parser.add_argument("--input_text_file", type=str, help="Path to input text file")
    parser.add_argument("--output_file", type=str, help="Path to output vcf file")
    args = parser.parse_args()

    df = pd.read_csv(args.input_text_file, sep="\t") 
    df = df.sort_values(by=["chrom", "pos"])
    with open(args.output_file, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##fileDate=20230523\n")
        f.write("##source=DenovoWEST\n")
        f.write("##reference=GRCh37\n")
        f.write("##phasing=none\n")
        f.write("##contig=<ID=1>\n")
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
        f.write("##contig=<ID=2>\n")
        f.write("##contig=<ID=20>\n")
        f.write("##contig=<ID=21>\n")
        f.write("##contig=<ID=22>\n")
        f.write("##contig=<ID=3>\n")
        f.write("##contig=<ID=4>\n")
        f.write("##contig=<ID=5>\n")
        f.write("##contig=<ID=6>\n")
        f.write("##contig=<ID=7>\n")
        f.write("##contig=<ID=8>\n")
        f.write("##contig=<ID=9>\n")
        f.write("##contig=<ID=X>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i, row in df.iterrows():
             f.write(f"{row['chrom']}\t{row['pos']}\t.\t{row['ref']}\t{row['alt']}\t.\t.\t.\n")