import pandas as pd 

if __name__ == "__main__":

    df = pd.read_csv("./ASD_vep_annotated_GRCh38.tsv", comment="#", sep="\t")
    print(df.columns)
    print(df["Feature"].value_counts())
    print(df[~(df["Feature"].str.contains("X") | df["Feature"].str.contains("NR"))]["Feature"].value_counts())
    # print(df["ENSP"].value_counts())
    # print(df["HGVSp"].value_counts()) 
    # # Removing variants with NaN HGVS annotations 
    # df["protein"] = df["ENSP"] 
    # diffs = df["HGVSp"].str.split(":p.").str[0] != df["protein"]
    # print(f"Number of differences between HGVSp and protein: {diffs.sum()}")
    # df_inframe_mapped_isoforms_valid = df[~diffs].copy()
    # df_inframe_mapped_isoforms_valid = df_inframe_mapped_isoforms_valid[df_inframe_mapped_isoforms_valid["protein"] != "-"]

    # # Dropping Terminating mutations 
    # df_inframe_mapped_clean = df_inframe_mapped_isoforms_valid[~df_inframe_mapped_isoforms_valid["HGVSp"].str.contains("Ter")].copy()

    # df = df_inframe_mapped_clean[df_inframe_mapped_clean["Consequence"].str.contains("inframe")]
    # df = df.drop_duplicates("protein")
    # df = df[~df["protein"].str.contains("X")]
    # print(df)