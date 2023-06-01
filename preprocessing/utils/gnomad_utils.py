def drop_invalid_hgvsp_and_add_variants(df):
    df = df[(df["BIOTYPE"] == "protein_coding")]
    # print(f"Number of protein consequences for protein coding gene variants: {df.index.nunique()}")
    df = df.dropna(subset="HGVSp")
    df["protein_variant"] = df["HGVSp"].str.extract("p\.([^\)]+)")
    df = df.dropna(subset="protein_variant")
    # print(f"Total number of protein variants with HGVSp annotations: {len(df)}")
    # print(f"Number of genetic variants with associated protein consequences that have HGVSp annotations: {df.index.nunique()}")
    return df 

def add_frequency_annotations(df, rare_threshold=0.01, common_threshold=0.05):
    def freq_check(allele_freq):
        if allele_freq < rare_threshold:
            return "rare" 
        elif allele_freq > common_threshold:
            return "common"
        else:
            return "neutral"
    df["frequency"] = df["AF_popmax"].apply(freq_check)

def add_mutation_type_annotations(df):
    df["is_frameshift"] = df["protein_variant"].str.contains("fs")
    df["protein_variant"] = df["protein_variant"].str.replace("%3D","=",regex=False)
    df["inframe_synon"] = df["protein_variant"].str.contains("=",regex=False)
    df["inframe_delins"] = df["protein_variant"].str.contains("delins")
    df["inframe_ins"] = df["protein_variant"].str.contains("ins") & ~df["inframe_delins"]
    df["inframe_del"] = df['protein_variant'].str.endswith("del")  # exclude mutants such as p.Ter117delextTer?
    df["inframe_dup"] = df["protein_variant"].str.contains("dup")
    df["terminating"] = df["protein_variant"].str.contains("Ter")
    df["unknown"] = df["protein_variant"].str.contains("\?")
    df['inframe_single_sub'] = df["protein_variant"].str.contains(".*[A-Z][a-z]{2}\d+[A-Z][a-z]{2}")
    # Cast all to bool
    df['inframe_delins'] = df['inframe_delins'].astype(bool)
    df['inframe_ins'] = df['inframe_ins'].astype(bool)
    df['inframe_del'] = df['inframe_del'].astype(bool)
    df['inframe_dup'] = df['inframe_dup'].astype(bool)
    
def get_duplicate_depth(mut, three_letter_abr=True):
    # Following http://varnomen.hgvs.org/recommendations/protein/variant/duplication/
    # Positions are inclusive e.g. S10_L11dup -> S10, L11, S12, L13
    assert mut.endswith("dup")
    # region = position amino acid or range of amino acids deleted = Arg123_Lys127
    mut = mut[:-3]
    # Checking mutant is valid (matches sequence)
    if "_" not in mut:
        return 1
    else:
        start_aa_pos, end_aa_pos = mut.split("_")
        if three_letter_abr:
            start_aa, start_pos = start_aa_pos[:3], int(start_aa_pos[3:])
            end_aa, end_pos = end_aa_pos[:3], int(end_aa_pos[3:])
        else:
            start_aa, start_pos = start_aa_pos[0], int(start_aa_pos[1:])
            end_aa, end_pos = end_aa_pos[0], int(end_aa_pos[1:])
        depth = end_pos - start_pos + 1 
        return depth 
    
def get_delins_depth(mut, three_letter_abr=True):
    # Following https://varnomen.hgvs.org/recommendations/protein/variant/delins/
    # A combination of the insert() and delete() functions
    
    # region = position amino acid or range of amino acids deleted = Arg123_Lys127
    region, inserted_aas = mut.split("delins")
    assert inserted_aas.isalpha()
    # Checking mutant is valid (matches sequence)
    if "_" not in region:
        # Single deletion e.g. S10insdel
        if three_letter_abr:
            start_aa, start_pos = region[:3], int(region[3:])
        else:
            start_aa, start_pos = region[0], int(region[1:])
        end_pos = start_pos
    else:
        start_aa_pos, end_aa_pos = region.split("_")
        if three_letter_abr:
            start_aa, start_pos = start_aa_pos[:3], int(start_aa_pos[3:])
            end_aa, end_pos = end_aa_pos[:3], int(end_aa_pos[3:])
        else:
            start_aa, start_pos = start_aa_pos[0], int(start_aa_pos[1:])
            end_aa, end_pos = end_aa_pos[0], int(end_aa_pos[1:])
    if three_letter_abr:
        aa_depth = len(inserted_aas)/3 
    else:
        aa_depth = len(inserted_aas)
    depth = max(start_pos - end_pos, aa_depth)     
    # Both positions exclusive as they're both deleted (recall that in the case of only one aa, start_pos=end_pos so only it is deleted)
    return depth

def get_deletion_depth(mut, three_letter_abr=True):
    # Similar to insertion, except now the positions may be > 1
    mut = mut[:-3]
    if "_" not in mut:
        return 1
    else:
        # e.g. S10_I11del deletes both S10 and I11
        # Check start and end AAs with wt, then insert
        start_aa_pos, end_aa_pos = mut.split("_")
        if three_letter_abr:
            start_aa, start_pos = start_aa_pos[:3], int(start_aa_pos[3:])
            end_aa, end_pos = end_aa_pos[:3], int(end_aa_pos[3:])
        else:
            start_aa, start_pos = start_aa_pos[0], int(start_aa_pos[1:])
            end_aa, end_pos = end_aa_pos[0], int(end_aa_pos[1:])
        return end_pos - start_pos 
    
# Getting mutation depth for gnomad variants (frequency of different kinds of mutations at different depths and frequency of mutations of a particular depth)
def get_insertion_depth(mut, three_letter_abr=True):
    "Using single-char AA notation and hgvs notation"
    # e.g. S10_I11insG
    start_stop, inserted_aas = mut.split("ins")
    # Check start and end AAs with wt, then insert
    start_aa_pos, end_aa_pos = start_stop.split("_")
    if three_letter_abr:
        start_aa, start_pos = start_aa_pos[:3], int(start_aa_pos[3:])
        end_aa, end_pos = end_aa_pos[:3], int(end_aa_pos[3:])
    else:
        start_aa, start_pos = start_aa_pos[0], int(start_aa_pos[1:])
        end_aa, end_pos = end_aa_pos[0], int(end_aa_pos[1:])
    # Check the insertion is between two contiguous AAs
    assert end_pos - start_pos == 1

    # Insert the aas
    # Start position and end positions are inclusive and retained
    assert inserted_aas.isalpha()
    # TODO check all valid amino acids with Bio.IUPAC dict
    if three_letter_abr:
        num_aas = len(inserted_aas)/3 
    else:
        num_aas = len(inserted_aas)
    return num_aas

def add_depth_annotation(row, three_letter_abr=True):
    if row["inframe_dup"]:
        depth = get_duplicate_depth(row["protein_variant"], three_letter_abr=three_letter_abr)
    elif row["inframe_ins"]:
        depth = get_insertion_depth(row["protein_variant"], three_letter_abr=three_letter_abr)
    elif row["inframe_del"]:
        depth = get_deletion_depth(row["protein_variant"])
    elif row["inframe_delins"]:
        depth = get_delins_depth(row["protein_variant"], three_letter_abr=three_letter_abr)
    elif row["inframe_synon"]:
        depth = 0
    elif row["inframe_single_sub"]:
        depth = 1
    else:
        depth = None 
    return depth 

def strip_gnomad_headers(in_file, out_file):
    from tqdm import tqdm
    import os
    import gzip
    t = tqdm(total=os.path.getsize(in_file), unit="B", unit_scale=True)
    with open(in_file, "rb") as f, gzip.open(out_file, "wt") as f_out:
        g = gzip.GzipFile(fileobj=f)
        last_pos = f.tell()
        for line in g:
            t.update(f.tell() - last_pos)
            last_pos = f.tell()
            line = line.decode()
            
            val = -1
            for i in range(0, 7):
                val = line.find("\t", val+1)
            pos_vep = line.find(";vep", val+1)
            pos_ac = line.find("AC=", val)
            pos_ac_end = line.find(";",pos_ac+1)
            pos_ac_raw = line.find(";AC_raw=", val+1)
            pos_ac_raw_end = line.find(";",pos_ac_raw+1)
            pos_ac_popmax = line.find(";AC_popmax=", val+1)
            pos_ac_popmax_end = line.find(";",pos_ac_popmax+1)
            pos_af = line.find(";AF=", val+1)
            pos_af_end = line.find(";",pos_af+1)
            pos_af_raw = line.find(";AF_raw=", val+1)
            pos_af_raw_end = line.find(";",pos_af_raw+1)
            pos_af_popmax = line.find(";AF_popmax=", val+1)
            pos_af_popmax_end = line.find(";",pos_af_popmax+1)
            if pos_vep == -1 or pos_ac_raw == -1 or pos_af_raw == -1:  # ignore missing pos_af and pos_af_popmax as some lines are missing an AF or AF_popmax annotation
                f_out.write(line)
            else:
                val_ac = "" if pos_ac == -1 else line[pos_ac:pos_ac_end+1]
                val_ac_raw = line[pos_ac_raw+1:pos_ac_raw_end+1]
                val_ac_popmax = "" if pos_ac_popmax == -1 else line[pos_ac_popmax+1:pos_ac_popmax_end+1]
                val_af = "" if pos_af == -1 else line[pos_af+1:pos_af_end+1]
                val_af_raw = line[pos_af_raw+1:pos_af_raw_end+1]
                val_af_popmax = "" if pos_af_popmax == -1 else line[pos_af_popmax+1:pos_af_popmax_end+1]
                line = line[:val] + "\t" + val_ac + val_ac_raw + val_ac_popmax + val_af + val_af_raw + val_af_popmax + line[pos_vep+1:]
                f_out.write(line)
    t.close()

if __name__ == "__main__":
    gnomad_latest_file = "/n/groups/marks/databases/gnomad/v2_exomes_GRCh38/latest_gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz"
    out_file = "/n/groups/marks/databases/gnomad/v2_exomes_GRCh38/latest_gnomad.exomes.r2.1.1.sites.liftover_grch38_AC_AF_vep_only.vcf.bgz"
    strip_gnomad_headers(gnomad_latest_file, out_file)

