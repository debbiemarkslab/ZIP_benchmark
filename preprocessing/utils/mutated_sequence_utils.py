hgvs_regex = " \(p\.(.+)\)$" # alskdfjlsakjf (p.<>)
refseq_regex = "^(NM_\d+.\d+)\("  # Optionally add .+ at the end to swallow rest

def insert(mut, SEQ):
    "Using single-char AA notation and hgvs notation"
    assert "_" in mut
    # e.g. S10_I11insG
    start_stop, inserted_aas = mut.split("ins")
    # Check start and end AAs with wt, then insert
    start_aa_pos, end_aa_pos = start_stop.split("_")
    start_aa, start_pos = start_aa_pos[0], int(start_aa_pos[1:])
    assert SEQ[start_pos-1] == start_aa

    end_aa, end_pos = end_aa_pos[0], int(end_aa_pos[1:])
    assert SEQ[end_pos-1] == end_aa

    # Check the insertion is between two contiguous AAs
    assert end_pos - start_pos == 1

    # Insert the aas
    # Start position and end positions are inclusive and retained
    assert inserted_aas.isalpha()
    # TODO check all valid amino acids with Bio.IUPAC dict
    
    # Both flanking positions are inclusive
    mutated_seq = SEQ[:start_pos+1-1] + inserted_aas + SEQ[end_pos-1:]
    assert len(mutated_seq) == len(SEQ) + len(inserted_aas)
    return mutated_seq

assert insert("S1_I2insG", "SIE") == "SGIE"
assert insert("S1_I2insGGS", "SIE") == "SGGSIE"

def delete(mut, SEQ):
    # Similar to insertion, except now the positions may be > 1
    if mut.endswith("del"):
        mut = mut[:-3]
    # mut = mut.removesuffix("del")
    if "_" not in mut:
        # Single deletion e.g. S10del
        aa, pos = mut[0], int(mut[1:])
        assert SEQ[pos-1] == aa, f"{SEQ[pos-1]=} != {aa=}; {pos=}, {mut=}"
        mutated_seq = SEQ[:pos-1] + SEQ[pos+1-1:]
        assert len(mutated_seq) == len(SEQ) - 1
        return mutated_seq
    else:
        # e.g. S10_I11del deletes both S10 and I11
        # Check start and end AAs with wt, then insert
        start_aa_pos, end_aa_pos = mut.split("_")
        start_aa, start_pos = start_aa_pos[0], int(start_aa_pos[1:])
        assert SEQ[start_pos-1] == start_aa

        end_aa, end_pos = end_aa_pos[0], int(end_aa_pos[1:])
        assert SEQ[end_pos-1] == end_aa

        mutated_seq = SEQ[:start_pos-1] + SEQ[end_pos+1-1:]  # Both positions exclusive as they're both deleted
        assert len(mutated_seq) == len(SEQ) - (end_pos-start_pos+1)
        return mutated_seq

assert delete("B2del", "ABCDE") == "ACDE", delete("B2del", "ABCDE")
assert delete("C3_D4del", "ABCDE") == "ABE"
assert delete("B2_E5del", "ABCDE") == "A"

def delins(mut, SEQ):
    # Following https://varnomen.hgvs.org/recommendations/protein/variant/delins/
    # A combination of the insert() and delete() functions
    
    # region = position amino acid or range of amino acids deleted = Arg123_Lys127
    region, inserted_aas = mut.split("delins")
    assert inserted_aas.isalpha()
    
    # Checking mutant is valid (matches sequence)
    if "_" not in region:
        # Single deletion e.g. S10insdel
        start_aa, start_pos = region[0], int(region[1:])
        end_pos = start_pos
        assert SEQ[start_pos-1] == start_aa, f"{SEQ[start_pos-1]=} != {start_aa=}; {start_pos=}, {mut=}"
    else:
        start_aa_pos, end_aa_pos = region.split("_")
        start_aa, start_pos = start_aa_pos[0], int(start_aa_pos[1:])
        assert SEQ[start_pos-1] == start_aa

        end_aa, end_pos = end_aa_pos[0], int(end_aa_pos[1:])
        assert SEQ[end_pos-1] == end_aa
    
    # Both positions exclusive as they're both deleted (recall that in the case of only one aa, start_pos=end_pos so only it is deleted)
    mutated_seq = SEQ[:start_pos-1] + inserted_aas + SEQ[end_pos+1-1:]  
    assert len(mutated_seq) == len(SEQ) - (end_pos-start_pos+1) + len(inserted_aas)
    return mutated_seq

assert delins("C2delinsTV", SEQ="ACDE") == "ATVDE", delins("C2delinsTV", SEQ="ACDE")
assert delins("C2_L3delinsT", SEQ="ACLDE") == "ATDE"
assert delins("C2_L3delinsT", SEQ="ACL") == "AT"

def duplicate(mut, SEQ):
    # Following http://varnomen.hgvs.org/recommendations/protein/variant/duplication/
    # Positions are inclusive e.g. S10_L11dup -> S10, L11, S12, L13
    assert mut.endswith("dup")
    # region = position amino acid or range of amino acids deleted = Arg123_Lys127
    if mut.endswith("dup"):
        mut = mut[:-3]
    # mut = mut.removesuffix("dup")
    
    # Checking mutant is valid (matches sequence)
    if "_" not in mut:
        # Single deletion e.g. S10dup
        start_aa, start_pos = mut[0], int(mut[1:])
        
        assert SEQ[start_pos-1] == start_aa, f"{SEQ[start_pos-1]=} != {start_aa=}; {start_pos=}, {mut=}"
        # start_pos is inclusive
        mutated_seq = SEQ[:start_pos+1-1] + start_aa + SEQ[start_pos+1-1:]
        assert len(mutated_seq) == len(SEQ) + 1, mutated_seq
        return mutated_seq

    else:
        start_aa_pos, end_aa_pos = mut.split("_")
        start_aa, start_pos = start_aa_pos[0], int(start_aa_pos[1:])
        assert SEQ[start_pos-1] == start_aa

        end_aa, end_pos = end_aa_pos[0], int(end_aa_pos[1:])
        assert SEQ[end_pos-1] == end_aa
        # Positions are inclusive and retained e.g. S10_L11dup -> S10, L11, S12, L13
        inserted_aas = SEQ[start_pos-1:end_pos+1-1]
        
        # Check again, we can't be too careful
        assert len(inserted_aas) == (end_pos-start_pos)+1, inserted_aas
        assert inserted_aas[0] == start_aa
        assert inserted_aas[-1] == end_aa  
        
        # Don't insert between span, insert after original substring of positions
        mutated_seq = SEQ[:end_pos+1-1] + inserted_aas + SEQ[end_pos+1-1:]  
        assert len(mutated_seq) == len(SEQ) + len(inserted_aas), mutated_seq
        return mutated_seq

# From reference site:
assert duplicate("A3dup", SEQ="MGARSSH") == "MGAARSSH"
assert duplicate("A3_S5dup", SEQ="MGARSSH") == "MGARSARSSH"
# Duplicate all
assert duplicate("A1_S3dup", SEQ="ACS") == "ACSACS"
# Duplicate last
assert duplicate("C3dup", SEQ="ABC") == "ABCC"
# Duplicate first
assert duplicate("A1dup", SEQ="ABC") == "AABC"
