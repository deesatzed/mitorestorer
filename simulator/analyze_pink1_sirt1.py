import math

def parse_cif_atoms(cif_path):
    atoms = []
    with open(cif_path, "r") as f:
        lines = f.readlines()
    in_atom = False
    headers = []
    for line in lines:
        s = line.strip()
        if s == "loop_":
            in_atom = False
            headers = []
            continue
        if s.startswith("_atom_site."):
            in_atom = True
            headers.append(s.split(".")[1])
            continue
        if in_atom:
            if s.startswith("_") or s.startswith("#") or not s:
                if headers:
                    in_atom = False
                continue
            parts = s.split()
            if len(parts) < len(headers):
                continue
            row = dict(zip(headers, parts))
            atoms.append(row)
    return atoms

def dist(a1, a2):
    return math.sqrt(
        (float(a1["Cartn_x"]) - float(a2["Cartn_x"]))**2 +
        (float(a1["Cartn_y"]) - float(a2["Cartn_y"]))**2 +
        (float(a1["Cartn_z"]) - float(a2["Cartn_z"]))**2
    )

def summarize_structure(cif_path, name):
    print(f"=== {name} ({cif_path}) ===")
    atoms = parse_cif_atoms(cif_path)
    chains = {}
    nonpoly = {}
    for a in atoms:
        ch = a["auth_asym_id"]
        seq = a["auth_seq_id"]
        comp = a["auth_comp_id"]
        group = a["group_PDB"]
        if group == "ATOM" or comp in ("SEP", "TPO", "ALY", "CSO", "PTR"):
            if ch not in chains:
                chains[ch] = {}
            chains[ch][int(seq)] = comp
        else:
            if comp != "HOH":
                nonpoly.setdefault((ch, seq, comp), []).append(a)
    for ch, resmap in sorted(chains.items()):
        seqs = sorted(resmap.keys())
        print(f"Polymer Chain {ch}: {len(seqs)} residues, auth_seq_id {min(seqs)}..{max(seqs)}")
        # print modified residues or key motifs
        mods = [(sq, c) for sq, c in resmap.items() if c not in (
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL")]
        if mods:
            print(f"  Modified residues on Chain {ch}: {mods}")
    for (ch, seq, comp), alist in sorted(nonpoly.items()):
        print(f"Non-polymer: Chain {ch} seq {seq} comp {comp} ({len(alist)} atoms)")
    return atoms, chains, nonpoly

atoms_6eqi, chains_6eqi, nonpoly_6eqi = summarize_structure(
    "/Users/o2satz/.gemini/antigravity/brain/e55fbd0d-7b92-471e-b9ba-63a180b98a1f/6EQI/pdb_00006eqi.cif",
    "6EQI: PINK1 + Ubiquitin"
)

# Check key residues in 6EQI:
# Ubiquitin Ser65 on Ubiquitin chain, and PINK1 Insertion-3, catalytic loop (HRD/DFG), Phospho-Serine, ATP pocket (KTP site)
ub_chain = [ch for ch, rm in chains_6eqi.items() if len(rm) < 100][0]
pink_chain = [ch for ch, rm in chains_6eqi.items() if len(rm) > 200][0]
print(f"PINK1 chain: {pink_chain}, Ubiquitin chain: {ub_chain}")

# Find Ubiquitin Ser65 and nearby PINK1 residues (< 6.0 A)
ub_s65_atoms = [a for a in atoms_6eqi if a["auth_asym_id"] == ub_chain and a["auth_seq_id"] == "65"]
print("Ubiquitin res 65:", set((a["auth_seq_id"], a["auth_comp_id"]) for a in ub_s65_atoms))
pink_near_s65 = {}
for a_ub in ub_s65_atoms:
    for a_pk in atoms_6eqi:
        if a_pk["auth_asym_id"] == pink_chain:
            d = dist(a_ub, a_pk)
            if d < 6.5:
                key = (int(a_pk["auth_seq_id"]), a_pk["auth_comp_id"])
                if key not in pink_near_s65 or d < pink_near_s65[key][0]:
                    pink_near_s65[key] = (d, a_ub["label_atom_id"], a_pk["label_atom_id"])
print("PINK1 residues within 6.5 A of Ub Ser65:")
for k, v in sorted(pink_near_s65.items()):
    print(f"  PINK1 {k[1]}{k[0]} <-> Ub Ser65: {v[0]:.2f} A ({v[1]} - {v[2]})")

# Print sequence around catalytic motifs in PINK1 (Chain A)
res_pink = chains_6eqi[pink_chain]
seq_str = "".join([f"{sq}:{c} " for sq, c in sorted(res_pink.items())])
print("PINK1 residues sample:")
for sq, c in sorted(res_pink.items()):
    if c in ("SEP", "TPO") or sq in (160, 161, 162, 186, 202, 204, 205, 290, 333, 334, 335, 355, 356, 357, 358, 359):
        print(f"  {sq}: {c}")

# Search for HRD / DFG / gatekeeper / hinge in PINK1
sorted_items = sorted(res_pink.items())
for i in range(len(sorted_items) - 2):
    s1, c1 = sorted_items[i]
    s2, c2 = sorted_items[i+1]
    s3, c3 = sorted_items[i+2]
    if (c1, c2, c3) in [("HIS","ARG","ASP"), ("ASP","PHE","GLY"), ("ASP","LEU","LYS")]:
        print(f"  Motif {c1}-{c2}-{c3} at {s1}-{s2}-{s3}")

atoms_4zzj, chains_4zzj, nonpoly_4zzj = summarize_structure(
    "/Users/o2satz/.gemini/antigravity/brain/e55fbd0d-7b92-471e-b9ba-63a180b98a1f/4ZZJ/pdb_00004zzj.cif",
    "4ZZJ: Human SIRT1 + CNA (NAD+) + 4TQ (Activator) + Acetylated Peptide"
)

# Find interactions in 4ZZJ around CNA, 4TQ, and acetylated lysine (ALY) on substrate chain
cna_atoms = [a for a in atoms_4zzj if a["auth_comp_id"] == "CNA" and a["auth_asym_id"] == "A"]
stac_atoms = [a for a in atoms_4zzj if a["auth_comp_id"] == "4TQ" and a["auth_asym_id"] == "A"]
aly_atoms = [a for a in atoms_4zzj if a["auth_comp_id"] == "ALY"]
print("ALY atoms chains/seqs:", set((a["auth_asym_id"], a["auth_seq_id"], a["auth_comp_id"]) for a in aly_atoms))

# Find SIRT1 Chain A residues near CNA, 4TQ, and ALY
near_cna = {}
for a_lig in cna_atoms:
    for a_p in atoms_4zzj:
        if a_p["auth_asym_id"] == "A" and a_p["group_PDB"] == "ATOM":
            d = dist(a_lig, a_p)
            if d < 4.0:
                k = (int(a_p["auth_seq_id"]), a_p["auth_comp_id"])
                if k not in near_cna or d < near_cna[k][0]:
                    near_cna[k] = (d, a_lig["label_atom_id"], a_p["label_atom_id"])
print("SIRT1 Chain A within 4.0 A of CNA (Carba-NAD+):")
for k, v in sorted(near_cna.items()):
    print(f"  SIRT1 {k[1]}{k[0]} <-> CNA: {v[0]:.2f} A ({v[1]} - {v[2]})")

near_stac = {}
for a_lig in stac_atoms:
    for a_p in atoms_4zzj:
        if a_p["group_PDB"] == "ATOM" or a_p["auth_comp_id"] == "ALY":
            d = dist(a_lig, a_p)
            if d < 4.5:
                k = (a_p["auth_asym_id"], int(a_p["auth_seq_id"]), a_p["auth_comp_id"])
                if k not in near_stac or d < near_stac[k][0]:
                    near_stac[k] = (d, a_lig["label_atom_id"], a_p["label_atom_id"])
print("Residues within 4.5 A of 4TQ (Allosteric STAC Activator):")
for k, v in sorted(near_stac.items()):
    print(f"  Chain {k[0]} {k[2]}{k[1]} <-> 4TQ: {v[0]:.2f} A ({v[1]} - {v[2]})")

# Check Glu230 in SIRT1 Chain A
e230_atoms = [a for a in atoms_4zzj if a["auth_asym_id"] == "A" and a["auth_seq_id"] == "230"]
near_e230 = {}
for a_e in e230_atoms:
    for a_p in atoms_4zzj:
        if a_p["auth_asym_id"] == "A" and a_p["auth_seq_id"] != "230" and a_p["group_PDB"] == "ATOM":
            d = dist(a_e, a_p)
            if d < 4.5:
                k = (int(a_p["auth_seq_id"]), a_p["auth_comp_id"])
                if k not in near_e230 or d < near_e230[k][0]:
                    near_e230[k] = (d, a_e["label_atom_id"], a_p["label_atom_id"])
print("SIRT1 Chain A within 4.5 A of Glu230 (Allosteric Hinge):")
for k, v in sorted(near_e230.items()):
    print(f"  SIRT1 {k[1]}{k[0]} <-> Glu230: {v[0]:.2f} A ({v[1]} - {v[2]})")
