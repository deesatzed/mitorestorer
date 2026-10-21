from analyze_pink1_sirt1 import parse_cif_atoms, dist

atoms_6eqi = parse_cif_atoms("/Users/o2satz/.gemini/antigravity/brain/e55fbd0d-7b92-471e-b9ba-63a180b98a1f/6EQI/pdb_00006eqi.cif")

# Check SEP202, SEP204, TPO305 contacts in PINK1 Chain C and Ubiquitin Chain A
for mod_seq in ("202", "204", "305"):
    mod_atoms = [a for a in atoms_6eqi if a["auth_asym_id"] == "C" and a["auth_seq_id"] == mod_seq]
    contacts = {}
    for am in mod_atoms:
        for ap in atoms_6eqi:
            if ap["auth_asym_id"] in ("A", "C") and ap["auth_seq_id"] != mod_seq:
                d = dist(am, ap)
                if d < 4.2:
                    k = (ap["auth_asym_id"], int(ap["auth_seq_id"]), ap["auth_comp_id"])
                    if k not in contacts or d < contacts[k][0]:
                        contacts[k] = (d, am["label_atom_id"], ap["label_atom_id"])
    print(f"Contacts (<4.2 A) of PINK1 Chain C res {mod_seq}:")
    for k, v in sorted(contacts.items()):
        print(f"  Chain {k[0]} {k[2]}{k[1]}: {v[0]:.2f} A ({v[1]} - {v[2]})")

# Check Ubiquitin Chain A contacts with PINK1 Chain C
ub_contacts = {}
for au in atoms_6eqi:
    if au["auth_asym_id"] == "A":
        for ap in atoms_6eqi:
            if ap["auth_asym_id"] == "C":
                d = dist(au, ap)
                if d < 3.8:
                    k = (int(au["auth_seq_id"]), au["auth_comp_id"], int(ap["auth_seq_id"]), ap["auth_comp_id"])
                    if k not in ub_contacts or d < ub_contacts[k][0]:
                        ub_contacts[k] = (d, au["label_atom_id"], ap["label_atom_id"])
print("Direct H-bonds/salt bridges (<3.8 A) between Ubiquitin (A) and PINK1 (C):")
for k, v in sorted(ub_contacts.items()):
    print(f"  Ub {k[1]}{k[0]} <-> PINK1 {k[3]}{k[2]}: {v[0]:.2f} A ({v[1]} - {v[2]})")

# Now check 4ZZJ ALY4 contacts and Glu230 long-range interaction with Catalytic Core (e.g. Arg446)
atoms_4zzj = parse_cif_atoms("/Users/o2satz/.gemini/antigravity/brain/e55fbd0d-7b92-471e-b9ba-63a180b98a1f/4ZZJ/pdb_00004zzj.cif")
aly_atoms = [a for a in atoms_4zzj if a["auth_asym_id"] == "B" and a["auth_seq_id"] == "4"]
aly_contacts = {}
for aa in aly_atoms:
    for ap in atoms_4zzj:
        if ap["auth_asym_id"] == "A":
            d = dist(aa, ap)
            if d < 4.2:
                k = (int(ap["auth_seq_id"]), ap["auth_comp_id"])
                if k not in aly_contacts or d < aly_contacts[k][0]:
                    aly_contacts[k] = (d, aa["label_atom_id"], ap["label_atom_id"])
print("4ZZJ Acetyl-Lysine (Chain B ALY4) contacts (<4.2 A) with SIRT1 Chain A:")
for k, v in sorted(aly_contacts.items()):
    print(f"  ALY4 <-> SIRT1 {k[1]}{k[0]}: {v[0]:.2f} A ({v[1]} - {v[2]})")

# Check distance between 4TQ (703) and Substrate Peptide (Chain B)
stac_atoms = [a for a in atoms_4zzj if a["auth_comp_id"] == "4TQ" and a["auth_asym_id"] == "A"]
for ab in [a for a in atoms_4zzj if a["auth_asym_id"] == "B"]:
    for ast in stac_atoms:
        d = dist(ab, ast)
        if d < 5.0:
            print(f"  4TQ <-> Substrate Chain B {ab['auth_comp_id']}{ab['auth_seq_id']} ({ab['label_atom_id']}): {d:.2f} A")
