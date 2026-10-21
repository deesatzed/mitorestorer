import json
from analyze_pink1_sirt1 import parse_cif_atoms, dist

atoms_2nsx = parse_cif_atoms("/Users/o2satz/.gemini/antigravity/brain/e55fbd0d-7b92-471e-b9ba-63a180b98a1f/2NSX/pdb_00002nsx.cif")

chains = {}
nonpoly = {}
for a in atoms_2nsx:
    ch = a["auth_asym_id"]
    seq = a["auth_seq_id"]
    comp = a["auth_comp_id"]
    if a["group_PDB"] == "ATOM":
        chains.setdefault(ch, {})[int(seq)] = comp
    elif comp != "HOH":
        nonpoly.setdefault((ch, seq, comp), []).append(a)

for ch, rmap in sorted(chains.items()):
    print(f"Chain {ch}: {len(rmap)} residues ({min(rmap)}..{max(rmap)})")
for k, v in sorted(nonpoly.items()):
    print(f"Ligand: Chain {k[0]} seq {k[1]} comp {k[2]} ({len(v)} atoms)")

# Check IFM (Isofagomine chaperone) on Chain B and Asn370 (rs76763715 N370S) on Chain B
ifm_atoms = [a for a in atoms_2nsx if a["auth_asym_id"] == "B" and a["auth_comp_id"] == "IFM"]
n370_atoms = [a for a in atoms_2nsx if a["auth_asym_id"] == "B" and a["auth_seq_id"] == "370"]

ifm_contacts = {}
for ai in ifm_atoms:
    for ap in atoms_2nsx:
        if ap["auth_asym_id"] == "B" and ap["group_PDB"] == "ATOM":
            d = dist(ai, ap)
            if d < 4.2:
                k = (int(ap["auth_seq_id"]), ap["auth_comp_id"])
                if k not in ifm_contacts or d < ifm_contacts[k][0]:
                    ifm_contacts[k] = (d, ai["label_atom_id"], ap["label_atom_id"])

print("GBA1 Chain B residues within 4.2 A of IFM (Pharmacological Chaperone Pocket):")
for k, v in sorted(ifm_contacts.items()):
    print(f"  GBA1 {k[1]}{k[0]} <-> IFM: {v[0]:.2f} A ({v[1]} - {v[2]})")

# Distance from Asn370 to Trp312 / Loop 1 (311-319) and Catalytic Glu235 / Glu340 / IFM on Chain B vs Chain A
for ch in ("A", "B"):
    n370_ch = [a for a in atoms_2nsx if a["auth_asym_id"] == ch and a["auth_seq_id"] == "370"]
    for target_sq in ("235", "340", "312", "313", "314", "378"):
        t_atoms = [a for a in atoms_2nsx if a["auth_asym_id"] == ch and a["auth_seq_id"] == target_sq]
        if t_atoms and n370_ch:
            min_d = min(dist(a1, a2) for a1 in n370_ch for a2 in t_atoms)
            print(f"[Chain {ch}] Min distance Asn370 <-> res {target_sq} ({t_atoms[0]['auth_comp_id']}): {min_d:.2f} A")

