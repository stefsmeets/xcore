#!/usr/bin/env python

import re
import pandas as pd
import unitcell
from multiplicity import calc_multiplicity

_split_cif_line = re.compile(r'\s+(?=[^"\']*?(?:\(|$))')

def parse_float(x):
    """Split numerical part from error"""
    if isinstance(x, (str, unicode)):
        return float(x.split("(")[0])
    else:
        return x

def read_cif(f):
    if isinstance(f, str):
        f = open(f, "r")
    d = {}
    incols = False
    inrows = False
    inblock = False
    for line in f:
        inp = re.split(_split_cif_line, line.strip())
        if not inp:
            continue
        if inp[0] == "loop_":
            inblock = True
            while inblock == True:
                incols = True
                cols = []
                rows = []
                while incols == True:
                    try:
                        line = f.next()
                    except StopIteration:
                        incols = False
                        inrows = False
                        inblock = False
                        break
                    inp = re.split(_split_cif_line, line.strip())
                    if not inp:
                        continue
                    if inp[0].startswith("_"):
                        cols.append(inp[0])
                        continue
                    else:
                        inrows = True
                        incols = False
                while inrows == True:
                    if not inp:
                        continue
                    if inp[0].startswith("_"):
                        inrows = False
                        inblock = False
                        break
                    elif inp[0] == "loop_":
                        incols == True
                        break
                    assert len(inp) == len(cols), str(keys) + str(inp)
                    rows.append(inp)
                    
                    try:
                        line = f.next()
                    except StopIteration:
                        incols = False
                        inrows = False
                        inblock = False
                        break
                    inp = re.split(_split_cif_line, line.strip())
                for i, key in enumerate(cols):
                    vals = [row[i] for row in rows]
                    d[key] = vals
            
        if len(inp) == 2:
            key, value = inp
            d[key] = value
        elif inp[0].startswith("data_"):
            d["data_"] = inp[0][5:]
    
    cif2simple = {
    '_atom_site_label': "label",
    '_atom_site_type_symbol': "symbol",
    '_atom_site_symmetry_multiplicity': 'm',
    '_atom_site_fract_x': "x",
    '_atom_site_fract_y': "y",
    '_atom_site_fract_z': "z",
    '_atom_site_occupancy': "occ",
    '_atom_site_B_iso_or_equiv': "biso" }
    
    a = parse_float( d.pop("_cell_length_a") )
    b = parse_float( d.pop("_cell_length_b") )
    c = parse_float( d.pop("_cell_length_c") )
    al = parse_float( d.pop("_cell_angle_alpha") )
    be = parse_float( d.pop("_cell_angle_beta") )
    ga = parse_float( d.pop("_cell_angle_gamma") )
    sg = d.pop("_symmetry_space_group_name_H-M")
    name = d.pop("data_")
    
    cell = unitcell.UnitCell((a,b,c,al,be,ga), sg, name=name)

    cols = [key for key in d.keys() if key.startswith("_atom")]
    vals = [d[key] for key in cols]
    
    # adps
    # adp_type = d["_atom_site_adp_type"] # "Biso", "Uiso", Uani"
    # Uiso: _atom_site_U_iso_or_equiv
    # Biso: _atom_site_B_iso_or_equiv
    # Uani: 
        # "_atom_site_aniso_label"
        # "_atom_site_aniso_U_11"
        # "_atom_site_aniso_U_22"
        # "_atom_site_aniso_U_33"
        # "_atom_site_aniso_U_23"
        # "_atom_site_aniso_U_13"
        # "_atom_site_aniso_U_12"
    
    atoms = pd.DataFrame([list(row) for row in zip(*vals)], columns=[cif2simple[key] for key in cols])
    
    for col in ["x", "y", "z", "occ", "biso"]:
        atoms[col] = atoms[col].apply(parse_float)
    
    if not "m" in atoms:
        atoms["m"] = calc_multiplicity(atoms, cell)
    else:
        atoms["m"] = atoms["m"].astype(int)
        assert all(atoms["m"] == calc_multiplicity(atoms, cell)), "multiplicities from cif and calculated are different"
    
    atoms = atoms[["label", "symbol", "m", "x", "y", "z", "occ", "biso"]]
  
    cell.info()
    print atoms
    
    return cell, atoms

def write_cif(cell, atoms):
	try:
		cell.write_cif_block()
		atoms.write_cif_block()
	except:
		raise NotImplementedError

if __name__ == '__main__':
	import sys
	for arg in sys.argv[1:]:
		cell, atoms = read_cif(arg)
		cell.info()
		print atoms

