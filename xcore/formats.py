#!/usr/bin/env python

import pandas as pd
import unitcell
import numpy as np
from multiplicity import calc_multiplicity

import shlex

def parse_float(x):
    """Split numerical part from error"""
    if isinstance(x, (str, unicode)):
        return float(x.split("(")[0])
    else:
        return x

def line_iterator(f):
    for line in f:
        line = line.strip()
        if line.startswith("#"):
            continue
        inp = shlex.split(line)
        if not inp:
            continue
        yield inp
    # stopiteration indicates end of file
    # return None to cleanly excit loop
    yield None

def read_cif(f, verbose=True):
    if isinstance(f, str):
        f = open(f, "r")
    d = {}
    incols = False
    inrows = False
    inblock = False
    fiter = line_iterator(f)
    for inp in fiter:
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
                        inp = fiter.next()
                    except StopIteration:
                        incols = False
                        inrows = False
                        inblock = False
                        break
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
                        break
                    if inp[0].startswith("_"):
                        inrows = False
                        inblock = False
                        break
                    elif inp[0] == "loop_":
                        incols == True
                        break
                    assert len(inp) == len(cols), str(cols) + " : " + str(inp)
                    rows.append(inp)
                    
                    try:
                        inp = fiter.next()
                    except StopIteration:
                        incols = False
                        inrows = False
                        inblock = False
                        break
                for i, key in enumerate(cols):
                    vals = [row[i] for row in rows]
                    d[key] = vals


        if not inp:
            continue
        elif inp[0].startswith("data_"):
            d["data_"] = inp[0][5:]
        elif len(inp) == 2:
            key, value = inp
            d[key] = value
        else:
            raise IOError("{} -> Could not read line: {}".format(f.name, inp))
    
    cif2simple = {
    '_atom_site_label': "label",
    '_atom_site_type_symbol': "symbol",
    '_atom_site_symmetry_multiplicity': 'm',
    '_atom_site_fract_x': "x",
    '_atom_site_fract_y': "y",
    '_atom_site_fract_z': "z",
    '_atom_site_occupancy': "occ",
    "_atom_site_adp_type": "adp_type",
    '_atom_site_U_iso_or_equiv': "uiso",
    '_atom_site_B_iso_or_equiv': "biso" }

    a = parse_float( d.pop("_cell_length_a") )
    b = parse_float( d.pop("_cell_length_b") )
    c = parse_float( d.pop("_cell_length_c") )
    al = parse_float( d.pop("_cell_angle_alpha") )
    be = parse_float( d.pop("_cell_angle_beta") )
    ga = parse_float( d.pop("_cell_angle_gamma") )
    sg = d.pop("_symmetry_space_group_name_H-M").replace(" ", "")
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
    
    for col in ["x", "y", "z", "occ", "biso", "uiso"]:
        if col in atoms:
            atoms[col] = atoms[col].apply(parse_float)
    
    if "uiso" in atoms:
        atoms["biso"] = 8*np.pi**2*atoms["uiso"]
    
    if not "m" in atoms:
        atoms["m"] = atoms.apply(calc_multiplicity, args=(cell,), axis=1)
    else:
        atoms["m"] = atoms["m"].astype(int)
        if not all(atoms["m"] == atoms.apply(calc_multiplicity, args=(cell,), axis=1)):
            mults = atoms.apply(calc_multiplicity, args=(cell,), axis=1)
            sel = atoms["m"] != mults
            print "From cif:",
            print atoms[sel]
            print 
            print "Calculated\n"
            print mults[sel]
            raise AssertionError("{} -> multiplicities from cif and calculated are different".format(f.name))
    
    atoms = atoms[["label", "symbol", "m", "x", "y", "z", "occ", "biso"]]
  
    cell.info()

    if verbose:
        print atoms
    else:
        print "{} atoms loaded".format(len(atoms))
    
    return cell, atoms


def write_cif(cell, atoms):
	try:
		cell.write_cif_block()
		atoms.write_cif_block()
	except:
		raise NotImplementedError


def read_hkl(fn):
    df = pd.read_table(fn, sep="\s*", engine="python", index_col=[0,1,2], header=None, names="h k l inty sigma".split())
    df.index = pd.Index(df.index)
    return df


def write_hkl(df, out=None):
    if isinstance(out, str):
        out = open(out, "w")
    if not "sigma" in df:
        df["sigma"] = 1.0
    for row in (h, k, l), row in merged.iterrows():
        print >> out, "{:4d} {:4d} {:4d} {:12.4f} {:12.4f}".format(h, k, l, row["inty"], row["sigma"])


def write_focus_inp(cell, df=None, out="foc.inp", key="amplitude"):
    if isinstance(out, str):
        out = open(out, "w")
    template = """
# Focus input generated by xcore.formats.write_focus_inp

Title {name}
SpaceGroup  {spgr}
UnitCell   {a} {b} {c} {al} {be} {ga}   # {volume}

# AtomType Use Class ScatFact #PerUnitCell OccDefault UisoDefault Radius
# AtomType  +  Node  Si  8
AtomType  +  Node        Si  {nSi}
AtomType  -  NodeBridge  O   {nO}

# Chemistry  MinDistance  Node  Si  Node  Si  2.6
Chemistry  MinDistance  Node        Si  Node        Si  2.6
Chemistry  MinDistance  Node        Si  NodeBridge  O   1.4
Chemistry  MinDistance  NodeBridge  O   NodeBridge  O   2.3

MaxPotentialAtoms  {mpa}
MaxRecycledAtoms   {mra}

FwSearchMethod  FwTracking
MaxPeaksFwSearch  {mpfs}
MaxPeaksFwFragmentSearch  {mpffs}
MinNodeDistance  2.6
MaxNodeDistance  3.7
MinSymNodes  0
MaxSymNodes  {mpa}
NodeType  4  *  -6 -3 -1 4 6
MinLoopSize  4
MaxLoopSize  24
EvenLoopSizesOnly  Off
Check3DimConnectivity  On
IdealT_NodeDistance  3.1
CheckTetrahedralGeometry  Normal

RandomInitialization  Time
FeedBackCycles  1 1 1 1 1 1 1 1 1 1
FeedBackBreakIf  PhaseDiff < 5.00 % and DeltaR < 1.00 %

Grid_xyz  {gridx} {gridy} {gridz}
PeakSearchLevel  3
eDensityCutOff  1 %
MinPfI  17
CatchDistance  0.5
eD_PeaksSortElement  Grid_eD

Lambda  {wl}
FobsMin_d  1
FobsScale  1 # Fabs = Fobs * FobsScale
SigmaCutOff  0
OverlapFactor  0.0
OverlapAction  NoAction
ReflectionUsage  200

# Site  Label  ScatFact  x  y  z  Occ  Uiso
# Site  Si4  Si4+   0.00000   0.00000   0.00000    1.00    0.035

#  h    k    l         Fobs      Sigma     FWHM
{data}
End
"""
    
    nSi = int(20 * cell.volume / 1000.0) / 2 * 2 # make divisible by 2
    nO  = 2*nSi
    
    mpa = int(nSi)
    mra = int(2*nSi)
    mpfs = int(1.5*(nSi + nO))
    mpffs = int(1.5*nSi)
    
    gridx = int(((cell.a*3) // 6 + 1) * 6)
    gridy = int(((cell.b*3) // 6 + 1) * 6)
    gridz = int(((cell.c*3) // 6 + 1) * 6)
    
    data = "# Data goes here\n# ...\n# ...\n# ..."
    if df is not None:
        if not key in df and "inty" in df:
            df["amplitude"] = df["inty"] ** 0.5
        if key in df:
            data = "\n".join(["{:4d} {:4d} {:4d} {:12.4f}".format(h, k, l, row["amplitude"]) for (h, k, l), row in merged.iterrows()])
        else:
            print " >> Data format not understood, provide a df with amplitude/inty"
                             
    print >> out, template.format(
        name=cell.name,
        spgr=cell.spgr_name,
        a=cell.a,
        b=cell.b,
        c=cell.c,
        al=cell.al,
        be=cell.be,
        ga=cell.ga,
        volume=cell.volume,
        gridx=gridx,
        gridy=gridy,
        gridz=gridz,
        nSi=nSi,
        nO=nO,
        mpa=mpa,
        mra=mra,
        mpfs=mpfs,
        mpffs=mpffs,
        wl=1.0,
        data=data
    )
    print "\n >> Wrote focus input file {}".format(out.name)


def make_focus_entry():
    params = raw_input("Enter cell parameters:\n [10 10 10 90 90 90] >> ") or "10 10 10 90 90 90"
    params = [float(val) for val in params.split()]
    spgr = raw_input("Enter spacegroup:\n [P1] >> ") or "P1"

    cell = unitcell.UnitCell(params, spgr)
    
    print
    cell.info()

    fn = raw_input("Enter hkl file (h k l Fobs sigma):\n [optional] >> ") or None

    if fn:
        df = read_hkl(fn)
    else:
        df = None

    write_focus_inp(cell, df=df)


def cif2hkl_entry():
    import sys
    for arg in sys.argv[1:]:
        cell, atoms = read_cif(arg)
        cell.info()
        print atoms


if __name__ == '__main__':
    cif2hkl_entry()
    # make_focus_entry()
