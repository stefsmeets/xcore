#!/usr/bin/env python

from spacegroup import generate_hkl_listing
import numpy as np
import pandas as pd

def gaussian(a, b, s):
    """General Gaussian"""
    return a * np.exp(-b * s**2)

def calc_formfactor_xray(data,   s):
    """Scattering factor function for xray/wk1995 table"""
    total = None
    for i in range(5):
        a, b, c = data[0+i], data[5+i], data[10]   # wk95
        y = gaussian(a, b, s)

        if total is None:
            total = y + c
        else:
            total += y
    return total

def calc_formfactor_electron(data,   s):
    """scattering factor function for electron/it432x table"""
    total = None
    for i in range(5):
        a, b, dZ = data[2+i], data[7+i], data[1]
        y = gaussian_fit(a, b, s)

        if total is None:
            total = y
        else:
            total += y
    return total


def calc_formfactor(element, sitl, table="xray"):
    if table == "xray":
        from scattering import wk1995
        table = wk1995.table
        func = calc_formfactor_xray
    elif table == "electron":
        from scattering import it_table_4322
        table = it_table_4322.table
        func = calc_formfactor_electron
    else:
        raise ValueError("Table {} not recognized".format(table))

    data = table[element]
    return func(data, sitl)


def calc_structure_factors(cell, atoms, table="xray"):
    hkl = generate_hkl_listing(cell)
    hkl = hkl[cell.is_absent_np(hkl) == False]
        
    df = pd.DataFrame(hkl, columns=["h", "k", "l"])
    df["sitl"] = 1/(2*cell.calc_dspacing_np(hkl))
    df.sort_values("sitl", ascending=True, inplace=True)
    
    elements = atoms.symbol.unique()

    for element in elements:
        df[element] = calc_formfactor(element, df["sitl"], table=table)

    df["F"] = 0
    
    for i, row in atoms.iterrows():
    #     print row.label, row.symbol, row.x, row.y, row.z, row.occ, row.m
        for i, symop in enumerate(cell.symmetry_operations):
            x, y, z = np.dot(symop.r, [row.x, row.y, row.z]) + symop.t.reshape(3,)
            n = row.occ * row.m / cell.order
            df["F"] += n * df[row.symbol] * np.exp(-row.biso * df.sitl**2) * np.exp(1j*2*np.pi*(x * df.h + y * df.k + z * df.l))
            
    for element in elements:
        del df[element]
    
    df["amplitude"] = (df["F"].real**2 + df["F"].imag**2)**0.5
    df["phase"] = np.arctan2(df["F"].imag, df["F"].real)
    del df["F"]
    del df["sitl"]
    df.phase = df.phase.round(4)
    df.amplitude = df.amplitude.round(4)
    # df.sitl = df.sitl.round(4)
    
    return df


def cif2hkl(fn):
    import cif
    cell, atoms = cif.read_cif(fn)
    df =  calc_structure_factors(cell, atoms)

    print df

    fout = open(fn.replace("cif", "hkl"), "w")

    print " >> Wrote {} reflections to file {}".format(len(df), fout.name)
    print >> fout, df.to_string(sparsify=False, header=False, index=False)


def cif2hkl_entry():
    import sys
    for arg in sys.argv[1:]:
        cif2hkl(arg)
  

if __name__ == '__main__':
    cif2hkl_entry()
