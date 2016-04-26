#!/usr/bin/env python

import numpy as np
import pandas as pd

def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

def calc_multiplicity(atoms, spgr, precision=5):
    # checking if it is a dataframe, so we have more than 1 atom
    # a single atom is supplied as a series object
    if isinstance(atoms, pd.DataFrame):
        return atoms.apply(calc_multiplicity, args=(spgr,precision), axis=1) # recursion like a boss
    arr = []
    for i, symop in enumerate(spgr.symmetry_operations):
        x, y, z = np.dot(symop.r, [atoms.x, atoms.y, atoms.z]) + symop.t.reshape(3,)
        arr.append((x, y, z))
    arr = np.array(arr).round(precision) % 1
    return len(unique_rows(arr))





