#!/usr/bin/env python

import os
import sys
import re

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

f = open(os.path.join(os.path.dirname(__file__), 'spacegroups.txt'), 'r')
lines = f.readlines()
f.close()


def getsymm(text):
    """Parses symmetry strings and returns a rotation matrix and translation vector."""
    rotmat = np.zeros([3, 3], dtype=int)
    transvec = np.zeros([3, 1], dtype=float)
    split = text.split(',')
    i = 0       # determines the row in the matrix
    for symeq in split:
        # separate x,y,z,+,-,float,int,fraction
        q = re.findall('[\+\-xXyYzZ]|[0-9]*[\.\/][0-9]+|[0-9]', symeq)
        coefficient = 1
        for r in q:
            if r == '-':
                coefficient = -1
            elif r == '+':
                coefficient = 1
            elif r.lower() == 'x':
                rotmat[i, 0] += coefficient
            elif r.lower() == 'y':
                rotmat[i, 1] += coefficient
            elif r.lower() == 'z':
                rotmat[i, 2] += coefficient
            elif '/' in r:          # check for fraction
                frac = re.findall('[0-9]+', r)
                num = float(frac[0])
                denom = float(frac[1])
                transvec[i, 0] = num/denom
            else:
                try:                # check for float or integer
                    r = float(r)
                    transvec[i, 0] = coefficient*r
                except:
                    pass
        i += 1
    return rotmat, transvec


def find_number(s, lines=lines):
    if not isinstance(s, str):
        s = str(s)
    for i, line in enumerate(lines):
        m = re.search(" "+s+"[: ]+", line)
        if m:
            return line
    raise ValueError("Cannot find space group {}".format(s))


def get_standard_setting(number, setting=None):
    line = find_number(number)
    number = line[:12].strip()
    try:
        number, setting = number.split(':')
    except ValueError:
        setting = ""
    return number, setting


def get_symmetry(number, setting=None):
    if not setting:
        number, setting = get_standard_setting(number)

    import importlib
    drc = 'spacegroup'
    fn = '{}.{}'.format(drc, number)
    spgr = importlib.import_module(fn, package='.')
    return spgr.d[setting]


def get_spacegroup_info(string):
    line = find_number(string)
    if not line:
        print "Could not find {}".format(string)
        return None
    number = line[:12].strip()
    schoenflies = line[12:26].strip()
    hm = line[26:54].strip()
    hall = line[54:].strip()

    try:
        number, setting = number.split(':')
    except ValueError:
        setting = ""
    finally:
        number = int(number)

    if number <= 2:
        crystal_system = "Triclinic"
    elif number <= 15:
        crystal_system = "Monoclinic"
    elif number <= 74:
        crystal_system = "Orthorhombic"
    elif number <= 142:
        crystal_system = "Tetragonal"
    elif number <= 167:
        crystal_system = "Trigonal"
    elif number <= 194:
        crystal_system = "Hexagonal"
    elif number <= 230:
        crystal_system = "Cubic"
    else:
        raise ValueError(
            "This should not happen, does space group {} not exist?".format(number))

    spgr = get_symmetry(number)
    spgr["crystal_system"] = crystal_system

    if setting:
        spgr["spgr"] = "{}:{}".format(number, setting)
    else:
        spgr["spgr"] = "{}".format(number)

    spgr["setting"] = setting
    spgr["number"] = number
    spgr["hall"] = hall
    spgr["hm"] = hm
    spgr["schoenflies"] = schoenflies

    return spgr


def main():
    for arg in sys.argv[1:]:
        spgr = get_spacegroup_info(arg)

        if not spgr:
            continue
        print "# {}\n".format(arg)

        print "Space group ", spgr["spgr"]
        print "             {hm},  {schoenflies},  {hall}".format(**spgr)
        print "Laue  group ", spgr['laue_group']
        print "Point group ", spgr['point_group']

        print spgr["crystal_system"]
        if spgr["centrosymmetric"]:
            print "Centrosymmetric"

        if spgr['unique_axis']:
            print "Unique axis", spgr['unique_axis']
        print
        print "Order    ", spgr["order"]
        print "Order P  ", spgr["order_p"]
        print
        cvecs = spgr["centering_vectors"]
        if len(cvecs) > 1:
            print "Centering vectors"
            for vec in cvecs:
                print "+({}  {}  {})".format(*vec)
            print

        print "Symmetry operations"
        for symop in spgr["symops"]:
            print symop


if __name__ == '__main__':
    main()

    # for x in (seq 1 230)
    #       echo "$x",
    #       sginfo $x | grep -E "Triclinic|Monoclinic|Orthorhombic|Tetragonal|Trigonal|Hexagonal|Cubic"
    #   end
    #
    # for x in (seq 1 230)
    #       echo "$x",
    #       sginfo $x | grep -E "Space Group"
    #       sginfo $x | grep -E "Point Group"
    #       sginfo $x | grep -E "Laue  Group"
    #   end
