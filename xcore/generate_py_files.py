#!/usr/bin/env python

import os
import subprocess as sp
import time
import re
from spacegroup import spacegrouptxt

from wyckoff import general_position_multiplicities, wyckofftable

get_centering = re.compile("\((.*)\)")


def get_transvec(string):
    transvec = []
    for i, val in enumerate(string.strip('()').split()):
        if '/' in val:
            frac = val.split('/')
            num = float(frac[0])
            denom = float(frac[1])
            transvec.append(num/denom)
        else:
            transvec.append(float(val))
    return transvec


def generate_py_files():
    exe = 'sginfo'
    drc = "spacegroup_tables"
    last_number = 0
    for i, row in enumerate(spacegrouptxt):
        # if i == 10:
            # exit()
        number = row[:12].strip()

        try:
            number, setting = number.split(':')
        except ValueError:
            setting = ""

        if last_number != number:
            fn = open(os.path.join(drc, "{}.py".format(number)), 'w')
            print "new file -", fn.name
            fn.write("# Generated {}".format(time.ctime()))
            fn.write("\nd = {}")
        print "{}:{}".format(number, setting)
        last_number = number

        spgr = '{}:{}'.format(number, setting)
        cmd = [exe, spgr, '-allxyz']
        p = sp.Popen(cmd, stdout=sp.PIPE)
        out, err = p.communicate()

        symops = []
        centering_vecs = []
        centrosymm = False
        save_symops = False
        uniq_axis = None
        chiral_spgr = False
        enantiomorphic = False
        off_origin = False
        obverse = False
        for line in out.split('\n'):
            line = line.strip()
            if line.startswith('Point Group'):
                pgr = line.split()[-1]
            elif line.startswith('Laue  Group'):
                lauegr = line.split()[-1]
            elif line.startswith('Order   '):
                nsym = int(line.split()[-1])
            elif line.startswith('Order P '):
                nsymp = int(line.split()[-1])
            elif line.startswith('Unique'):
                uniq_axis = line.split()[-1]
            elif line.startswith("Chiral space group"):
                chiral_spgr = True
            elif line.startswith("Enantiomorphic"):
                enantiomorphic = True
            elif line.startswith("Obverse"):
                obverse = True
            elif line.startswith("Note: Inversion operation off origin"):
                off_origin = True

            m = re.search(get_centering, line)
            if m:
                centering = get_transvec(m.group())
                if centering not in centering_vecs:
                    centering_vecs.append(centering)
            if "Inversion-Flag = 1" in line:
                centrosymm = True

            if '(0 0 0)' in line or line.startswith('x, y, z'):
                save_symops = True
            if save_symops and line.startswith('#'):
                save_symops = False
            if save_symops and not line:
                save_symops = False

            if save_symops:
                symops.append(line)

        spgr = '{}:{}'.format(number, setting)
        cmd = [exe, spgr, '-Conditions']
        p = sp.Popen(cmd, stdout=sp.PIPE)
        out, err = p.communicate()

        print "number", number
        mult = general_position_multiplicities[int(number)]
        wyckoffraw = wyckofftable.get(int(number), ())

        wyckoff_positions = reversed([repr((mult, "x, y, z"))] + [repr(w) for w in wyckoffraw])

        reflection_conditions = []
        phase_restrictions = []
        enhanced_reflections = []
        save_refl_cond = False
        save_refl_phase = False
        save_refl_enh = False
        for line in out.split('\n'):
            if line.startswith("Reflection conditions"):
                save_refl_cond = True
                continue
            if line.startswith("Reflections with phase restriction"):
                save_refl_phase = True
                continue
            if line.startswith("Systematically enhanced reflections"):
                save_refl_enh = True
                continue
            
            if not line.split():
                save_refl_cond = False
                save_refl_phase = False
                save_refl_enh = False
            
            if save_refl_phase:
                phase_restrictions.append(repr(line.strip()))
            if save_refl_cond:
                reflection_conditions.append(repr(line.strip()))
            if save_refl_enh:
                if len(line.split("=")[-1]) <=4:
                    enhanced_reflections.append(repr(line.strip()))
                else:
                    print "skip", line

        if len(phase_restrictions) == 0:
            phase_restrictions = ["'hkl: No Condition'"]

        if not centering_vecs:
            centering_vecs = [[0.0, 0.0, 0.0]]
        centering_vecs = [repr(vec) for vec in centering_vecs]

        symops = [repr(symop) for symop in symops]
        fn.write("""
d[{setting}] = {{
    'number': {number},
    'setting': {setting},
    'point_group': {pgr},
    'laue_group': {lauegr},
    'order': {nsym},
    'order_p': {nsymp},
    'unique_axis': {uniq_axis},
    'centrosymmetric': {centrosymmetric},
    'enantiomorphic': {enantiomorphic},
    'chiral': {chiral},
    'obverse': {obverse},
    'reflection_conditions': (
        {reflection_conditions},
    ),
    'enhanced_reflections': (
        {enhanced_reflections},
    ),
    'phase_restrictions': (
        {phase_restrictions},
    ),
    'centering_vectors': (
        {centering_vectors},
    ),
    'symops': (
        {symops},
    ),
    'wyckoff_positions' : (
        {wyckoff_positions},
    )
}}
""".format(number=number,
           setting=repr(setting),
           pgr=repr(pgr),
           lauegr=repr(lauegr),
           nsym=nsym,
           nsymp=nsymp,
           uniq_axis=repr(uniq_axis),
           centrosymmetric=repr(centrosymm),
           enantiomorphic=repr(enantiomorphic),
           chiral=repr(chiral_spgr),
           obverse=repr(obverse),
           reflection_conditions=',\n        '.join(reflection_conditions),
           enhanced_reflections=',\n        '.join(enhanced_reflections),
           phase_restrictions=',\n        '.join(phase_restrictions),
           centering_vectors=',\n        '.join(centering_vecs),
           symops=',\n        '.join(symops),
           wyckoff_positions=',\n        '.join(wyckoff_positions)))


if __name__ == '__main__':
    generate_py_files()
