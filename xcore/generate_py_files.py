#!/usr/bin/env python

import os
import subprocess as sp
import time
import re
from spacegroup import lines

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
    drc = "spacegroup"
    last_number = 0
    for row in lines:
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
        for line in out.split('\n'):
            if line.startswith('Point Group'):
                pgr = line.split()[-1]
            if line.startswith('Laue  Group'):
                lauegr = line.split()[-1]
            if line.startswith('Order   '):
                nsym = int(line.split()[-1])
            if line.startswith('Order P '):
                nsymp = int(line.split()[-1])
            if line.startswith('Unique'):
                uniq_axis = line.split()[-1]

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

        reflection_conditions = []
        save_refl = False
        for line in out.split('\n'):
        	if "Reflection conditions" in line:
        		save_refl = True
        		continue
        	if save_refl and not line.split():
        		break
        	if save_refl:
        		reflection_conditions.append(repr(line.strip()))

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
    'reflection_conditions': [
    	{reflection_conditions}
    ],
    'centering_vectors': [
        {centering_vectors}
    ],
    'symops': [
        {symops}
    ]
}}
""".format(number=number,
           setting=repr(setting),
           pgr=repr(pgr),
           lauegr=repr(lauegr),
           nsym=nsym,
           nsymp=nsymp,
           uniq_axis=repr(uniq_axis),
           centrosymmetric=repr(centrosymm),
           reflection_conditions=',\n        '.join(reflection_conditions),
           centering_vectors=',\n        '.join(centering_vecs),
           symops=',\n        '.join(symops)))


if __name__ == '__main__':
    generate_py_files()
