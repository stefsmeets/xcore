#!/usr/bin/env python

import os
import subprocess as sp
import time
import re
from spgr import lines

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


# All reflection conditions
  # h-hl: 2h-l=4n
  # h-hl: 2h-l=6n
  # hhl: 2h+l=4n
  # hk-h: h=2n
  # hk0: k=2n
  # hkl: No Condition
  # hk-k: h-2k=4n
  # hkl: k+l=2n
  # h0l: h=2n
  # h00: h=4n
  # h00: h=2n
  # hkh: 2h+k=4n
  # 00l: l=4n
  # 00l: l=2n
  # hk0: h=2n
  # 00l: l=6n
  # hhl: h=2n
  # 0k0: k=4n
  # 0k0: k=2n
  # hk-k: h=2n
  # h-hl: h=2n
  # h-hl: 2h+l=4n
  # h0l: h+l=2n
  # 0kl: 2k-l=6n
  # hkh: h=2n
  # hkl: h+l=2n
  # hk-h: 2h+k=4n
  # h-hl: l=2n
  # hkk: h+2k=4n
  # hhl: 2h-l=4n
  # 0kl: l=2n
  # h-2hl: l=2n
  # h0l: l=2n
  # hkl: -h+k+l=3n
  # hhl: l=2n
  # hkl: h+k=2n
  # 00l: l=3n
  # hk-h: k=2n
  # hkk: h=2n
  # h0l: 2h+l=6n
  # hk0: h+k=4n
  # hkl: h+k+l=2n
  # -2kkl: l=2n
  # 0kl: k=2n
  # hkh: k=2n
  # 0kl: k+l=2n
  # hk0: h+k=2n
  # h0l: h+l=4n
  # 0kl: k+l=4n