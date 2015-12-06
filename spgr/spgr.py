#!/usr/bin/env python

import os
import sys
import re

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

f = open(os.path.join(os.path.dirname(__file__),'spacegroups.txt'),'r')
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
	for i,line in enumerate(lines):
		m = re.search(" "+s+"[: ]+", line)
		if m:
			return line
	raise ValueError("Cannot find space group {}".format(i))


def get_standard_setting(number, setting=None):
	line = find_number(number)
	number = line[:12].strip()
	try:
		number, setting = number.split(':')
	except ValueError:
		setting = ""
	return number, setting


def generate_py_files():
	import subprocess as sp
	import time
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
			first = True
			fn = open(os.path.join(drc, "sg{}.py".format(number)), 'w')
			print "new file -", fn.name
			fn.write("# Generated {}".format(time.ctime()))
			fn.write("\nd = {}")
		print "{}:{}".format(number, setting)
		last_number = number

		spgr = '{}:{}'.format(number, setting)
		cmd = [exe, spgr, '-allxyz']	
		p = sp.Popen(cmd, stdout=sp.PIPE)
		out,err = p.communicate()

		symops = []
		save_symops = False
		uniq_axis = None
		for line in out.split('\n'):
			if line.startswith('Point Group'):
				pgr = line.split()[-1]
			if line.startswith('Laue  Group'):
				lauegr = line.split()[-1]
			if line.startswith('Order   '):
				nsym  = int(line.split()[-1])
			if line.startswith('Order P '):
				nsymp = int(line.split()[-1])
			if line.startswith('Unique'):
				uniq_axis = line.split()[-1]

			if '(0 0 0)' in line or line.startswith('x, y, z'):
				save_symops = True
			if not line.split():
				save_symops = False

			if save_symops:
				symops.append(line)

		symops = [repr(symop) if not symop.startswith('#') else symop.replace(',','') for symop in symops]
		fn.write( """
d[{setting}] = {{
	'number': {number},
	'setting': {setting},
	'point_group': {pgr},
	'laue_group': {lauegr},
	'order': {nsym},
	'order_p': {nsymp},
	'unique_axis': {uniq_axis},
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
					symops=',\n'.join(symops)) )


def get_symmetry(number, setting=None):
	if not setting:
		number, setting = get_standard_setting(number)

	import importlib
	drc = 'spacegroup'
	print __name__
	fn = '{}.sg{}'.format(drc,number)
	spgr = importlib.import_module(fn, package='.')
	return spgr.d[setting]


def main():
	for arg in sys.argv[1:]:
		line = find_number(arg)
		if not line:
			print "Could not find {}".format(arg)
			continue
		# inp = line.split('|')
		# inp = [s.strip() for s in inp]
		# number, schoenflies, hm, hall = inp
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
		# assert number == arg, "{} != {}".format(number, arg)
		if setting:
			print "Space Group {}:{}".format(number, setting)
		else:
			print "Space Group {}".format(number)
		print "Schoenflies", schoenflies
		print "Hermann-Mauguin", hm
		print "Hall", hall

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
			raise ValueError("This should not happen, does space group {} not exist?".format(number))

		print "Crystal system:", crystal_system

		spgr = get_symmetry(number)

		print "Laue  group", spgr['laue_group']
		print "Point group", spgr['point_group']
		if spgr['unique_axis']:
			print "Unique axis", spgr['unique_axis']
		print
		print "Order", spgr["order"]
		print "Order P", spgr["order_p"]
		print
		for symop in spgr["symops"]:
			print symop
			# print getsymm(symop)
		# print spgr





if __name__ == '__main__':
	# generate_py_files()
	main()


	# for x in (seq 1 230)
 #       echo "$x",
 #       sginfo $x | grep -E "Triclinic|Monoclinic|Orthorhombic|Tetragonal|Trigonal|Hexagonal|Cubic"
 #   end

 #   	for x in (seq 1 230)
 #       echo "$x",
 #       sginfo $x | grep -E "Space Group"
 #       sginfo $x | grep -E "Point Group"
 #       sginfo $x | grep -E "Laue  Group"
 #   end