#!/usr/bin/env python

import os
import sys
import re

import numpy as np
import pandas as pd


sys.path.insert(0, os.path.dirname(__file__))

f = open(os.path.join(os.path.dirname(__file__), 'spacegroups.txt'), 'r')
spacegrouptxt = f.readlines()
f.close()


from IPython.terminal.embed import InteractiveShellEmbed
InteractiveShellEmbed.confirm_exit = False
ipshell = InteractiveShellEmbed(banner1='')

reflection_conditions = {
    "00l": lambda h, k, l: abs(h)+abs(k) == 0,
    "0k0": lambda h, k, l: abs(h)+abs(l) == 0,
    "h00": lambda h, k, l: abs(k)+abs(l) == 0,
    "0kl": lambda h, k, l: h == 0,
    "h0l": lambda h, k, l: k == 0,
    "hk0": lambda h, k, l: l == 0,
    "hkl": lambda h, k, l: np.ones_like(h),
    "hhl": lambda h, k, l: (h-k) == 0,
    "hkh": lambda h, k, l: (h-l) == 0,
    "hkk": lambda h, k, l: (k-l) == 0,
    "h-hl": lambda h, k, l: (h+k) == 0,
    "hk-h": lambda h, k, l: (h+l) == 0,
    "hk-k": lambda h, k, l: (k+l) == 0,
    "-2kkl": lambda h, k, l: (h == -2*k),
    "h-2hl": lambda h, k, l: (k == -2*h),
    "l=6n": lambda h, k, l: (l % 6) == 0,
    "l=4n": lambda h, k, l: (l % 4) == 0,
    "l=3n": lambda h, k, l: (l % 3) == 0,
    "l=2n": lambda h, k, l: (l % 2) == 0,
    "k=4n": lambda h, k, l: (k % 4) == 0,
    "k=2n": lambda h, k, l: (k % 2) == 0,
    "k+l=4n": lambda h, k, l: ((k+l) % 4) == 0,
    "k+l=2n": lambda h, k, l: ((k+l) % 2) == 0,
    "h=4n": lambda h, k, l: (h % 4) == 0,
    "h=2n": lambda h, k, l: (h % 2) == 0,
    "h-2k=4n": lambda h, k, l: ((h-2*k) % 4) == 0,
    "h+l=4n": lambda h, k, l: ((h+l) % 4) == 0,
    "h+l=2n": lambda h, k, l: ((h+l) % 2) == 0,
    "h+k=4n": lambda h, k, l: ((h+k) % 4) == 0,
    "h+k=2n": lambda h, k, l: ((h+k) % 2) == 0,
    "h+k+l=2n": lambda h, k, l: ((h+k+l) % 2) == 0,
    "h+2k=4n": lambda h, k, l: ((h+2*k) % 4) == 0,
    "2k-l=6n": lambda h, k, l: ((2*k-l) % 6) == 0,
    "2h-l=6n": lambda h, k, l: ((2*h-l) % 6) == 0,
    "2h-l=4n": lambda h, k, l: ((2*h-l) % 4) == 0,
    "2h+l=6n": lambda h, k, l: ((2*h+l) % 6) == 0,
    "2h+l=4n": lambda h, k, l: ((2*h+l) % 4) == 0,
    "2h+k=4n": lambda h, k, l: ((2*h+k) % 4) == 0,
    "-h+k+l=3n": lambda h, k, l: ((-h+k+l) % 3) == 0,
    "No Condition": lambda h, k, l: True
}

# laue_merge = {
#     'triclinic' : lambda (h,k,l): (abs(h),k,l),
#     'monoclinic' : lambda (h,k,l): (abs(h),abs(k),l),
#     'monoclinicb' : lambda (h,k,l): (abs(h),abs(k),l),
#     'monoclinicc' : lambda (h,k,l): (abs(h), l,abs(l)),
#     'orthorhombic' : lambda (h,k,l): (abs(h),abs(k),abs(l))
#     #'tetragonal' : lambda h,k,l: raise ValueError,
#     #'hexagonal' : lambda h,k,l: raise ValueError,
#     #'rhombohedral1' : lambda h,k,l: raise ValueError,
#     #'rhombohedral2' : lambda h,k,l: raise ValueError,
#     #'cubic' : lambda h,k,l: raise ValueError
# }

def asym_1((h,k,l)):
    ret = False
    if (h >= 0):
        ret = True

    if (k < 0) and (h == 0):
        ret = False
    if (l < 0) and (h == k == 0):
        ret = False
    return ret


def asym_2a((h,k,l)):
    ret = False
    if (h >= 0) and (l >= 0):
        ret = True
    
    if (l < 0) and (k == 0):
        ret = False
    return ret

def asym_2b((h,k,l)):
    ret = False
    if (k >= 0) and (l >= 0):
        ret = True

    if (h < 0) and (l == 0):
        ret = False
    return ret

def asym_2c((h,k,l)):
    ret = False
    if (k >= 0) and (l >= 0):
        ret = True

    if (h < 0) and (k == 0):
        ret = False
    return ret

def asym_3((h,k,l)):
    ret = False
    if (h >= 0) and (k >= 0) and (l >= 0):
        ret = True
    return ret

def asym_4a((h,k,l)):
    ret = False
    if (h >= 0) and (k >= 0) and (l >= 0):
        ret = True

    if (k > 0) and (h == 0):
        ret = False
    return ret

def asym_4b((h,k,l)):
    ret = False
    if (k >= h >= 0) and (l >= 0):
        ret = True
    return ret

def asym_5a((h,k,l)):
    ret = False
    if (k >= 0) and (l >= 0) and (h <= k) and (h <= l):
        ret = True

    if (h < k) and (h == l):
        ret = False
    if (l < 0) and (h ==0):
        ret = False
    return ret

def asym_5b((h,k,l)):
    ret = False
    if (h >= 0) and (k >= 0):
        ret = True

    if (l < 0) and (k ==0):
        ret = False
    if (l <= 0) and (h == 0):
        ret = False
    return ret

def asym_5c((h,k,l)):
    ret = False
    if (l >= k >= 0) and (k >= h):
        ret = True

    if (h < 0) and (abs(h) > l) and (l == 0):
        ret = True
    return ret

def asym_5d((h,k,l)):
    ret = False
    if (k >= h >= 0):
        ret = True

    if (l < 0) and (h == 0):
        ret = False
    return ret

def asym_5e((h,k,l)):
    ret = False
    if (h >= 0) and (k >= 0) and (l >= 0):
        ret = True
    
    if (h < k) and (l == 0):
        ret = False
    return ret

def asym_6a((h,k,l)):
    ret = False
    if (h >= 0) and (k >= 0) and (l >= 0):
        ret = True

    if (k > 0) and (h == 0):
        ret = False
    return ret

def asym_6b((h,k,l)):
    ret = False
    if (k >= h >= 0) and (l >= 0):
        ret = True
    return ret

def asym_7a((h,k,l)):
    ret = False
    if (k >= h >= 0) and (l >= h >= 0):
        ret = True

    if (h < k) and (h == l):
        ret = False
    return ret

def asym_7b((h,k,l)):
    ret = False
    if (l >= k >= h >= 0):
        ret = True
    return ret

# http://scripts.iucr.org/cgi-bin/paper?se0073
laue_asymmetry = {
    # Triclinic
    '-1': asym_1,  # 1/2
    # Monoclinic
    '2/m:a': asym_2a,  # 1/4
    '2/m:b': asym_2b,  # 1/4
    '2/m:c': asym_2c,  # 1/4
    # Orthorhombic
    'mmm': asym_3,  # 1/8
    # Tetragonal
    '4/m': asym_4a,  # 1/8
    '4/mmm': asym_4b,  # 1/16
    # Trigonal, Rhombohedral setting
    '-3:R': asym_5a,  # 1/6
    '-3m:R': asym_5b,  # 1/12
    # Trigonal, Hexagonal setting
    '-3:H': asym_5c,  # 1/6
    '-31m:H': asym_5d,  # 1/12
    '-3m1:H': asym_5e,  # 1/12
    # Hexagonal
    '6/m': asym_6a,  # 1/12
    '6/mmm': asym_6b,  # 1/24
    # Cubic
    'm-3': asym_7a,  # 1/24
    'm-3m': asym_7b  # 1/48
}

laue_fraction = {
    '-1': "1/2",
    '2/m:a': "1/4",
    '2/m:b': "1/4",
    '2/m:c': "1/4",
    'mmm': "1/8",
    '4/m': "1/8",
    '4/mmm': "1/16",
    '-3:R': "1/6",
    '-3m:R': "1/12",
    '-3:H': "1/6",
    '-31m:H': "1/12",
    '-3m1:H': "1/12",
    '6/m': "1/12",
    '6/mmm': "1/24",
    'm-3': "1/24",
    'm-3m': "1/48"
}


def find_number(s, spacegrouptxt=spacegrouptxt):
    if isinstance(s, int):
        s = str(s)
    assert isinstance(
        s, str), "Internal Error: variable s is of wrong type {}".format(type(s))
    for i, line in enumerate(spacegrouptxt):
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


def get_symmetry(number, setting):
    if not setting:
        number, setting = get_standard_setting(number)

    import importlib
    drc = 'spacegroup'
    fn = '{}.{}'.format(drc, number)
    spgr = importlib.import_module(fn, package='.')
    return spgr.d[setting]


def get_random_cell(spgr):
    import random
    a = float(random.randrange(500, 2000)) / 100
    b = float(random.randrange(500, 2000)) / 100
    c = float(random.randrange(500, 2000)) / 100
    al = float(random.randrange(60, 120))
    be = float(random.randrange(60, 120))
    ga = float(random.randrange(60, 120))

    system = spgr.crystal_system
    lauegr = spgr.laue_group
    setting = spgr.setting

    if system == "Triclinic":
        return (a, b, c, al, be, ga)
    elif system == "Monoclinic":
        return (a, b, c, be)
    elif system == "Orthorhombic":
        return (a, b, c)
    elif system == "Tetragonal":
        return (a, c)
    elif system == "Trigonal":
        if setting == "R":
            return (a, al)
        else:
            return (a, c)
    elif system == "Hexagonal":
        return (a, c)
    elif system == "Cubic":
        return (a,)
    else:
        raise ValueError("Invalid system {}".format(system))



def get_spacegroup_info(string, as_dict=False):
    if string == "random":
        import random
        line = random.choice(spacegrouptxt)
    else:
        try:
            line = find_number(string)
        except ValueError:
            raise ValueError("Could not find space group {}".format(string))
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

    spgr = get_symmetry(number, setting)

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

    if as_dict:
        return spgr
    else:
        return SpaceGroup(spgr)


class CVec(tuple):

    """Small class for storing and representing Centering Vectors"""

    def __init__(self, items):
        super(CVec, self).__init__()

    def __repr__(self):
        return "+({}  {}  {})".format(*self)


def symm2str(r,t):
    xyz = "xyz"
    string_all = []
    for i in range(3):
        string_row = ""
        for j in range(3):
            rval = r[i,j]
            if rval == -1:
                string_row += "-{}".format(xyz[j])
            if rval == 1:
                if string_row != "":
                    string_row += "+"
                string_row += "{}".format(xyz[j])
        tval = t[i,0]
        if tval%1 == 0.5:
            tval = "+1/2"
        elif tval%1 == 0.25:
            tval = "+1/4"
        elif tval%1 == 0.75:
            tval = "+3/4"
        elif tval%1 == 1.0/3.0:
            tval = "+1/3"
        elif tval%1 == 2.0/3.0:
            tval = "+2/3"
        elif tval%1 == 1.0/6.0:
            tval = "+1/6"
        elif tval%1 == 5.0/6.0:
            tval = "+5/6"
        else:
            tval = ""
        string_row += "{}".format(tval)
        string_all.append(string_row)
    return ", ".join(string_all)


class SymOp(object):

    """Generate symmetry operations from string"""

    def __init__(self, s):
        super(SymOp, self).__init__()
        self._s = s
        self.r, self.t = self.getsymm(self._s)

        string = symm2str(self.r, self.t)
        assert string == s, "got {}, need {}".format(string, s)

    def __repr__(self):
        return "{}".format(self._s)

    def __get__(self):
        return self.r, self.t

    def getsymm(self, op=None):
        """Parses symmetry strings and returns a rotation matrix and translation vector."""
        rotmat = np.zeros([3, 3], dtype=int)
        transvec = np.zeros([3, 1], dtype=float)
        if not op:
            op = self._s
        split = op.split(',')
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

    def inverse(self):
        r = np.dot(self.r, -np.eye(3,3))
        t = self.t
        return SymOp(symm2str(r,t))



class ReflCond(object):

    """Generate reflection conditions from string
    Uses lookup table to find corresponding operations"""

    def __init__(self, s):
        super(ReflCond, self).__init__()
        self._s = s
        self._ref, self._cond = self.build_reflection_condition(self._s)

    def __repr__(self):
        return self._s

    def build_reflection_condition(self, string=None):
        if not string:
            string = self._s

        ref, cond = string.split(":")
        ref = reflection_conditions[ref.strip()]
        cond = reflection_conditions[cond.strip()]

        return ref, cond

    def test(self, index):
        # return self._ref(*index) and self._cond(*index)
        return self._ref(*index) & self._cond(*index)

    def is_absent(self, index):
        # return self._ref(*index) and not self._cond(*index)
        return self._ref(*index) & ~self._cond(*index)


class SpaceGroup(object):

    """SpaceGroup class that takes a spgr dict (get_spacegroup_info)
    Stores all space group info that could be extracted from sginfo"""

    def __init__(self, kwargs):
        super(SpaceGroup, self).__init__()

        self.hall = kwargs["hall"]
        self.hermann_mauguin = kwargs["hm"]
        self.schoenflies = kwargs["schoenflies"]

        self.point_group = kwargs["point_group"]
        self.laue_group = kwargs["laue_group"]

        self.crystal_system = kwargs["crystal_system"]

        self.number = kwargs["number"]
        self.setting = kwargs["setting"]

        self.name = self.hermann_mauguin.split(" = ")[0]

        self.centering_vectors = [CVec(cv)
                                  for cv in kwargs["centering_vectors"]]
        self.symmetry_operations = [SymOp(op) for op in kwargs["symops"]]

        self.reflection_conditions = [
            ReflCond(rc) for rc in kwargs["reflection_conditions"]]

        self.unique_axis = kwargs["unique_axis"]

        assert self.order_p == kwargs["order_p"], "{} {} {}".format(
            self.space_group, self.order_p, kwargs["order_p"])
        assert self.order == kwargs["order"], "{} {} {}".format(
            self.space_group, self.order, kwargs["order"])
        assert self.is_centrosymmetric == kwargs["centrosymmetric"]

    def __repr__(self):
        return self.name

    @property
    def space_group(self):
        if self.setting:
            return "{}:{}".format(self.number, self.setting)
        else:
            return "{}".format(self.number)

    @property
    def centering(self):
        if self.hall[0] == "-":
            return self.hall[1]
        else:
            return self.hall[0]

    @property
    def is_centrosymmetric(self):
        return self.hall[0] == "-"

    @property
    def order_p(self):
        n = 2 if self.is_centrosymmetric else 1
        return len(self.symmetry_operations) * n

    @property
    def order(self):
        return self.order_p * len(self.centering_vectors)

    def info(self):
        print "Space group ", self
        print
        print "    Number      ", self.space_group
        print "    Schoenflies ", self.schoenflies
        print "    Hall        ", self.hall
        print "    H-M symbol  ", self.hermann_mauguin
        print
        print "Laue  group ", self.laue_group
        print "Point group ", self.point_group

        print self.crystal_system
        if self.is_centrosymmetric:
            print "Centrosymmetric"

        if self.unique_axis:
            print "Unique axis", self.unique_axis
        print
        print "Order    ", self.order
        print "Order P  ", self.order_p

        print "\nCentering vectors ({})".format(self.centering)
        for cvec in self.centering_vectors:
            print cvec

        print "\nSymmetry operations"
        for symop in self.symmetry_operations:
            print symop
        if self.is_centrosymmetric:
            print "Inversion symmetry"
            for symop in self.symmetry_operations:
                print symop.inverse()

        print "\nReflection conditions"
        for rc in self.reflection_conditions:
            print rc

    def is_absent(self, index):
        """Expects tuple/list of 3 elements

        return True/False"""
        return any(c.is_absent(index) for c in self.reflection_conditions)
    
    def is_absent_np(self, index):
        """Efficient function to run on (n by 3) numpy arrays

        Return boolean array"""
        h = index[:, 0]
        k = index[:, 1]
        l = index[:, 2]
        return np.any([c.is_absent((h, k, l)) for c in self.reflection_conditions], axis=0)

    def is_absent_pd(self, index):
        """Efficient function to run on pandas index objects

        Return boolean array"""
        return index.map(self.is_absent)

    def unique_set(self, index):
        """Take list of reflections, use laue symmetry to merge to unique set"""
        raise NotImplementedError

    def filter_systematic_absences(self, df):
        """Takes a reflection list and filters reflections that are absent"""
        conditions = self.reflection_conditions
        try:
            index = df.index
        except AttributeError:
            index = df
        sel = self.is_absent_pd(index)
        return df[~sel]  # use binary not operator ~

    def merge(self, df):
        return merge(df, self)

    def completeness(self, df):
        return completeness(df, self)


def is_absent(index, conditions):
    """Expects tuple/list of 3 elements

    return True/False"""
    return any(c.is_absent(index) for c in conditions)


def is_absent_np(index, conditions):
    """Efficient function to run on (n by 3) numpy arrays

    Return boolean array"""
    h = index[:, 0]
    k = index[:, 1]
    l = index[:, 2]
    return np.any([c.is_absent((h, k, l)) for c in conditions], axis=0)


def filter_systematic_absences(index, conditions):
    """Takes a reflection list and filters reflections that are absent"""
    sel = is_absent_np(index, conditions)
    return index[~sel]  # use binary not operator ~


def test_sysabs_against_cctbx():
    """Test reflection conditions for all space groups, and verify
    them with CCTBX"""
    from cctbx import sgtbx
    from sysabs import refset
    import time

    time_numpy = 0
    time_cctbx_np = 0
    time_python = 0
    time_cctbx_loop = 0

    for i, row in enumerate(spacegrouptxt):
        print i

        number = row[:12].strip()
        spgr = get_spacegroup_info(number)

        sg = sgtbx.space_group(spgr["hall"])
        is_sys_absent = sg.is_sys_absent

        refconds = spgr["reflection_conditions"]
        refconds = [ReflCond(s) for s in refconds]

        t1 = time.time()
        a1 = is_absent_np(refset, refconds)
        t2 = time.time()
        a2 = np.apply_along_axis(is_sys_absent, axis=1, arr=refset)
        t3 = time.time()

        if not a1.shape == a2.shape:
            print refconds

        assert np.all(a1 == a2)
        assert sum(a1) == sum(a2)

        time_numpy += t2-t1
        time_cctbx_np += t3-t2

        t1 = time.time()
        a3 = [is_absent(ref, refconds) for ref in refset]
        t2 = time.time()
        a4 = [is_sys_absent(ref) for ref in refset]
        t3 = time.time()

        assert np.all(np.array(a3) == np.array(a4))
        assert sum(a3) == sum(a4)

        time_python += t2-t1
        time_cctbx_loop += t3-t2

    print "Numpy: {} s".format(time_numpy)
    print "CCTBX np.apply_along_axis: {} s".format(time_cctbx_np)
    print "Pure Python: {} s".format(time_python)
    print "CCTBX loop: {} s".format(time_cctbx_loop)


def test_print_all():
    """Parse and print all space groups"""
    for i, row in enumerate(spacegrouptxt):
        print "====="*16
        number = row[:12].strip()
        spgr = get_spacegroup_info(number)
        spgr.info()


def generate_hkl_listing_old(cell, dmin=1.0):
    """Generate hkllisting up to the specified dspacing.

    Based on the routine described by:
    Le Page and Gabe, J. Appl. Cryst. 1979, 12, 464-466
    """

    # if not cell.is_centrosymmetric:
    #     raise RuntimeError("Only centric structures can be used")

    lauegr = cell.laue_group
    system = cell.crystal_system
    uniq_axis = cell.unique_axis
    reflconds = cell.reflection_conditions
    setting = cell.setting

    if system == "Triclinic":  # "-1"
        segments = np.array([[[ 0, 0, 0], [ 1, 0, 0], [ 0, 1, 0], [ 0, 0, 1]],
                             [[-1, 0, 1], [-1, 0, 0], [ 0, 1, 0], [ 0, 0, 1]],
                             [[-1, 1, 0], [-1, 0, 0], [ 0, 1, 0], [ 0, 0,-1]],
                             [[ 0, 1,-1], [ 1, 0, 0], [ 0, 1, 0], [ 0, 0,-1]]])
    elif system == "Monoclinic":  # "2/m"
        if uniq_axis == "x":
            segments = np.array([[[ 0, 0, 0], [ 0, 1, 0], [1, 0, 0], [ 0, 0, 1]],
                                 [[ 0,-1, 1], [ 0,-1, 0], [1, 0, 0], [ 0, 0, 1]]])
        elif uniq_axis == "y":
            segments = np.array([[[ 0, 0, 0], [ 1, 0, 0], [ 0, 1, 0], [ 0, 0, 1]],
                                 [[-1, 0, 1], [-1, 0, 0], [ 0, 1, 0], [ 0, 0, 1]]])
        elif uniq_axis == "z":
            segments = np.array([[[ 0, 0, 0], [ 1, 0, 0], [ 0, 0, 1], [ 0, 1, 0]],
                                 [[-1, 1, 0], [-1, 0, 0], [ 0, 0, 1], [ 0, 1, 0]]])
    elif system == "Orthorhombic":  # mmm
        segments = np.array([[[ 0, 0, 0], [ 1, 0, 0], [ 0, 1, 0], [ 0, 0, 1]]])
    elif system == "Tetragonal":
        if lauegr == "4/mmm":
            segments = np.array([[[ 0, 0, 0], [ 1, 0, 0], [ 1, 1, 0], [ 0, 0, 1]]])
        elif lauegr == "4/m":
            segments = np.array([[[ 0, 0, 0], [ 1, 0, 0], [ 1, 1, 0], [ 0, 0, 1]],
                                 [[ 1, 2, 0], [ 1, 1, 0], [ 0, 1, 0], [ 0, 0, 1]]])
    elif system == "Trigonal":
        if setting == "R":
            if lauegr == "-3m":
                segments = np.array([[[ 0, 0, 0], [ 1, 0, 0], [ 1, 0,-1], [ 1, 1, 1]],
                                     [[ 1, 1, 0], [ 1, 0,-1], [ 0, 0,-1], [ 1, 1, 1]]])
            elif lauegr == "-3":
                segments = np.array([[[ 0, 0, 0], [ 1, 0, 0], [ 1, 0,-1], [ 1, 1, 1]],
                                     [[ 1, 1, 0], [ 1, 0,-1], [ 0, 0,-1], [ 1, 1, 1]],
                                     [[ 0,-1,-2], [ 1, 0, 0], [ 1, 0,-1], [-1,-1,-1]],
                                     [[ 1, 0,-2], [ 1, 0,-1], [ 0, 0,-1], [-1,-1,-1]]])
        else:
            if lauegr == "-3m1":
                segments = np.array([[[ 0, 0, 0], [ 1, 0, 0], [ 1, 1, 0], [ 0, 0, 1]],
                                     [[ 0, 1, 1], [ 0, 1, 0], [ 1, 1, 0], [ 0, 0, 1]]])
            elif lauegr == "-31m":
                segments = np.array([[[ 0, 0, 0], [ 1, 0, 0], [ 1, 1, 0], [ 0, 0, 1]],
                                     [[ 1, 1,-1], [ 1, 0, 0], [ 1, 1, 0], [ 0, 0,-1]]])
            elif lauegr == "-3":
                segments = np.array([[[ 0, 0, 0], [ 1, 0, 0], [ 1, 1, 0], [ 0, 0, 1]],
                                     [[ 1, 2, 0], [ 1, 1, 0], [ 0, 1, 0], [ 0, 0, 1]],
                                     [[ 0, 1, 1], [ 0, 1, 0], [-1, 1, 0], [ 0, 0, 1]]])
    elif system == "Hexagonal":
        if lauegr == "6/mmm":
            segments = np.array([[[ 0, 0, 0], [ 1, 0, 0], [ 1, 1, 0], [ 0, 0, 1]]])
        elif lauegr == "6/m":
            segments = np.array([[[ 0, 0, 0], [ 1, 0, 0], [ 1, 1, 0], [ 0, 0, 1]],
                                 [[ 1, 2, 0], [ 0, 1, 0], [ 1, 1, 0], [ 0, 0, 1]]])
    elif system == "Cubic":
        # TODO: Paper states lauegr m3m and 3m, typo? difference to m-3m/-3m??
        # if lauegr == "m3m":
        if lauegr == "m-3m":
            segments = np.array(
                [[[ 0, 0, 0], [ 1, 0, 0], [ 1, 1, 0], [ 1, 1, 1]]])
        # elif lauegr == "m3":
        elif lauegr == "m-3":
            segments = np.array([[[ 0, 0, 0], [ 1, 0, 0], [ 1, 1, 0], [ 1, 1, 1]],
                                 [[ 1, 2, 0], [ 0, 1, 0], [ 1, 1, 0], [ 1, 1, 1]]])
    else:
        raise ValueError, "Could not find crystal system {}".format(system)

    indices = []

    for row in segments:
        loop_h = True
        loop_k = True
        loop_l = True
        apex = row[0, :]

        index_stored = apex

        if sum(np.abs(apex)) == 0:
            dsp = 0.0  # prevent RuntimeWarning divide by zero
        else:
            dsp = cell.calc_dspacing(apex)

        while loop_l:
            while loop_k:
                while loop_h:
                    if dsp >= dmin:
                        indices.append(list(apex))
                    index_new = apex + row[1, :]

                    dsp = cell.calc_dspacing(index_new)

                    if dsp >= dmin:
                        apex = index_new
                    else:
                        loop_h = False

                apex[0] = index_stored[0]
                apex += row[2, :]
                index_new = apex
                dsp = cell.calc_dspacing(index_new)
                if dsp < dmin:
                    loop_k = False
                loop_h = True

            apex[1] = index_stored[1]
            apex += row[3, :]
            index_new = apex
            dsp = cell.calc_dspacing(index_new)
            if dsp < dmin:
                loop_l = False
            loop_k = True

    indices = np.array(indices)
    assert indices.dtype.kind == 'i', "Wrong datatype {}, need 'i'".format(x.dtype.kind)
    return indices


def generate_hkl_listing(cell, dmin=1.0, full_sphere=False, as_type='pd.Index'):
    """Generate hkllisting up to the specified dspacing.

    Based on the routine described by:
    Le Page and Gabe, J. Appl. Cryst. 1979, 12, 464-466
    """

    # if not cell.is_centrosymmetric:
    #     raise RuntimeError("Only centric structures can be used")

    lauegr = cell.laue_group
    system = cell.crystal_system
    uniq_axis = cell.unique_axis
    reflconds = cell.reflection_conditions
    setting = cell.setting

    segments = np.array([[[ 0, 0, 0], [ 1, 0, 0], [ 0, 1, 0], [ 0, 0, 1]],
                         [[-1, 0, 1], [-1, 0, 0], [ 0, 1, 0], [ 0, 0, 1]],
                         [[-1, 1, 0], [-1, 0, 0], [ 0, 1, 0], [ 0, 0,-1]],
                         [[ 0, 1,-1], [ 1, 0, 0], [ 0, 1, 0], [ 0, 0,-1]]])

    if system == "Trigonal":
        if setting == "R":
            lauegr += ":R"
        else:
            lauegr += ":H"

    if lauegr == "2/m":
        if uniq_axis == "y":
            lauegr += ":b"
        elif uniq_axis == "z":
            lauegr += ":c"
        elif uniq_axis == "x":
            lauegr += ":a"

    func = laue_asymmetry[lauegr]
    fraction = laue_fraction[lauegr]

    if full_sphere:
        func = lambda (h,k,l): True

    total = 0
    indices = []

    for row in segments:
        loop_h = True
        loop_k = True
        loop_l = True
        apex = row[0, :]

        index_stored = apex

        if sum(np.abs(apex)) == 0:
            dsp = np.inf  # prevent RuntimeWarning divide by zero
        else:
            dsp = cell.calc_dspacing(apex)

        while loop_l:
            while loop_k:
                while loop_h:
                    if dsp >= dmin:
                        if func(apex):
                            indices.append(list(apex))
                        else:
                            total += 1
                        if func(-apex):
                            indices.append(list(-apex))
                        else:
                            total += 1
                    index_new = apex + row[1, :]

                    dsp = cell.calc_dspacing(index_new)

                    if dsp >= dmin:
                        apex = index_new
                    else:
                        loop_h = False

                apex[0] = index_stored[0]
                apex += row[2, :]
                index_new = apex
                dsp = cell.calc_dspacing(index_new)
                if dsp < dmin:
                    loop_k = False
                loop_h = True

            apex[1] = index_stored[1]
            apex += row[3, :]
            index_new = apex
            dsp = cell.calc_dspacing(index_new)
            if dsp < dmin:
                loop_l = False
            loop_k = True

    indices = np.array(indices[2:]) # ignore double [0, 0, 0]
    # print indices
    # print "Generated {} indices, {} in asymmetric unit for laue group {} ({:.1f}%)".format(total, len(indices), lauegr, 100*len(indices)/float(total))
    # print "Kept {} / {} indices for lauegr {} ({:.1f}%). Fraction {}".format(len(indices), 
                                                                             # len(indices)+total,
                                                                             # lauegr, 
                                                                             # 100.0*float(len(indices))/(len(indices)+total), 
                                                                             # fraction)

    assert indices.dtype.kind == 'i', "Wrong datatype {}, need 'i'".format(indices.dtype.kind)
    if as_type == "pd.Index":
        return pd.Index([tuple(idx) for idx in indices])
    if as_type == "pd.DataFrame":
        return pd.DataFrame(index=[tuple(idx) for idx in indices])
    else:
        return indices

def get_laue_symops(key):
    from laue_symops import symops
    return (SymOp(op) for op in symops[key])

def get_merge_dct(df, cell):
    if 'd' not in df:
        df['d'] = df.index.map(cell.calc_dspacing)

    dmin = df['d'].min()
    unique_set = generate_hkl_listing(cell, dmin=dmin)

    lauegr = cell.laue_group
    setting = cell.setting
    uniq_axis = cell.unique_axis

    if lauegr == "2/m":
        if uniq_axis == "y":
            lauegr += ":b"
        elif uniq_axis == "z":
            lauegr += ":c"
        elif uniq_axis == "x":
            lauegr += ":a"
    symops = get_laue_symops(lauegr)
    
    merge_dct = {}
    for op in symops:
        for idx in unique_set:
            new = tuple(np.dot(idx, op.r))
            if new not in merge_dct:
                merge_dct[new] = tuple(idx)
    return merge_dct

def expand(df, cell):
    lauegr = cell.laue_group
    setting = cell.setting
    uniq_axis = cell.unique_axis

    if lauegr == "2/m":
        if uniq_axis == "y":
            lauegr += ":b"
        elif uniq_axis == "z":
            lauegr += ":c"
        elif uniq_axis == "x":
            lauegr += ":a"
    symops = get_laue_symops(lauegr)

    ret = []
    for op in symops:
        for idx in df.index:
            new = tuple(np.dot(idx, op.r))
            if new not in ret:
                ret.append(idx)
    return ret

def merge(df, cell):
    merge_dct = get_merge_dct(df, cell)

    new = df.groupby(merge_dct).mean()

    print "Merged {} to {} reflections".format(len(df), len(new))
    return new

def completeness(df, cell, dmin=None):
    if not dmin:
        if 'd' not in df:
            df['d'] = df.index.map(cell.calc_dspacing)
        dmin = df['d'].min()

    unique_set = generate_hkl_listing(cell, dmin=dmin)
    unique_set = cell.filter_systematic_absences(unique_set)

    len_a = (df['d'] > dmin).sum()
    len_b = len(unique_set)

    return len_a / float(len_b)


if __name__ == '__main__':
    # main()

    # test_sysabs_against_cctbx()
    # test_print_all()

    # exit()

    # arg = "P2/c:a1"
    arg = "random"

    spgr = get_spacegroup_info(arg)

    from unitcell import UnitCell
    cell = get_random_cell(spgr)
    cell = UnitCell(cell, spgr.space_group)

    cell.info()

    z = generate_hkl_listing(cell, dmin=1.0, as_type="pd.DataFrame")
    full = generate_hkl_listing(cell, dmin=1.0, full_sphere=True, as_type="pd.DataFrame")
    new = merge(full, cell)

    print "\nCompleteness {}".format(completeness(z, cell))
    print
    print "\nCompleteness {}".format(completeness(full, cell))

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
