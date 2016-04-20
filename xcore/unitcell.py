#!/usr/bin/env python

import numpy as np
from math import radians, cos, sin

try:
    # raise ImportError
    import tinyarray as ta
except ImportError:
    import numpy as ta
    TINYARRAY = False
else:
    TINYARRAY = True

import spacegroup
from spacegroup import SpaceGroup, get_spacegroup_info


from IPython.terminal.embed import InteractiveShellEmbed
InteractiveShellEmbed.confirm_exit = False
ipshell = InteractiveShellEmbed(banner1='')


def get_unitcell(parameters, spgr):
    spgr = get_spacegroup_info(spgr, as_dict=True)
    cell = UnitCell(parameters, spgr)
    return cell


def dict2uc(d):
    cell_params = [d[key] for key in "a", "b", "c", "al", "be", "ga"]
    spgr = d["spgr"]
    name = d.get("name", "")
    return UnitCell(cell_params=cell_params, spgr=spgr, name=name)


class UnitCell(SpaceGroup):

    """Class for unit cell/space group functions"""

    def __init__(self, cell_params, spgr, name=""):
        if isinstance(spgr, str):
            spgr = get_spacegroup_info(spgr, as_dict=True)
        super(UnitCell, self).__init__(spgr)

        self.name = name
        
        if len(cell_params) != 6:
            cell_params = self.parse_cellparams(cell_params)
        
        self.parameters = tuple(float(par) for par in cell_params)

        if not self.is_valid_cell():
            print "\n >> Warning: Unit cell parameters do not fit with space group {}".format(self.space_group)

    def __repr__(self):
        if self.name:
            return "{}: {} - {}".format(self.name, str(self.parameters), self.spgr_name)
        else:
            return "{} - {}".format(str(self.parameters), self.spgr_name)

    def __iter__(self):
        for par in self.parameters:
            yield par

    @property
    def a(self):
        return self.parameters[0]

    @property
    def b(self):
        return self.parameters[1]

    @property
    def c(self):
        return self.parameters[2]

    @property
    def al(self):
        return self.parameters[3]

    @property
    def be(self):
        return self.parameters[4]

    @property
    def ga(self):
        return self.parameters[5]

    def as_dict(self):
        return {"name": self.name, 
                "spgr": self.spgr_name, 
                "a": self.a, 
                "b": self.b, 
                "c": self.c, 
                "al": self.al, 
                "be": self.be, 
                "ga": self.ga }

    def info(self):
        print "Cell {}".format(self.name)
        print "   a {:12.4f}       al {:12.2f}".format(self.a, self.al)
        print "   b {:12.4f}       be {:12.2f}".format(self.b, self.be)
        print "   c {:12.4f}       ga {:12.2f}".format(self.c, self.ga)
        print "Vol. {:10.2f}".format(self.volume)
        print "Spgr {}".format(self.spgr_name)
        print


    def metric_tensor(self, inverse=False):
        """Returns the metric tensor

        http://reference.iucr.org/dictionary/Metric_tensor

        Dunitz, 1979, p227"""

        a, b, c, al, be, ga = self.parameters

        al = radians(al)
        be = radians(be)
        ga = radians(ga)

        vol = self.volume

        if inverse:
            m11 = (b*c*sin(al)/vol)**2
            m22 = (c*a*sin(be)/vol)**2
            m33 = (a*b*sin(ga)/vol)**2

            m12 = a*b*(c/vol)**2 * (cos(al)*cos(be)-cos(ga))
            m23 = b*c*(a/vol)**2 * (cos(be)*cos(ga)-cos(al))
            m13 = a*c*(b/vol)**2 * (cos(ga)*cos(al)-cos(be))

            mati = ta.array([[m11, m12, m13],
                             [m12, m22, m23],
                             [m13, m23, m33]])
        else:
            mat = ta.array([[a*a,         a*b*cos(ga), a*c*cos(be)],
                            [a*b*cos(ga),         b*b, b*c*cos(al)],
                            [a*c*cos(be), b*c*cos(al),         c*c]])

        return mat

    def orthogonalization_matrix(self, inverse=False):
        """orthogonalization matrix for crystal to cartesian coordinates
        not to be confused with the unit cell orthogonalization matrix, which is the transpose of this one

        Dunitz convention, Dunitz, 1979, p237"""

        a, b, c, al, be, ga = self.parameters

        al = radians(al)
        be = radians(be)
        ga = radians(ga)

        vol = self.volume

        if inverse:
            mat = ta.array([[1/a, (-1*cos(ga)) / (a*sin(ga)), (cos(ga) * cos(al) - cos(be)) / (a*vol * sin(ga))],
                            [0,            1 /
                                (b*sin(ga)), (cos(ga) * cos(be) - cos(al)) / (b*vol * sin(ga))],
                            [0,                          0,                 (a*b*sin(ga)) / (vol)]])
        else:
            mat = ta.array([[a, b*cos(ga),                           c*cos(be)],
                            [0, b*sin(ga),
                             c*(cos(al)-cos(be)*cos(ga))/sin(ga)],
                            [0,         0,                   vol/(a*b*sin(ga))]])

        return mat

    def calc_dspacing(self, idx):
        """Calc dspacing at given index (i.e. idx= (1,0,0)

        Calculates d-spacing based on given parameters.
        a,b,c,al,be,ge are given as floats
        al,be,ga can be given as ndarrays or floats
        kind specifies the type of cell -> triclinic works for the general case, but is a bit slower
        although still fast enough

        Tested: orthorhombic cell on (orthorhombic, monoclinic, triclinic)
        Tested: triclinic cell with dvalues from topas
        """

        kind = self.crystal_system
        a, b, c, al, be, ga = self.parameters
        h = idx[0]
        k = idx[1]
        l = idx[2]

        if kind == 'Cubic':
            idsq = (h**2 + k**2 + l**2) / a**2

        elif kind == 'Tetragonal':
            idsq = (h**2 + k**2) / a**2 + l**2 / c**2

        elif kind == 'Orthorhombic':
            idsq = h**2 / a**2 + k**2 / b**2 + l**2 / c**2

        elif kind == "Trigonal":
            if self.setting == "R":
                al = radians(al)
                num = (h**2 + k**2 + l**2) * sin(al)**2 + 2*(h*k + k*l + h*l)*(cos(al)**2 - cos(al))
                denom = a**2 * (1 - 3*cos(al)**2 + 2*cos(al)**3)
                idsq = num / denom
            else:
                idsq = (4.0/3.0) * (h**2 + h*k + k**2) / (a**2) + l**2 / c**2

        elif kind == 'Hexagonal':
            idsq = (4.0/3.0) * (h**2 + h*k + k**2) / (a**2) + l**2 / c**2

        elif kind == 'Monoclinic':
            be = radians(be)
            idsq = (1/sin(be)**2) * (h**2/a**2 + k**2 * sin(be)**2 /
                                     b**2 + l**2/c**2 - (2*h*l*cos(be)) / (a*c))

        elif kind == 'Triclinic':
            V = self.volume

            al = radians(al)
            be = radians(be)
            ga = radians(ga)

            idsq = (1/V**2) * (
                h**2 * b**2 * c**2 * sin(al)**2
                + k**2 * a**2 * c**2 * sin(be)**2
                + l**2 * a**2 * b**2 * sin(ga)**2
                + 2*h*k*a*b*c**2 * (cos(al) * cos(be) - cos(ga))
                + 2*k*l*b*c*a**2 * (cos(be) * cos(ga) - cos(al))
                + 2*h*l*c*a*b**2 * (cos(al) * cos(ga) - cos(be))
            )
        else:
            raise ValueError("Unknown crystal system {}, fallback to Triclinic".format(kind))

        # if idsq == 0:
            # prevent RuntimeWarning: divide by zero
            # return np.inf
        # else:
        return idsq**-0.5

    def calc_dspacing_np(self, idx):
        h = idx[:,0]
        k = idx[:,1]
        l = idx[:,2]
        return self.calc_dspacing((h,k,l))

    @property
    def volume(self):
        """Returns volume for the general case from cell parameters"""
        if hasattr(self, "_volume"):
            return self._volume
        a, b, c, al, be, ga = self.parameters
        al = radians(al)
        be = radians(be)
        ga = radians(ga)
        vol = a*b*c * \
            ((1+2*cos(al)*cos(be)*cos(ga)-cos(al)**2-cos(be)**2-cos(ga)**2)
             ** .5)
        self._volume = vol
        return vol


    def is_valid_cell(self):
        a,b,c,al,be,ga = self.parameters
        system = self.crystal_system
        setting = self.setting
        if system == "Triclinic":
            return True
        elif system == "Monoclinic":
            if self.unique_axis == "y":
                return al == ga == 90.0
            elif self.unique_axis == "x":
                return be == ga == 90.0
            elif self.unique_axis == "z":
                return al == be == 90.0
        elif system == "Orthorhombic":
            return al == be == ga
        elif system == "Tetragonal":
            return (a == b) and (al == be == ga == 90.0)
        elif system == "Trigonal":
            if setting == "R":
                return (a == b == c) and (al == be == ga)
            else:
                return (a == b) and (al == be == 90.0) and (ga == 120.0)
        elif system == "Hexagonal":
            return (a == b) and (al == be == 90.0) and (ga == 120.0)
        elif system == "Cubic":
            return (a == b == c) and (al == be == ga == 90.0)
        else:
            raise ValueError("Unknown crystal system ".format(system))

    def parse_cellparams(self, parameters):
        system = self.crystal_system
        if system == "Triclinic":
            assert len(parameters) == 6, "Expect 6 cell parameters"
        elif system == "Monoclinic":
            assert len(parameters) == 4, "Expect 4 cell parameters"
            a, b, c, angle = parameters
            if self.unique_axis == "y":
                parameters = [a, b, c, 90.0, angle, 90.0]
            elif self.unique_axis == "x":
                parameters = [a, b, c, angle, 90.0, 90.0]
            elif self.unique_axis == "z":
                parameters = [a, b, c, 90.0, 90.0, angle]
        elif system == "Orthorhombic":
            assert len(parameters) == 3, "Expect 3 cell parameters"
            a, b, c = parameters
            parameters = [a, b, c, 90.0, 90.0, 90.0]
        elif system == "Tetragonal":
            assert len(parameters) == 2, "Expect 2 cell parameters"
            a, c = parameters
            parameters = [a, a, c, 90.0, 90.0, 90.0]
        elif system == "Trigonal":
            if self.setting == "R":
                assert len(parameters) == 2, "Expect 2 cell parameters"
                a, al = parameters
                parameters = [a, a, a, al, al, al]
            else:
                assert len(parameters) == 2, "Expect 2 cell parameters"
                a, c = parameters
                parameters = [a, a, c, 90.0, 90.0, 120.0]
        elif system == "Hexagonal":
            assert len(parameters) == 2, "Expect 2 cell parameters"
            a, c = parameters
            parameters = [a, a, c, 90.0, 90.0, 120.0]
        elif system == "Cubic":
            assert len(parameters) == 1, "Expect 1 cell parameters"
            a = parameters[0]
            parameters = [a, a, a, 90.0, 90.0, 90.0]
        else:
            raise ValueError("Unknown crystal system ".format(system))

        assert len(parameters) == 6, "Expect 6 cell parameters"
        return parameters

    def get_dmin(self, indices):
        ipshell();exit()
        return np.min(self.calc_dspacing_np(indices))
if __name__ == '__main__':
    arg = "random"

    spgr =  spacegroup.get_spacegroup(arg)

    from unitcell import UnitCell
    cell = spacegroup.get_random_cell(spgr)
    cell = UnitCell(cell, spgr.space_group)

    cell.info()

    z = spacegroup.generate_hkl_listing(cell, dmin=1.0, as_type="pd.DataFrame")
    full = spacegroup.generate_hkl_listing(cell, dmin=1.0, full_sphere=True, as_type="pd.DataFrame")
    new = spacegroup.merge(full, cell)

    print "\nCompleteness {}".format(spacegroup.completeness(z, cell))
    print
    print "\nCompleteness {}".format(spacegroup.completeness(full, cell))
