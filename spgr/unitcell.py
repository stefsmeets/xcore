#!/usr/bin/env python

import numpy as np
from math import radians, cos, sin

from IPython.terminal.embed import InteractiveShellEmbed
InteractiveShellEmbed.confirm_exit = False
ipshell = InteractiveShellEmbed(banner1='')

class UnitCell(object):

    """Class for unit cell/space group functions"""

    def __init__(self, a, b, c, al=90.0, be=90.0, ga=90.0):
        self.parameters = (float(a), float(b), float(c), float(al), float(be), float(ga))

    def __repr__(self):
        return str(self.parameters)

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

    @profile
    def calc_dspacing(self, idx, kind="Triclinic"):
        """Calc dspacing at given index (i.e. idx= (1,0,0)

        Calculates d-spacing based on given parameters.
        a,b,c,al,be,ge are given as floats
        al,be,ga can be given as ndarrays or floats
        kind specifies the type of cell -> triclinic works for the general case, but is a bit slower
        although still fast enough

        Tested: orthorhombic cell on (orthorhombic, monoclinic, triclinic)
        Tested: triclinic cell with dvalues from topas
        """

        a, b, c, al, be, ga = self.parameters
        h = idx[0]
        k = idx[1]
        l = idx[2]

        if kind == 'Cubic':
            print '\n** Warning: cubic dspacing calculation unverified!! **\n'
            idsq = (h**2 + k**2 + l**2) / a**2

        elif kind == 'Tetragonal':
            print '\n** Warning: tetragonal dspacing calculation unverified!! **\n'
            idsq = (h**2 + k**2) / a**2 + l**2 / c**2

        elif kind == 'Orthorhombic':
            idsq = h**2 / a**2 + k**2 / b**2 + l**2 / c**2

        elif kind == 'Hexagonal':
            print '\n** Warning: hexagonal dspacing calculation unverified!! **\n'
            idsq = (4.0/3.0) * (h**2 + h*k + k**2) * (a**-2) + l**2 / c**2

        elif kind == 'Monoclinic':
            print '\n** Warning: monoclinic dspacing calculation unverified!! **\n'
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
            print "Unknown crystal system {}, fallback to Triclinic"
            return self.calc_dspacing(idx, kind="Triclinic")

        if idsq == 0:
            # prevent RuntimeWarning: divide by zero
            return np.inf
        else:
            return idsq**-0.5

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

