#!/usr/bin/env python

from spgr import *
from unitcell import *

import numpy as np

import argparse

__version__ = "2015-12-10"


def main():
    description = """"""

    epilog = 'Updated: {}'.format(__version__)

    parser = argparse.ArgumentParser(description=description,
                                     epilog=epilog,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser = argparse.ArgumentParser()

    parser.add_argument("args",
                        type=str, metavar="spgr", nargs='*',
                        help="Space group")

    parser.add_argument("-l", "--hkllisting",
                        action="store_true", dest="hkllisting",
                        help="Print hkl listing to file out.hkl")

    parser.add_argument("-d", "--dmin",
                        action="store", type=float, dest="dmin",
                        help="Resolution to generate hkl indices with (default = 1.0)")

    parser.add_argument("-k", "--maxhkl",
                        action='store', type=int, dest="maxhkl",
                        help="Maximum hkl values to generate indices with")

    parser.add_argument("-c", "--cell",
                        action='store', type=float, nargs="+", dest="cell",
                        help="Unit cell")

    parser.add_argument("-s", "--spgr",
                        action="store", type=str, dest="spgr",
                        help="Space group to use. Tell the program to interpret args as filenames, and use one of the commands below.")

    parser.add_argument("-m", "--merge",
                        action="store_true", dest="merge",
                        help="Merge symmetry equivalent reflectiosn in data files")

    parser.add_argument("-p", "--completeness",
                        action="store_true", dest="completeness",
                        help="Calculate completeness of given data files")

    parser.add_argument("-f", "--filter",
                        action="store_true", dest="filter",
                        help="Filter systematic absences from data files.")

    parser.set_defaults(hkllisting=False,
                        dmin=1.0,
                        maxhkl=10,
                        cell=None,
                        # file specific
                        spgr=None,
                        merge=False,
                        completeness=False,
                        filtersymequiv=False)

    options = parser.parse_args()
    args = options.args

    if not options.spgr:
        for arg in args:
            spgr = get_spacegroup_info(arg)
            if not spgr:
                continue
            print "# {}\n".format(arg)
            spgr.info()

            if options.hkllisting:
                if options.cell:
                    cell = options.cell
                    if cell[0] < 0:
                        cell = get_random_cell(spgr)
                    cell = get_unitcell(cell, spgr.space_group)
                    print "\nUnit cell:", cell
                elif options.maxhkl:
                    # bit hacky, but works! :)
                    cell = (1, 1, 1, 90, 90, 90)
                    cell = get_unitcell(cell, spgr.space_group)
                    options.dmin = 1/float(options.maxhkl)

                indices = generate_hkl_listing(cell, dmin=options.dmin)
                if len(args) > 1:
                    out = 'out_{}.hkl'.format(cell.spgr)
                else:
                    out = 'out.hkl'
                np.savetxt(out, indices, fmt="%4d")

    if options.spgr:
        for arg in args:
            if options.merge:
                pass
            if options.completeness:
                pass
            if options.filter:
                pass

if __name__ == '__main__':
    main()
