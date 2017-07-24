# xcore
Xcore is a crystallographic space group library written in Python.

## Usage

    >>> from xcore import UnitCell, SpaceGroup
    
    >>> spgr = SpaceGroup("Pnma")
    >>> spgr.info()
    Space group Pnma

        Number       62
        Schoenflies  D2h^16
        Hall         -P 2ac 2n
        H-M symbol   Pnma

    Laue  group  mmm
    Point group  mmm
    Orthorhombic
    Centrosymmetric

    Order     8
    Order P   8

    # +(0.0 0.0 0.0), Inversion Flag = 0
    x,y,z
    -x+1/2,-y,z+1/2
    x+1/2,-y+1/2,-z+1/2
    -x,y+1/2,-z
    # +(0.0 0.0 0.0), Inversion Flag = 1
    -x,-y,-z
    x+1/2,y,-z+1/2
    -x+1/2,y+1/2,z+1/2
    x,-y+1/2,z
    
or:
    
    >>> cell = UnitCell((19.9020, 10.1561, 24.6892, 105.88), "P21/c", name="SSO", composition="Si80 O163")
    >>> cell.info()
    Cell SSO (Si80 O163)
       a      19.9020       al        90.00
       b      10.1561       be       105.88
       c      24.6892       ga        90.00
    Vol.    4799.90
    Spgr P121/c1


## Installation

    pip install xcore

## Requirements
    
 - Python2.7
 - Numpy
 - Pandas

Inspired by (and built on) [sginfo](http://cci.lbl.gov/sginfo/).