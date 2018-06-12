#!/usr/bin/env python
#===============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu 
#  June 1, 2016
#  Script to make simple POSCAR files
#  Edited by Emily Hwang eyhwang@hmc.edu
# ============================================================================
from Cell import *
import math

# ============================================================================
# POSCAR Building
# ============================================================================


def makePOSCAR(element, structure, a0, style=5, dest='.'):
    """ Make pure element POSCAR file with specified structure """
    header = '%s_%s' % (element, structure)
    name = 'POSCAR_%s' % header
    print 'Making %s...\n' % name
    if structure == 'FCC':
        makeFCCPOSCAR(header, element, a0, style, dest, name)
    elif structure == 'BCC':
        makeBCCPOSCAR(header, element, a0, style, dest, name)
    elif structure == 'DHCP':
        makeDHCPPOSCAR(header, element, a0, style, dest, name)
    else:
        makeHCPPOSCAR(header, element, a0, style, dest, name)

def makeBinaryPOSCAR(element1, element2, structure, a0, style=5, dest='.'):
    """Make a binary POSCAR file with specified structure """
    header = '%s_%s' % (element1+element2, structure)
    name = 'POSCAR_%s' % header
    print 'Making %s...\n' % name
    if structure == 'B2':
        makeB2POSCAR(header, element1, element2, a0, style, dest, name)

def makeFCCPOSCAR(header, element, a0, style=5, dest='.', name='POSCAR'):
    """ Make pure element FCC POSCAR """
    a = [0.0, 0.5, 0.5]
    b = [0.5, 0.0, 0.5]
    c = [0.5, 0.5, 0.0]
    cell = Cell()
    cell.setHeader(header)
    cell.setElements([element])
    cell.setA0(a0)  # overall scaling
    cell.setCoordinateSystem('Direct')  # change to direct
    cell.setLatticeVectors([a, b, c])
    for p in [[0.0, 0.0, 0.0]]:
        cell.addSite(0, Site(p))  # add sites
    cell.sendToPOSCAR('%s/%s' % (dest, name), style=style)


def makeBCCPOSCAR(header, element, a0, style=5, dest='.', name='POSCAR'):
    """ Make pure element BCC POSCAR """
    a = [0.5, 0.5, -0.5]
    b = [-0.5, 0.5, 0.5]
    c = [0.5, -0.5, 0.5]
    cell = Cell()
    cell.setHeader(header)
    cell.setElements([element])
    cell.setA0(a0)  # overall scaling
    cell.setCoordinateSystem('Direct')  # change to direct
    cell.setLatticeVectors([a, b, c])
    for p in [[0.0, 0.0, 0.0]]:
        cell.addSite(0, Site(p))  # add sites
    cell.sendToPOSCAR('%s/%s' % (dest, name), style=style)

def makeB2POSCAR(header,element1, element2,a0,style=5,dest=".",name='POSCAR'):
    """ Make binary B2 POSCAR """
    a = [1, 0, 0]
    b = [0, 1, 0]
    c = [0, 0, 1]
    cell = Cell()
    cell.setHeader(header) 
    cell.setElements([element1, element2])
    cell.setA0(a0)
    cell.setCoordinateSystem('Direct')
    cell.setLatticeVectors([a,b,c])
    for p in [[0.0,0.0,0.0], [0.5, 0.5, 0.5]]:
        cell.addSite(0,Site(p))
    cell.elementCounts = [1, 1]
    cell.sendToPOSCAR('%s/%s' % (dest,name), style=style)

def makeHCPPOSCAR(header, element, a0, style=5, dest='.', name='POSCAR'):
    """ Make pure element HCP POSCAR """
    a = [1, 0, 0]
    b = [-0.5, math.sqrt(3) / 2, 0]
    c = [0, 0, math.sqrt(8.0 / 3.0)]  # ideal c/a ratio
    cell = Cell()
    cell.setHeader(header)
    cell.setElements([element])
    cell.setA0(a0)  # overall scaling
    cell.setCoordinateSystem('Direct')  # change to direct
    cell.setLatticeVectors([a, b, c])
    for p in [[0.0, 0.0, 0.0], [1.0 / 3.0, 2.0 / 3.0, 0.5]]:
        cell.addSite(0, Site(p))  # add sites
    cell.sendToPOSCAR('%s/%s' % (dest, name), style=style)


def makeDHCPPOSCAR(header, element, a0, style=5, dest='.', name='POSCAR'):
    """ Make pure element DHCP POSCAR """
    a = [1, 0, 0]
    b = [-0.5, math.sqrt(3) / 2, 0]
    c = [0, 0, 2 * math.sqrt(8.0 / 3.0)]  # ideal c/a ratio
    cell = Cell()
    cell.setHeader(header)
    cell.setElements([element])
    cell.setA0(a0)  # overall scaling
    cell.setCoordinateSystem('Direct')  # change to direct
    cell.setLatticeVectors([a, b, c])
    for p in [[0.0, 0.0, 0.0], [1.0 / 3.0, 2.0 / 3.0, 0.25],
              [2.0 / 3.0, 1.0 / 3.0, 0.5], [1.0 / 3.0, 2.0 / 3.0, 0.75]]:
        cell.addSite(0, Site(p))  # add sites
    cell.sendToPOSCAR('%s/%s' % (dest, name), style=style)

# ============================================================================
#  Main Program
# ============================================================================
# Get number of elements from user
binary = raw_input('Binary structure? (Y/N): ')
binary = binary.strip().upper()

# Get structure from user
structure = str(raw_input('Structure (FCC, BCC, B2, HCP): '))
structure = structure.strip().upper()
if structure and structure[0:2] == 'BC':
    structure = 'BCC'
elif structure and structure[0:2] == 'B2':
    structure = 'B2'
elif structure and structure[0] == 'H':
    structure = 'HCP'
elif structure and structure[0] == 'D':
    structure = 'DHCP'
else:
    structure = 'FCC'
print structure, '\n'

# Get element from user
if binary == 'Y':
    element1 = str(raw_input('Element 1: '))
    print element1, '\n'
    element2 = str(raw_input('Element 2: '))
    print element2, '\n'
else:
    element = str(raw_input('Element: '))
    print element, '\n'

# Get lattice constant from user
a0 = float(raw_input('Lattice parameter: '))
print a0, '\n'

# Get VASP style
style = raw_input('VASP style 4 or 5? (default 5): ')
if style and style == '4':
    style = int(style)
else:
    style = 5
print style, '\n'

# Make POSCAR in PWD
if binary == 'Y':
    makeBinaryPOSCAR(element1, element2, structure, a0, style)
else: 
    makePOSCAR(element, structure, a0, style)
