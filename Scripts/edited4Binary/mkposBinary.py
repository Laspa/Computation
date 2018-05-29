#!/usr/bin/env python
# ============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  June 1, 2016
#  Script to make simple POTCAR files
# ============================================================================
from Cell2 import *
import math

# ============================================================================
# POSCAR Building
# ============================================================================


def makePOSCAR(compStr, compoundL, structure, a0, style=5, dest='.'):
    """ Make compound POSCAR file with specified structure """
    header = '%s_%s' % (compStr, structure)
    name = 'POSCAR_%s' % header
    print 'Making %s...\n' % name
    if structure == 'FCC':
        makeFCCPOSCAR(header, compoundL, a0, style, dest, name)
    elif structure == 'BCC':
        makeBCCPOSCAR(header, compoundL, a0, style, dest, name)
    elif structure == 'DHCP':
        makeDHCPPOSCAR(header, compoundL, a0, style, dest, name)
    else:
        makeHCPPOSCAR(header, compoundL, a0, style, dest, name)


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


def makeBCCPOSCAR(header, compoundL, a0, style=5, dest='.', name='POSCAR'):
    """ Make pure element BCC POSCAR """
    a = [1, 0, 0]
    b = [0, 1, 0]
    c = [0, 0, 1]
    cell = Cell()
    cell.setHeader(header)
    cell.setElements(compoundL)
    cell.setA0(a0)  # overall scaling
    cell.setCoordinateSystem('Direct')  # change to direct
    cell.setLatticeVectors([a, b, c])
    pos = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    for i in range(0, len(pos)):
        cell.addSite(i, Site(pos[i]))  # add sites
    cell.sendToPOSCAR('%s/%s' % (dest, name), style=style)


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
# Get element from user
temp = str(raw_input('Element: '))
compoundL = [x.strip() for x in temp.split(',')]
compStr = temp.replace(", ", "")
print compStr, '\n'

# Get structure from user
structure = raw_input('Structure (FCC, BCC, HCP): ')
structure = structure.strip().upper()
if structure and structure[0] == 'B':
    structure = 'BCC'
elif structure and structure[0] == 'H':
    structure = 'HCP'
elif structure and structure[0] == 'D':
    structure = 'DHCP'
else:
    structure = 'FCC'
print structure, '\n'

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
makePOSCAR(compStr, compoundL, structure, a0, style)