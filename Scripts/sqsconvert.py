#!/usr/bin/env python
#=========================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  January 31, 2016
#  Script to batch convert bestsqs.out to POSCARs
#=========================================================================
from Cell import *

# Load from SQS
sqsFile = raw_input('Path to structure file: ')

scale = raw_input('Scale by first coordinate vector? (y/n): ')
if 'Y' in scale.upper():
    scale = True
else:
    scale = False
cell = Cell().loadFromSQS(sqsFile, scale=scale)

# Reassign elements, sort
print 'Elements are ' + ','.join(cell.elements)
elements = raw_input('New elements (press ENTER to keep): ')
if elements:
	elements = elements.split(',')
	cell.setElements(elements)
	cell.sortElements()

# Set a0

a0 = raw_input('Scaling factor: ')
if a0:
	cell.setA0(float(a0))

# Send to POSCAR
style = raw_input('VASP style 4 or 5? (default 5): ')
if style and style == '4':
    style = int(style)
else:
    style = 5
cell.sendToPOSCAR('POSCAR', style=style)
