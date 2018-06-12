#!/usr/bin/env python
#===============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu 
#  June 9, 2016
#  Script to make KPOINTS file
#  Edited by Emily Hwang eyhwang@hmc.edu
# ============================================================================
import subprocess as sp
import numpy as np
from Cell import *

# ============================================================================
# KPOINTS Generation
# ============================================================================


def KPOINTS(points, header='', dest='.', gamma=True):
    """" Make KPOINTS in dest file with specified number of KPOINTS, center """
    print 'Making KPOINTS file...'
    if gamma:
        center = 'Gamma'
    else:
        center = 'Monkhorst'
    header = str(header)
    
    s = 'Automatic mesh %s' % header
    s += '\n0' 
    s += '\n%s' % center
    s += '\n%d %d %d' % (points,points,points)
    s += '\n0 0 0'
   
    f = open('%s/KPOINTS' % dest, 'w+')
    f.write(s)
    f.close()


# ============================================================================
#  Main Program
# ============================================================================
center = raw_input('Center (Gamma or Monkhorst): ')
gamma = True
if center and center[0].upper() == 'M':
    center = 'Monkhorst'
    gamma = False
else:
    center = 'Gamma'
print center, '\n'

points = int(raw_input('Number of kpoints: '))
print points, '\n'

header = str(raw_input('Header: '))

# Make file
KPOINTS(points, header, gamma=gamma)
print
