#!/usr/bin/env python
# ==============================
#  Emily Hwang
#  May 28, 2018
#  Script to make a simple INCAR
# ==============================
"""
The script generates a simple and generic INCAR file to the working directory.
"""
from Cell import *

def makeINCAR(dest='.'):
    """ Creates a simple INCAR file in dest """
    s = 'ENCUT = %s' % ENCUT
    s += '\nPREC = Accurate'
    s += '\nISMEAR = 1'
    s += '\nSIGMA = 0.1'
    s += '\nNSW = 51'
    s += '\nIBRION = %s' % IBRION
    s += '\nISIF = %s' % ISIF
    s += '\nNPAR = %s' % NPAR
    s += '\nKPAR = %s' % KPAR
    s += '\nLWAVE = .FALSE.'
    s += '\nLCHARG = .FALSE.'
    
    f = open('INCAR', 'w+')
    f.write(s)
    f.close()


# ========================
#  Main Program
# ========================

# Get parameters from user
ENCUT = raw_input('ENCUT: ')
print ENCUT, '\n'

valid = False
while not valid:
    IBRION = int(raw_input('IBRION (default 2): '))
    if not IBRION:
        IBRION = 2
        valid = True
    elif IBRION < 9 or IBRION == 44:
        IBRION = int(IBRION)
        valid = True
    else:
        print 'Not a valid value for IBRION'
print IBRION,'\n'

correct = False
while not correct:    
    ISIF = int(raw_input('ISIF (default 3): '))
    if not ISIF:
        ISIF = 3
        correct = True
    elif ISIF < 8:
        ISIF = int(ISIF)
        correct = True
    else:
        print 'Not a valid value for ISIF'
print ISIF,'\n'

NPAR = raw_input('NPAR: ')
print NPAR, '\n'

KPAR = raw_input('KPAR: ')
print KPAR, '\n', '\n'

print 'Creating INCAR file...\n' 
makeINCAR()
