#!/usr/bin/env python
# ============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  May 28, 2016
#  Script to generate files for making SQSs using screen
#  Supports running multiple independent processes of mcsqs
# ============================================================================
import math
import subprocess as sp
import os

# ============================================================================
# SFE Structure Building
# ============================================================================


def makeLat(structure, a0, concLine, dest='.'):
    """ Makes lat.in file for given structure """
    if structure == 'hcp':
        a = [i * a0 for i in [1, 0, 0]]
        b = [i * a0 for i in [-0.5, math.sqrt(3) / 2, 0]]
        c = [i * a0 for i in [0, 0, math.sqrt(8.0 / 3.0)]]
        u = [1, 0, 0]
        v = [0, 1, 0]
        w = [0, 0, 1]
        sites = [[0.0, 0.0, 0.0], [1.0 / 3.0, 2.0 / 3.0, 0.5]]
    elif structure == 'dhcp':
        a = [i * a0 for i in [1, 0, 0]]
        b = [i * a0 for i in [-0.5, math.sqrt(3) / 2, 0]]
        c = [i * a0 for i in [0, 0, 2 * math.sqrt(8.0 / 3.0)]]
        u = [1, 0, 0]
        v = [0, 1, 0]
        w = [0, 0, 1]
        sites = [[0.0, 0.0, 0.0], [1.0 / 3.0, 2.0 / 3.0, 0.25],
                 [2.0 / 3.0, 1.0 / 3.0, 0.5], [1.0 / 3.0, 2.0 / 3.0, 0.75]]
    elif structure == 'fcc':  # Primitive
        a = [i * a0 for i in [1, 0, 0]]
        b = [i * a0 for i in [0, 1, 0]]
        c = [i * a0 for i in [0, 0, 1]]
        u = [0.0, 0.5, 0.5]
        v = [0.5, 0.0, 0.5]
        w = [0.5, 0.5, 0.0]
        sites = [[0.0, 0.0, 0.0]]
    else:  # structure == 'sfe'
        a = [i * a0 for i in [math.sqrt(6) / 4, -math.sqrt(2) / 4, 0]]
        b = [i * a0 for i in [math.sqrt(6) / 4, math.sqrt(2) / 4, 0]]
        c = [i * a0 for i in [0, 0, math.sqrt(3)]]
        u = [1, 0, 0]
        v = [0, 1, 0]
        w = [0, 0, 1]
        sites = [[0.0, 0.0, 0.0], [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
                 [2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0]]
    s = ' '.join([str(x) for x in a])
    s += '\n' + ' '.join([str(x) for x in b])
    s += '\n' + ' '.join([str(x) for x in c])
    s += '\n' + ' '.join([str(x) for x in u])
    s += '\n' + ' '.join([str(x) for x in v])
    s += '\n' + ' '.join([str(x) for x in w])
    for site in sites:
        s += '\n' + ' '.join([str(x) for x in site]) + ' ' + concLine
    f = open(dest + '/lat.in', 'w+')
    f.write(s)
    f.close()


def makeSQSCell(structure, a0, nx, ny, nz, dest='.'):
    """ Makes sqscell.out file for given structure """
    if structure == 'fcc':
        a = [i * nx for i in [0.0, 0.5, 0.5]]
        b = [i * ny for i in [0.5, 0.0, 0.5]]
        c = [i * nz for i in [0.5, 0.5, 0.0]]
    else:
        a = [i * nx for i in [1.0, 0.0, 0.0]]
        b = [i * ny for i in [0.0, 1.0, 0.0]]
        c = [i * nz for i in [0.0, 0.0, 1.0]]
    s = '1\n'
    s += '\n' + ' '.join([str(x) for x in a])
    s += '\n' + ' '.join([str(x) for x in b])
    s += '\n' + ' '.join([str(x) for x in c])
    f = open(dest + '/sqscell.out', 'w+')
    f.write(s)
    f.close()


# ============================================================================
#  Script Generation
# ============================================================================


def genSQSScript(ranges, wr, wn, wd, temp, nTotal, nSQS, nice, dest):
    """Make script for running mcsqs"""
    corr = 'nice -n %d mcsqs -l=lat.in' % nice
    for r in range(len(ranges)):
        corr += ' -%d=%f' % (r + 2, ranges[r])
    cmd = 'nice -n %d mcsqs -l=lat.in -n=%d -rc' % (nice, nTotal)
    cmd += ' -wr=%f -wn=%f -wd=%f -T=%f' % (wr, wn, wd, temp)
    s = '#!/bin/bash'
    s += '\n' + corr
    for i in range(nSQS):
        name = '%s_%d' % (dest, i)
        ipcmd = cmd + ' -ip=%d' % i
        s += '\nscreen -dm -S %s %s' % (name, ipcmd)
    outFile = '%s/mcsqs.sh' % dest
    f = open(outFile, 'w+')
    f.write(s)
    f.close()
    sp.call(['chmod', '+x', outFile])


def getCompName(elements, comp):
    """Return string for composition name"""
    compName = ''
    for i in range(len(elements)):
        if comp[i] != 0:
            compName += elements[i]
            if comp[i] != 1:
                compName += str(comp[i])
    return compName


def getConcLine(elements, comp):
    """Return concentrations string for ATAT"""
    concLine = []
    concTotal = float(sum(comp))
    for i in range(len(elements)):
        if comp[i] != 0:
            concLine += ['%s=%f' % (elements[i], comp[i] / concTotal)]
    return ','.join(concLine)

# ============================================================================
#  Main Program
# ============================================================================
# Get user input
elements = raw_input('Elements (e.g. A, B, C): ').split(',')
elements = sorted([e.strip() for e in elements])  # strip whitespace
print 'Elements are ' + ' '.join(elements)
comp = raw_input('Composition (e.g. 1, 2, 3): ').split(',')
comp = map(int, comp)  # convert to int
structure = raw_input('Structure (fcc, hcp, dhcp or sfe): ').strip().lower()
# add arbitary stacking sequence option
if structure[0] == 'f':
    structure = 'fcc'
    nPer = 1
elif structure[0] == 'h':
    structure = 'hcp'
    nPer = 2
elif structure[0] == 'd':
    structure = 'dhcp'
    nPer = 4
else:  # structure == sfe
    nPer = 3
print '%s, %d atoms per unit cell' % (structure, nPer)
periodicities = raw_input('Periodicities (e.g. 1, 1, 2): ').split(',')
periodicities = map(int, periodicities)
a0 = float(raw_input('Lattice parameter (fcc): '))
if 'cp' in structure:
    a0 = a0 / math.sqrt(2)
ranges = raw_input('Cluster ranges (e.g. 4, 4, 2): ').split(',')
ranges = map(int, ranges)
nSQS = int(raw_input('Number of SQSs: '))
nice = int(raw_input('Nice value: '))
wr = raw_input('Perfect match weight: ')
if not wr:
    wr = 1
else:
    wr = float(wr)
wn = raw_input('Cluster point decrease: ')
if not wn:
    wn = 1
else:
    wn = float(wn)
wd = raw_input('Cluster diameter decay: ')
if not wd:
    wd = 0
else:
    wd = float(wd)
temp = raw_input('Temperature: ')
if not temp:
    temp = 1
else:
    temp = float(temp)

# timeout option

# Make name, elements list, concentration line
na = periodicities[0]
nb = periodicities[1]
nc = periodicities[2]

compName = getCompName(elements, comp)
concLine = getConcLine(elements, comp)
nTotal = nPer * na * nb * nc
path = '%s_%s_%d' % (compName, structure, nTotal)

if not os.path.exists(path):
    print 'Creating directory %s...' % path
    os.makedirs(path)
makeLat(structure, a0, concLine, path)
makeSQSCell(structure, a0, na, nb, nc, path)
genSQSScript(ranges, wr, wn, wd, temp, nTotal, nSQS, nice, path)
