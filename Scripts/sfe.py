#!/usr/bin/env python
# ========================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  June 22, 2015
#  Script to calculate generalized stacking fault energy
# ========================================================================

from Cell import *
import math
import numpy as np
import os
import subprocess as sp
EMAIL = '...@hmc.edu'  # Change this to your own!
ALLOCATION = 'TG-DMR140093'


def makeSFECell(a0, element, nLayers=9, POSCAR='POSCAR'):
    """"""
    na = 1
    nb = 1
    cell = Cell()
    cell.setHeader(element)
    cell.setA0(a0)
    a = [na * math.sqrt(6) / 4.0, -na * math.sqrt(2) / 4.0, 0]
    b = [nb * math.sqrt(6) / 4.0, nb * math.sqrt(2) / 4.0, 0]
    c = [0, 0, nLayers * math.sqrt(3) / 3.0]
    cell.setLatticeVectors([a, b, c])
    cell.elements = [element]
    cell.setCoordinateSystem('Direct')
    for k in range(nLayers):
        z = float(k) / nLayers
        for j in range(nb):
            y = ((k % 3) / 3.0 + float(j)) / nb
            for i in range(na):
                x = ((k % 3) / 3.0 + float(i)) / na
                cell.addSite(0, Site([x, y, z]))
    cell.sendToPOSCAR(POSCAR)


def makeShiftedPOSCAR(inFile, disp, layerChoice, destDir):
    """ """
    na, nb, nc = getPeriodicities(inFile)
    cell = Cell().loadFromPOSCAR(inFile)

    # shift all positions by layer choice
    for i in range(cell.numberOfElements()):
        for j in range(cell.numberOfAtomsOfElement(i)):
            x, y, z = cell.sites[i][j].position
            z += layerChoice / float(nc)
            if z > 1.0:
                z += -1.0
            cell.sites[i][j].move([x, y, z])

    inCutoff = (nc / 2 - 0.5)
    exCutoff = inCutoff + 1.0

    a, b, c = cell.latticeVectors
    adjdisp = disp / 3.0
    c = np.array(c) + adjdisp * (np.array(a) / na + np.array(b) / nb)
    c = c.tolist()
    cell.setLatticeVectors([a, b, c])

    for i in range(cell.numberOfElements()):
        for j in range(cell.numberOfAtomsOfElement(i)):
            x, y, z = cell.sites[i][j].position
            if z > 1.0:
                z += -1.0
            if z < 0.0:
                z += 1.0
            # Correct for lattice vector
            x += -adjdisp * z / na
            y += -adjdisp * z / nb

            # shift for SFE (ISF vs ESF, general case?)
            # check z value
            if z > inCutoff / nc:
                if disp > 1.0:
                    if z > exCutoff / nc:
                        x += adjdisp / na
                        y += adjdisp / nb
                    else:
                        x += (1.0 / 3.0) / na
                        y += (1.0 / 3.0) / nb
                else:
                    x += adjdisp / na
                    y += adjdisp / nb

            cell.sites[i][j].move([x, y, z])
    #         cell.sites[i][j].setFree(False, False, False)
    #         if z > (inCutoff - 1)/nc:
    #             if disp > 1.0:
    #                 if z < (exCutoff + 1)/nc:
    #                     cell.sites[i][j].setFree(True, True, True)
    #             else:
    #                 if z < (inCutoff + 1)/nc:
    #                     cell.sites[i][j].setFree(True, True, True)

    #

    # cell.SelectiveDynamics = True
    cell.setSiteMobilities(False, False, True)  # allow only z relaxations
    cell.sendToPOSCAR('%s/POSCAR' % destDir)


def getPeriodicities(POSCAR='POSCAR'):
    """Read periodicities of SFE cell"""
    cell = Cell().loadFromPOSCAR(POSCAR)
    a, b, c = cell.latticeVectors
    na = np.linalg.norm(a) / (math.sqrt(2) / 2.0)
    nb = np.linalg.norm(b) / (math.sqrt(2) / 2.0)
    nc = np.linalg.norm(c) / (math.sqrt(3) / 3.0)
    na = int(round(na))
    nb = int(round(nb))
    nc = int(round(nc))
    return na, nb, nc

# ========================================================================
# VASP Files
# ========================================================================


def makeINCAR(system, static, dest, index, isif=2):
    s = 'SYSTEM = %s' % system
    s += '\nPREC = %s' % PREC
    s += '\nENCUT = %d' % ENCUT
    s += '\nISTART = 0'
    s += '\nICHARG = 2'
    if static:
        s += '\nNSW = 0'
        s += '\nISMEAR = -5'
        s += '\nLREAL = .FALSE.'
    else:
        s += '\nNSW = 50'
        s += '\nISIF = %d' % isif
        s += '\nIBRION = 1'
        s += '\nADDGRID = .TRUE.'
        s += '\nEDIFFG = -0.02'
        s += '\nISMEAR = 1'
        s += '\nSIGMA = 0.2'
        s += '\nLCHARG = .FALSE.'
        s += '\nLREAL = %s' % LREAL
    s += '\nNPAR = %d' % NPAR
    s += '\nALGO = %s' % ALGO
    s += '\nLWAVE = .FALSE.'
    f = open('%s/INCAR%d' % (dest, index), 'w+')
    f.write(s)
    f.close()


# ========================================================================
# KPOINTS Generation
# ========================================================================
def makeKPOINTS(subdivisions, header='', dest='.', gamma=True):
    """ Make KPOINTS in dest file with specified subdivisions, center """
    print 'Making KPOINTS file...'
    if gamma:
        center = 'Gamma'
    else:
        center = 'Monkhorst'
    header = str(header)

    s = 'Automatic mesh %s' % header          # header
    s += '\n0'                              # 0 -> automatic generation scheme
    s += '\n%s' % center                      # grid center
    s += '\n%d %d %d' % tuple(subdivisions)   # subdivisions along recip. vect.
    s += '\n0 0 0'                          # optional shift

    f = open('%s/KPOINTS' % dest, 'w+')
    f.write(s)
    f.close()


def autoSubdivisions(length, a0=0, POSCAR='POSCAR'):
    """ Calculate subdivisions automatically from POSCAR """
    # Load POSCAR
    cell = Cell().loadFromPOSCAR(POSCAR)
    if a0 == 0:
        a0 = cell.a0
    nAtoms = sum(cell.elementCounts)

    # Calculate reciprocal lattice vectors
    a1, a2, a3 = cell.latticeVectors
    b1 = np.cross(a2, a3) / (np.dot(a1, np.cross(a2, a3))) / a0
    b2 = np.cross(a3, a1) / (np.dot(a2, np.cross(a3, a1))) / a0
    b3 = np.cross(a1, a2) / (np.dot(a3, np.cross(a1, a2))) / a0

    bNorms = [np.linalg.norm(b) for b in [b1, b2, b3]]

    print bNorms
    # Calculate subdivision as per
    # http://cms.mpi.univie.ac.at/vasp/vasp/Automatic_k_mesh_generation.html
    subdivisions = [1] * 3
    for i in [0, 1, 2]:
        subdivisions[i] = int(max(1, ((length * bNorms[i]) + 0.5)))
    KPPRA = int(np.prod(subdivisions) * nAtoms)  # k-points per recip. atom

    print 'Subdivisions are %d %d %d' % tuple(subdivisions)
    print 'KPPRA = %d' % KPPRA

    return subdivisions

# ========================================================================
# File Preparation
# ========================================================================


def genSubScript(name, nSteps, dest):
    """ Creates a submission script for Stampede's SLURM queueing system """
    # Set SLURM tags
    hrs = RUNTIME / 60
    mins = RUNTIME % 60
    s = '#!/bin/bash'
    s += '\n#SBATCH -J %s' % name               # specify job name
    s += '\n#SBATCH -o %s_' % name + '%j'           # write output to this file
    s += '\n#SBATCH -n %d' % NCORES                 # request cores
    s += '\n#SBATCH -N %d' % NNODES                 # request nodes
    s += '\n#SBATCH -p normal'                  # send to normal queue
    s += '\n#SBATCH -t %02d:%02d:00' % (hrs, mins)  # set maximum wall time
    s += '\n#SBATCH --mail-user=%s' % EMAIL     # set email
    s += '\n#SBATCH --mail-type=all'            # send all emails
    s += '\n#SBATCH -A %s' % ALLOCATION         # specify project
    s += '\ncp INCAR0 INCAR'
    s += '\nmodule load vasp'                   # load vasp module
    s += '\nibrun vasp_std > vasp_output'
    for i in range(nSteps - 1):
        for file in ['POSCAR', 'vasp_output', 'OSZICAR', 'OUTCAR']:
            s += '\nmv %s %s%d' % (file, file, i)
        s += '\ncp INCAR%d INCAR' % (i + 1)
        s += '\ncp CONTCAR POSCAR'
        s += '\nibrun vasp_std > vasp_output'
    # Write file
    f = open('%s/%s_submit' % (dest, name), 'w+')
    f.write(s)
    f.close()


def genBatchScript(name, directories, scripts, dest='.'):
    """ Creates script to quickly send all jobs to the queue """
    s = '#!/bin/bash'
    for i in range(len(directories)):
        s += '\ncd %s' % directories[i]  # enter working directory
        s += '\nsbatch %s' % scripts[i]  # submit job
        s += '\ncd ..'  # exit working directory
    f = open('%s/%s_batch' % (dest, name), 'w+')
    f.write(s)
    f.close()


def prepareDirectories(system, dList, allLayers):
    na, nb, nc = getPeriodicities()
    if allLayers:
        choices = range(nc)
    else:
        choices = [0]

    for chc in choices:
        nDir = str(chc)
        """ Prepare directories for each value """
        directories = []
        scripts = []
        for disp in dList:
            dDir = nDir + '/%.5f' % disp
            header = '%s_%d_%.5f' % (system, chc, disp)
            if not os.path.exists(dDir):
                os.makedirs(dDir)
            makeINCAR(header, False, dDir, 0)
            makeINCAR(header, True, dDir, 1)
            nSteps = 2
            for f in ['KPOINTS', 'POTCAR']:
                sp.call(['cp', f, '%s/%s' % (dDir, f)])

            makeShiftedPOSCAR('POSCAR', disp, chc, dDir)

            # Make submission script
            genSubScript(header, nSteps, dDir)
            directories += ['%.5f' % disp]
            scripts += ['%s_submit' % header]
        # Create batch script in parent directory
        header = header = '%s_%d' % (system, chc)
        genBatchScript(header, directories, scripts, str(chc))
        sp.call(['chmod', '+x', '%d/%s_batch' % (chc, header)])

# ========================================================================
#  Main Program
# ========================================================================
print 'Add POSCAR and POTCAR files to the working directory\n'

# VASP Settings
ENCUT = 450
PREC = 'Accurate'

system = raw_input('System: ')
print system, '\n'

make = raw_input('Make POSCAR? (y/n): ')
if make:
    if make.upper()[0] == 'Y':
        make = True
    else:
        make = False
else:
    make = False

print make

if make:
    nLayers = int(raw_input('Number of layers in cell: '))
    a0 = float(raw_input('FCC lattice parameter (A): '))
    element = raw_input('Element: ')
    makeSFECell(a0, element, nLayers)

klength = int(raw_input('Length for automatic k-mesh: '))
if klength > 0:
    subdivisions = autoSubdivisions(klength)
    makeKPOINTS(subdivisions, gamma=True)

cell = Cell().loadFromPOSCAR('POSCAR')
nAtoms = sum(cell.elementCounts)
print '%d atoms detected' % nAtoms
if nAtoms >= 20:
    LREAL = 'Auto'
    ALGO = 'Fast'
else:
    LREAL = '.FALSE.'
    ALGO = 'Normal'
print 'Using LREAL = %s' % LREAL, '\n'
print 'Using ALGO = %s' % ALGO, '\n'

dList = [0.5, 1.0, 1.5, 2.0]

allLayers = raw_input('Shift at every layer choice? (y/n): ')
if allLayers:
    if allLayers.upper()[0] == 'Y':
        allLayers = True
    else:
        allLayers = False
else:
    allLayers = False

valid = False
while not valid:
    NCORES = raw_input('Number of cores per job (default 16): ')
    if not NCORES:
        NCORES = 16
        valid = True
    elif int(NCORES) % 16 == 0:
        NCORES = int(NCORES)
        valid = True
    else:
        print 'Number of cores must be a multiple of 16'
print NCORES, '\n'

NPER = raw_input('Number of nodes per 16 cores (1 or 2, default 1): ')
if NPER and int(NPER) == 2:
    NNODES = NCORES * int(NPER)
else:
    NPER = 1
    NNODES = NCORES / 16
print NPER, '\n'

if NCORES == 64:
    NPAR = 8
else:
    NPAR = 4

RUNTIME = int(raw_input('Single job run time in minutes: '))
print RUNTIME, '\n'


# Prepare directories
print 'Creating directories...\n'
prepareDirectories(system, dList, allLayers)
