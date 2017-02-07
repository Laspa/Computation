#!/usr/bin/env python
# ============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  May 30, 2016
#  Script to run set of fixed volume calculations for on Stampede
# ============================================================================
"""
Add POSCAR, POTCAR and KPOINTS files to the working directory.
Make sure to change the INCAR parameters to your liking.
ENCUT, PREC, etc. are set in the Main Program section
"""
from Cell import *
import os
import subprocess as sp
import math
import numpy as np
EMAIL = '...@hmc.edu'  # Change this to your own!
ALLOCATION = 'TG-DMR140093'

# ============================================================================
# POSCAR Building
# ============================================================================


def makePOSCAR(header, a0, dest):
    """ Creates updated POSCAR with new lattice parameter """
    cell = Cell().loadFromPOSCAR('POSCAR')
    cell.setHeader(header)
    cell.setA0(a0)
    cell.sendToPOSCAR('%s/POSCAR' % dest)

# ============================================================================
# VASP Files
# ============================================================================


def makeINCAR(system, static, dest, index, isif=2):
    """ Make INCAR file """
    s = 'SYSTEM = %s' % system
    s += '\nPREC = %s' % PREC
    s += '\nENCUT = %d' % ENCUT
    s += '\nISTART = 0'
    s += '\nICHARG = 2'
    if static:
        s += '\nNSW = 0'
        s += '\nISMEAR = -5'
    else:
        s += '\nNSW = 50'
        s += '\nISIF = %d' % isif
        s += '\nIBRION = 1'
        s += '\nISMEAR = 1'
        s += '\nSIGMA = 0.2'
    s += '\nNPAR = %d' % NPAR
    f = open('%s/INCAR%d' % (dest, index), 'w+')
    f.write(s)
    f.close()

# ============================================================================
# File Preparation
# ============================================================================


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


def prepareDirectories(system, aList, relax):
    """ Prepare directories for each value """
    directories = []
    scripts = []

    # Make subdirectories
    for a in aList:
        aDir = '%.5f' % a
        header = '%s_%.5f' % (system, a)
        if not os.path.exists(aDir):
            os.makedirs(aDir)

        if relax:
            makeINCAR(header, False, aDir, 0)
            makeINCAR(header, True, aDir, 1)
            nSteps = 2
        else:
            makeINCAR(header, True, aDir, 0)
            nSteps = 1

        for f in ['KPOINTS', 'POTCAR']:
            sp.call(['cp', f, '%s/%s' % (aDir, f)])

        # Make POSCAR in starting directory
        makePOSCAR(header, a, aDir)

        # Make submission script
        genSubScript(header, nSteps, aDir)
        directories += [aDir]
        scripts += ['%s_submit' % header]

    # Create batch script in parent directory
    genBatchScript(system, directories, scripts)
    sp.call(['chmod', '+x', '%s_batch' % system])

# ============================================================================
#  Main Program
# ============================================================================
print 'Add POSCAR, POTCAR and KPOINTS files to the working directory\n'

# VASP Settings
ENCUT = 450
PREC = 'Accurate'
NPAR = 4

# Get user inputs
system = raw_input('System name: ')
print system, '\n'

guess = float(raw_input('Equilibrium lattice parameter guess (A): '))
print guess, '\n'
aRange = float(raw_input('Lattice parameter range (A): '))
print aRange, '\n'
nVal = int(raw_input('Number of values: '))
print nVal, '\n'

relax = raw_input('Relax ions? (y/n): ')
if relax and relax[0].upper() == 'Y':
    relax = True
else:
    relax = False
print relax, '\n'

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

# NPER = raw_input('Number of nodes per 16 cores (1 or 2, default 1): ')
# if NPER and int(NPER) == 2:
#     NNODES = NCORES * int(NPER)
# else:
#     NPER = 1
#     NNODES = NCORES / 16
# print NPER, '\n'
NPER = 1

RUNTIME = int(raw_input('Single job run time in minutes: '))
print RUNTIME, '\n'

aMin = guess - 0.5 * aRange
aMax = guess + 0.5 * aRange
aList = np.linspace(aMin, aMax, nVal).tolist()
print 'a values: %s\n' % ', '.join([str(a) for a in aList])

# Prepare directories
print 'Creating directories...\n'
prepareDirectories(system, aList, relax)
