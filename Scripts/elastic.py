#!/usr/bin/env python
# ============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  June 1, 2016
#  IN PROGRESS Script to run elastic constants calculations on Stampede
#  Feel free to modify this script if you need to run these calculations
# ============================================================================
"""
Add POSCAR, POTCAR and KPOINTS files to the working directory.
Make sure to change the INCAR parameters to your liking.
"""
from Cell import *
import os
import subprocess as sp
import math
import numpy as np

EMAIL = '...@hmc.edu'  # email address for SLURM
ALLOCATION = 'TG-DMR140093'  # allocation on Stampede

# ============================================================================
# POSCAR Building
# ============================================================================


def makePOSCAR(header, a0, strainTensor, dest):
    """ Creates FCC POSCAR under the given strain tensor """
    cell = Cell().loadFromPOSCAR('POSCAR')
    cell.setHeader(header)
    cell.setA0(a0)
    # Update each lattice vector
    e1, e2, e3, e4, e5, e6 = strainTensor
    for i in range(len(cell.latticeVectors)):
        x, y, z = cell.latticeVectors[i]
        newVector = [(1 + e1) * x + (e6 * y) / 2 + (e5 * z) / 2,
                     (e6 * x) / 2 + (1 + e2) * y + (e4 * z) / 2,
                     (e5 * x) / 2 + (e4 * y) / 2 + (1 + e3) * z]
        cell.latticeVectors[i] = newVector
    # Send to file
    cell.sendToPOSCAR('%s/POSCAR' % dest)

# ============================================================================
# VASP Files
# ============================================================================


def makeStaticINCAR(system, continuation, dest):
    """ Makes INCAR file for static calculation """
    s = 'SYSTEM = %s' % system
    s += '\nPREC = %s' % PREC
    s += '\nENCUT = %d' % ENCUT
    s += '\nNSW = 0'
    s += '\nISMEAR = -5'
    s += '\nSIGMA = 0.2'
    if continuation:
        s += '\nISTART = 1'
        s += '\nICHARG = 0'
    else:
        s += '\nISTART = 0'
        s += '\nICHARG = 2'
    s += '\nNPAR = %d' % NPAR
    s += '\nLREAL = %s' % LREAL
    f = open('%s/INCAR' % dest, 'w+')
    f.write(s)
    f.close()


def makeRelaxINCAR(system, continuation, isif, dest):
    """ Makes INCAR file for relaxation calculation """
    s = 'SYSTEM = %s' % system
    s += '\nPREC = %s' % PREC
    s += '\nENCUT = %d' % ENCUT
    s += '\nNSW = 50'
    s += '\nISIF = %d' % isif
    s += '\nIBRION = 2'
    s += '\nISMEAR = 1'
    s += '\nSIGMA = 0.2'
    if continuation:
        s += '\nISTART = 1'
        s += '\nICHARG = 0'
    else:
        s += '\nISTART = 0'
        s += '\nICHARG = 2'
    s += '\nNPAR = %d' % NPAR
    s += '\nLREAL = %s' % LREAL
    f = open('%s/INCAR' % dest, 'w+')
    f.write(s)
    f.close()

# ============================================================================
# File Preparation
# ============================================================================


def genSubScript(name, rList, dest):
    """ Creates a submission script for Stampede's SLURM queueing system """
    # Set SLURM tags
    hrs = RUNTIME / 60
    mins = RUNTIME % 60
    s = '#!/bin/bash'
    s += '\n#SBATCH -J %s' % name      			# specify job name
    s += '\n#SBATCH -o %s_' % name + '%j' 			# write output to this file
    s += '\n#SBATCH -n %d' % NCORES 			  	# request cores
    s += '\n#SBATCH -N %d' % NNODES 			  	# request nodes
    s += '\n#SBATCH -p normal'                	# send to normal queue
    s += '\n#SBATCH -t %02d:%02d:00' % (hrs, mins)  # set maximum wall time
    s += '\n#SBATCH --mail-user=%s' % EMAIL    		# set email
    s += '\n#SBATCH --mail-type=all'            # send all emails
    s += '\n#SBATCH -A %s' % ALLOCATION       		# specify project
    s += '\nmodule load vasp'                  	# load vasp module

    # Multiple relaxation steps
    if len(rList) > 1:
        s += '\ncd %s' % rList[0]							# enter first directory
        s += '\nibrun vasp_std > vasp_output.out'		# run vasp
        for r in rList[1:]:
            s += '\ncp CONTCAR ../%s/POSCAR' % r			# copy CONTCAR
            s += '\nmv WAVECAR ../%s/WAVECAR' % r 		# copy WAVECAR
            s += '\ncd ../%s' % r 						# change directories
            s += '\nibrun vasp_std > vasp_output.out'  # run vasp

    # Single calculation
    else:
        s += '\nibrun vasp_std > vasp_output.out'		# run vasp

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


def prepareDirectories(system, a0, Cqq, tList, rList, isif=2):
    """ Prepare directories for each value """
    directories = []
    scripts = []

    # Make parent directory
    if not os.path.exists(Cqq):
        os.makedirs(Cqq)

    # Make subdirectories
    for t in tList:
        tDir = '%.5f' % t
        path = '%s/%s' % (Cqq, tDir)
        header = '%s_%.5f' % (system, t)
        if not os.path.exists(path):
            os.makedirs(path)

        # Make relaxation directories, if necessary
        continuation = False
        if len(rList) > 1:
            for r in rList:
                rDir = '%s/%s' % (path, r)
                if not continuation:
                    startDir = rDir
                if not os.path.exists(rDir):
                    os.makedirs(rDir)

                # Add INCAR, KPOINTS, POTCAR
                if 'R' in r:
                    makeRelaxINCAR(header, continuation, isif, rDir)
                else:
                    makeStaticINCAR(header, continuation, rDir)
                continuation = True
                for f in ['KPOINTS', 'POTCAR']:
                    sp.call(['cp', f, '%s/%s' % (rDir, f)])

        # Add INCAR, KPOINTS, POTCAR
        else:
            startDir = path
            makeStaticINCAR(header, False, path)
            for f in ['KPOINTS', 'POTCAR']:
                sp.call(['cp', f, '%s/%s' % (path, f)])

        # Choose appropriate strain tensor
        if Cqq == 'Cpr':
            # volume-conserving orthorhombic deformation
            strainTensor = [t, -t, (1 - t**2)**(-1) - 1, 0, 0, 0]
        else:  # Cqq == 'C44'
            # volume-conserving monoclinic deformation
            strainTensor = [0, 0, (1 - t**2)**(-1) - 1, 0, 0, 2 * t]

        # Make POSCAR in starting directory
        makePOSCAR(header, a0, strainTensor, startDir)

        # Make submission script
        genSubScript(header, rList, path)
        directories += [tDir]
        scripts += ['%s_submit' % header]

    # Create batch script in parent directory
    genBatchScript(system, directories, scripts, Cqq)
    sp.call(['chmod', '+x', '%s_batch' % system])

# ============================================================================
#  Main Program
# ============================================================================
print '\nAdd POSCAR, POTCAR and KPOINTS files to working directory.\n'

# Get user inputs
system = raw_input('System name: ')
print system, '\n'

# VASP Settings
ENCUT = 450
PREC = 'Accurate'
NPAR = 4
cell = Cell().loadFromPOSCAR('POSCAR')
nAtoms = sum(cell.elementCounts)
print '%d atoms detected' % nAtoms
if nAtoms > 20:
    LREAL = 'Auto'
else:
    LREAL = '.FALSE.'
print 'Using LREAL = %s' % LREAL, '\n'

# Get strains from user
perc = raw_input('Maximum percent strain (default 5): ')
if perc:
    perc = float(perc)
else:
    perc = 5.0
print str(perc) + ' %', '\n'
tmax = perc / 100.0
nt = raw_input('Number of strains to run (default 5): ')
if nt:
    nt = int(nt)
else:
    nt = 5
print nt, '\n'
tList = np.linspace(perc / nt / 100.0, tmax, nt).tolist()
print 'Strains: ' + ', '.join([str(t) for t in tList])
print 'Zero strain will be added to C44 directories\n'

# Prepare files
a0 = float(raw_input('FCC lattice parameter: '))
print a0, '\n'

# Relaxation settings
rList = []
relax = raw_input('Relax? (y/n): ')
if relax and relax[0].upper() == 'Y':
    relax = True
    rList += ['R']
else:
    relax = False
print relax, '\n'
isif = 2
rList += ['S']

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

RUNTIME = int(raw_input('Single job run time in minutes: '))
print RUNTIME, '\n'

# Prepare directories
print 'Creating directories...\n'
prepareDirectories(system, a0, 'C44', [0.0] + tList, rList, isif)
prepareDirectories(system, a0, 'Cpr', tList, rList, isif)
