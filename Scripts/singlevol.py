#!/usr/bin/env python
#=========================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  August 17, 2016
#  Script to run a single fixed volume calculation on Stampede
#=========================================================================
"""
Add POSCAR, POTCAR and KPOINTS files to the working directory.
Make sure to change the INCAR parameters to your liking.
"""
from Cell import *
import os
import subprocess as sp
import math
import numpy as np
EMAIL = 'jlkaufman@hmc.edu'  # Change this to your own!
ALLOCATION = 'TG-DMR140093'

#=========================================================================
# POSCAR Building
#=========================================================================


def makePOSCAR(header, a0, dest):
    """ Creates updated POSCAR with new lattice parameter """
    cell = Cell().loadFromPOSCAR('POSCAR')
    cell.setHeader(header)
    cell.setA0(a0)
    cell.sendToPOSCAR('%s/POSCAR' % dest)

#=========================================================================
# VASP Files
#=========================================================================

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
    f = open('%s/INCAR%d'%(dest, index),'w+')
    f.write(s)
    f.close()


#===============================================================================
# KPOINTS Generation
#===============================================================================
def makeKPOINTS(subdivisions, header='', dest='.', gamma=True):
    """ Make KPOINTS in dest file with specified subdivisions, center """
    print 'Making KPOINTS file...'
    if gamma:
        center = 'Gamma'
    else:
        center = 'Monkhorst'
    header = str(header)
    
    s = 'Automatic mesh %s'%header          # header
    s += '\n0'                              # 0 -> automatic generation scheme 
    s += '\n%s'%center                      # grid center
    s += '\n%d %d %d'%tuple(subdivisions)   # subdivisions along recip. vect.
    s += '\n0 0 0'                          # optional shift 
    
    f = open('%s/KPOINTS'%dest,'w+')
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
    a1,a2,a3 = cell.latticeVectors
    b1 = np.cross(a2,a3)/(np.dot(a1,np.cross(a2,a3)))/a0
    b2 = np.cross(a3,a1)/(np.dot(a2,np.cross(a3,a1)))/a0
    b3 = np.cross(a1,a2)/(np.dot(a3,np.cross(a1,a2)))/a0

    bNorms = [np.linalg.norm(b) for b in [b1,b2,b3]]

    print bNorms
    # Calculate subdivision as per
    # http://cms.mpi.univie.ac.at/vasp/vasp/Automatic_k_mesh_generation.html
    subdivisions = [1]*3
    for i in [0,1,2]:
        subdivisions[i] = int(max(1,((length*bNorms[i])+0.5)))
    KPPRA = int(np.prod(subdivisions)*nAtoms) # k-points per recip. atom
    
    print 'Subdivisions are %d %d %d'%tuple(subdivisions)
    print 'KPPRA = %d'%KPPRA
    
    return subdivisions

#=========================================================================
# File Preparation
#=========================================================================


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

#=========================================================================
#  Main Program
#=========================================================================
print 'Add POSCAR, POTCAR and KPOINTS files to the working directory\n'

# VASP Settings
ENCUT = 450
PREC = 'Accurate'
NPAR = 4
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

# Get user inputs
system = raw_input('System name: ')
print system, '\n'

a0 = float(raw_input('Lattice parameter (A): '))
print a0, '\n'

klength = int(raw_input('Length for automatic k-mesh: '))
if klength > 0:
    subdivisions = autoSubdivisions(klength, a0=a0)
    makeKPOINTS(subdivisions,gamma=True)   

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
print 'Creating files...\n'
# Make subdirectories
aDir = '.'
header = system
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
makePOSCAR(header, a0, aDir)

# Make submission script
genSubScript(header, nSteps, aDir)
