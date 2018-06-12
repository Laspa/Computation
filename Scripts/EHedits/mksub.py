#!/usr/bin/env python
# =================================
#   Emily Hwang eyhwang@hmc.edu
#   May 25, 2018
#   Script to generate a submit file
# =================================

from Cell import *
EMAIL = 'eyhwang@hmc.edu'
ALLOCATION = 'TG-DMR140093'

# ================================
#  File Preparation
# ================================

def genSubmit(name, dest='.'):
    """ Creates a SLURM submission script in dest """
    hrs = RUNTIME / 60
    mins = RUNTIME % 60
    s = '#!/bin/bash'
    s += '\n#SBATCH -J %s' % name
    s += '\n#SBATCH -o %s' % name + '%j'
    s += '\n#SBATCH -n %d' %NCORES
    s += '\n#SBATCH -N %d' %NNODES
    s += '\n#SBATCH -p normal'
    s += '\n#SBATCH -t %02d:%02d:00' % (hrs,mins)
    s += '\n#SBATCH --mail-user=%s' % EMAIL
    s += '\n#SBATCH --mail-type=all'
    s += '\n#SBATCH -A %s' % ALLOCATION
    #s += '\nmodule load vasp'
    s += '\nibrun vasp_std > vasp.output'
    
    f = open('%s/%s_submit' % (dest, name), 'w+')
    f.write(s)
    f.close()
# ==============================
#   Main Program
# ==============================
system = raw_input('System name: ')
print system, '\n'

valid = False
while not valid:
    NCORES = raw_input('Number of cores (default 16): ')
    if not NCORES:
        NCORES = 16
        valid = True
    elif int(NCORES) % 16 == 0:
        NCORES = int(NCORES)
        valid = True
    else:
        print 'Number of cores must be a multiple of 16'
print NCORES, '\n'

correct = False
while not correct:
    NNODES = int(raw_input('Number of nodes? (default 1): '))
    if not NNODES:
        NNODES = 1
        correct = True
    elif NCORES/64 > NNODES:
        NCORES = NNODES * 64
        correct = True
    elif NNODES <= 8:
        NNODES =int(NNODES)
        correct = True
    else:
        print 'Number of nodes must be less than 8'

print 'Number of nodes: ' + str(NNODES)
print 'Number of cores: ' + str(NCORES), '\n'

RUNTIME = int(raw_input('Single job runtime in minutes: '))
print RUNTIME, '\n'

print 'Creating a submission file...\n'
genSubmit(system)
