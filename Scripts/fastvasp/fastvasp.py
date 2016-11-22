# =============================================================================
# Jonas Kaufman jlkaufman@hmc.edu
# June 17, 2016
# Fast and flexible script VASP calculations on different platforms.
# =============================================================================
import os
import subprocess as sp
import collections
import math
import numpy as np

from Cell import *
from Vinputs import *

# =============================================================================
# User Information
# =============================================================================
EMAIL = 'jlkaufman@g.hmc.edu'
ALLOCATION = 'TG-DMR140093'
PBE = '$HOME/PBE'

# =============================================================================
# User Input
# =============================================================================


def getYesOrNo(message, default=True):
    """ Asks y/n question and returns user's answer """
    answer = raw_input(message).lower()
    if answer:
        if answer[0] == 'y':
            return True
        elif answer[0] == 'n':
            return False
        else:
            return default
    else:
        return default


def getInput(message, inType, isList=False, default=None):
    """
    Prompts user for input of a particular type and returns it
    Supports ints, floats, strings or lists of those specified inType and
    isList
    """
    answer = raw_input(message)
    if not answer:
        print default
        return default
    else:
        answer = answer.split(',')
        try:
            if inType[0] == 'i':
                answer = [int(a) for a in answer]
            elif inType[0] == 'f':
                answer = [float(a) for a in answer]
            else:
                answer = [str(a) for a in answer]
        except ValueError:
            print default
            return default
        if isList:
            return answer
        else:
            return answer[0]

# =============================================================================
# Job Submission
# =============================================================================

def makeSLURMScript(job, cores, nodes, time, allocation, commands,
                    email, send=True, priority='normal', dest='./submit.sh'):
    """ Makes submission script for SLURM queueing system """
    # Set SLURM tags
    hrs = time / 60
    mins = time % 60
    s = '#!/bin/bash'
    s += '\n#SBATCH -J %s' % job                  # specify job name
    s += '\n#SBATCH -o %s_' % job + '%j'          # write output to this file
    s += '\n#SBATCH -n %d' % cores                # request cores
    s += '\n#SBATCH -N %d' % nodes                # request nodes
    s += '\n#SBATCH -t %02d:%02d:00' % (hrs, mins) # set maximum wall time
    s += '\n#SBATCH -A %s' % allocation           # specify project
    s += '\n#SBATCH -p %s' % priority             # request queue
    if send:
        s += '\n#SBATCH --mail-user=%s' % email   # set email
        s += '\n#SBATCH --mail-type=all'        # send all emails
    # Specify commands
    for cmd in commands:
        s += '\n%s' % cmd
    # Write file
    f = open(dest, 'w+')
    f.write(s)
    f.close()


def makeBatchScript(scripts):
    """ Makes script to submit multiple SLURM scripts to the queue """
    # navigate to each directory, submit each script, return, etc.
    return


def makeScreenScript(job, commands, timeout=100000000, dest='./screen.sh'):
    """ Makes script to run multiple commands in separate screens """

    cmd = 'nice -n %d timeout %d mcsqs -l=lat.in -n=%d -rc' % (
        NICE, time, nTotal)
    s = '#!/bin/bash'
    newdirs = []
    for d in dirs:
        for i in range(nSQS):
            newdir = '%s_%d' % (d, i)
            if i == 0:
                s += '\nmv %s %s' % (d, newdir)
            else:
                s += '\ncp -r %s_0 %s' % (d, newdir)
            newdirs += [newdir]
    for d in newdirs:
        s += '\nscreen -dm -t %s bash -c \"cd %s ; %s\" ' % (d, d, cmd)
    f = open('mcsqs.sh', 'w+')
    f.write(s)
    f.close()
    sp.call(['chmod', '+x', 'mcsqs.sh'])

# !/bin/bash
# SBATCH -J test_40
# SBATCH -o test_40_%j
# SBATCH -n 16
# SBATCH -N 1
# SBATCH -p normal
# SBATCH -t 02:00:00
# SBATCH --mail-user=jlkaufman@hmc.edu
# SBATCH --mail-type=all
# SBATCH -A TG-DMR140093
# module load vasp
# ibrun vasp_std > vasp_output.out

def readINCAR(filename='INCAR'):
    """ Reads INCAR file and returns dictionary of tags """
    f = open(filename,'r')
    tags = collections.OrderedDict()
    while True:
        nextLine = f.readline()
        if not nextLine:
            break
        elif nextLine != '\n':
            words = nextLine.split()
            key = words[0]
            value = ' '.join(words[2:])
            tags[key] = value
    return tags

def makeINCAR(tags, dest='.', filename='INCAR'):
    """ Takes in ordered dictionary of tags and makes INCAR file in dest """
    s = ''
    for k, v in tags.items():
        s += '%s = %s\n'%(k,str(v))
    f = open('%s/%s'%(dest,filename),'w+')
    f.write(s)
    f.close()


def makePOTCAR(elements, dest='.', potentials=PBE):
    """Make POTCAR file for given elements."""
    name = 'POTCAR_%s'%'-'.join(elements)
    file = '%s/%s'%(dest, name)
    f = open(file, 'w+')
    for e in elements:
        p = open('%s/%s/POTCAR'%(PATHTOPBE, e), 'r')
        f.write(p.read())
        p.close()
    f.close()

def checkForFile(file):
    """ Checks if file exists in particular directory """
    # pathToFile
    if os.path.exists(file):
        return True
    else:
        return False

# def getRadii():
#     radii = {}
#     f=open('radii.txt', 'r')
#     while True:
#         nextLine = f.readline()
#         if not nextLine:
#             break
#         else:
#             nextLine = nextLine.split()
#             element = nextLine[0]
#             radius = nextLine[1]
#             radii[element] = float(radius)
#     return radii
    
# def getLat(composition, structure):

#     radii = getRadii()
#     for e in composition.keys():

#     if structure == 'HCP':
#         a = 2*radius
#     elif structure == 'BCC'
#         a = 4*radius/math.sqrt(3)
#     else: # structure == 'FCC'
#         a = 4*radius/math.sqrt(2)
#     return a

# def getComposition(POSCAR):
#     """Get composition from a POSCAR file"""

# =============================================================================
# Calculation Initialization
# =============================================================================

def initializeCell():
    """Return cell from existing POSCAR or user input."""
    isPOSCAR = False
    if checkForFile('POSCAR'):
        cell = Cell().loadFromPOSCAR('POSCAR')
        print 'Found POSCAR:\n%s'%cell
        if getYesOrNo('Use this POSCAR? (y/n): ', default=True):
            isPOSCAR = True
    else:
        print 'No POSCAR found'

    if isPOSCAR:
        return cell
    else:
        print 'Making new POSCAR'
        cell = Cell()
        
        element = getInput('Element: ', 'str', default='Cu')
        cell.setElements([element])
        structure = getInput('Structure (FCC, BCC, HCP, SFE): ',
            'str', default='FCC').upper()
        cell.setHeader('%s_%s'%(element, structure))
        print structure
        
        if structure == 'SFE':
            nz = getInput('Number of 111 layers (14, 18, 22; 0 for custom): ',
                'int', default=14)
            if nz == 0:
                sequence = getInput('Custom sequence: ', 'str', default='ABC')
            elif nz == 18: 
                sequence = 'ABCABACBACBACBCABC'
            elif nz == 22: 
                sequence = 'ABCABCBACBACBACBABCABC' 
            else: 
                sequence = 'ABCACBACBACABC' # nz == 14 
            layers = list(sequence)

            radius = math.sqrt(2)/4
            a = [2*radius,0,0]
            b = [-radius,math.sqrt(3)*radius,0]
            c = [0,0,nz*math.sqrt(3)/3]
            z = 0.0
            sites = []
            while layers:
                next = layers.pop(0)
                if next == 'A':
                    shift = 0.0
                elif next == 'B':
                    shift = 1.0
                else: # next == 'C'
                    shift = 2.0
                x = shift/3.0
                y = shift*2.0/3.0
                sites += [[x,y,z]]
                z += (1.0/nz)

        elif structure == 'BCC':
            a = [0.5,0.5,-0.5]
            b = [-0.5,0.5,0.5]
            c = [0.5,-0.5,0.5]
            sites = [[0.0,0.0,0.0]]

        elif structure == 'HCP':
            a = [1,0,0]
            b = [-0.5,math.sqrt(3)/2,0]
            c = [0,0,math.sqrt(8.0/3.0)] # ideal c/a ratio
            site = [[0.0,0.0,0.0],[1.0/3.0,2.0/3.0,0.5]]

        else: # structure == 'FCC'
            a = [0.0,0.5,0.5]
            b = [0.5,0.0,0.5]
            c = [0.5,0.5,0.0]
            sites = [[0.0,0.0,0.0]]

        cell.setLatticeVectors([a,b,c])
        for s in sites:
            cell.addSite(0,Site(s)) # add sites
        a0 = getInput('Lattice parameter?: ', 'float', default=3.636)
        cell.setA0(a0)
        return cell
        cell.sendToPOSCAR('POSCAR')

def initializeKPOINTS():
    """Read and or make KPOINTS file."""
    isKPOINTS = False
    if checkForFile('KPOINTS'):
        f = open('KPOINTS', 'r')
        print 'Found KPOINTS file:\n%s'%f.read()
        f.close()
        if getYesOrNo('Use this KPOINTS? (y/n): ', default=True):
            isKPOINTS = True
    else:
        print 'No KPOINTS found'

    if isKPOINTS:
        print ''
    else:
        print 'Making new KPOINTS'

def autoSubdivisions(self, length, a0=0, POSCAR='POSCAR'):
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
        
        return subdivisions


# add emailing for Tera jobs
# def writeEmail(outFile)
#     s = 'Subject: %s'%subject
#     s += '\n'
#     for l in lines
#         s += '\n%s'%l
#     f = open(outFile, 'w+')
#     f.write(s)
#     f.close()


    # sendmail user@example.com  < /tmp/email.txt



# =============================================================================
# Main Program
# =============================================================================

def main():
    """Execute main functionality of the script."""
    # cell = initializeCell()
    # print cell
    # initializeKPOINTS()
    kpts = KPOINTS()
    print kpts



    # POTCAR

    # VASPVERSION = getInput('VASP Version?: ', 'int', 5)

    # elements = []
    # # makePOTCAR(elements)
    # LREAL = '.FALSE.'
    # nAtoms = cell.numberOfAtoms
    # print '%d atoms detected'%nAtoms
    # if nAtoms > 20:
    #     print 'Number of atoms is greater than 20'
    #     autoReal = getYesOrNo('Use LREAL = \'Auto\'? (y/n):', 
    #         default=True)
    #     if autoReal: LREAL = 'Auto'

    """
    detect elements and make potcar 
    (both 4/5)

    suggest incar parameters 
    allow user to set in car tags from terminal 

    magnetic moments 

    database of lattice constants and averaging suggestions 

    vasp calculation class 

    stages 
    in car 
    k points 
    pot car 

    queuing system? 
    daisy chain jobs option

    Would you like to run jobs in series or parallel?
    running in screen...
    """



if __name__ == '__main__':
    main()   




# incorporate spin-polarization (user supplied or default moments)
