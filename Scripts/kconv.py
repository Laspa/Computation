#!/usr/bin/env python
#===============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu 
#  January 25, 2016
#  Script to run k-point convergence calculations on Stampede
#===============================================================================
""" 
Add POSCAR, POTCAR and INCAR files to the working directory. 
The script generates k-points based on a single length quantity which is scaled 
by the reciprocal lattice vectors to produce subdivisions
"""
import subprocess as sp
import numpy as np
from Cell import *
EMAIL = 'jlkaufman@hmc.edu' # Change this to your own!
ALLOCATION = 'TG-DMR140093'

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
	
	s = 'Automatic mesh %s'%header			# header
	s += '\n0' 								# 0 -> automatic generation scheme 
	s += '\n%s'%center 						# grid center
	s += '\n%d %d %d'%tuple(subdivisions)	# subdivisions along recip. vect.
	s += '\n0 0 0' 							# optional shift 
	
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
	
	# Calculate subdivision as per
	# http://cms.mpi.univie.ac.at/vasp/vasp/Automatic_k_mesh_generation.html
	subdivisions = [1]*3
	for i in [0,1,2]:
		subdivisions[i] = int(max(1,((length*bNorms[i])+0.5)))
	KPPRA = int(np.prod(subdivisions)*nAtoms) # k-points per recip. atom
	
	print 'Subdivisions are %d %d %d'%tuple(subdivisions)
	print 'KPPRA = %d'%KPPRA
	
	return subdivisions


#===============================================================================
# File Preparation
#===============================================================================
def genSubScriptSimple(name, dest='.'):
	""" Creates a SLURM submission script name_submit in dest """
	hrs = RUNTIME/60
	mins = RUNTIME%60
	s = '#!/bin/bash'
	s += '\n#SBATCH -J %s'%name      			# specify job name
	s += '\n#SBATCH -o %s_'%name+'%j' 			# write output to this file
	s += '\n#SBATCH -n %d'%NCORES 			  	# request cores
	s += '\n#SBATCH -N %d'%NNODES 			  	# request nodes
	s += '\n#SBATCH -p normal'                	# send to normal queue
	s += '\n#SBATCH -t %02d:%02d:00'%(hrs,mins) # set maximum wall time
	s += '\n#SBATCH --mail-user=%s'%EMAIL    	# set email
	s += '\n#SBATCH --mail-type=all'            # send all emails
	s += '\n#SBATCH -A %s'%ALLOCATION       	# specify project
	s += '\nmodule load vasp'                  	# load vasp module
	s += '\nibrun vasp_std > vasp_output.out'	# run vasp
	
	f = open('%s/%s_submit'%(dest,name),'w+')
	f.write(s)
	f.close()

def prepareKPOINTSDirectories(system, lengths, center):
	""" Prepares directory for each length value """
	directories = []
	scripts = []
	subList = []
	lList = []

	gamma = True
	if center[0] == 'M': gamma = False
	# Go through lengths, calculate subdivisions
	for l in lengths:
		print 'Length %s'%l
		dest = str(l)
		name = '%s_%s'%(system,l)
		subdivisions = autoSubdivisions(l)

		# Check if subdiv. already used
		if not subdivisions in subList:	
			lList += [l]	
			subList += [subdivisions]
			sp.call(['mkdir',dest])						# Make directory		
			directories += [dest]
			makeKPOINTS(subdivisions,l,dest,gamma)		# Make KPOINTS
			for file in ['POSCAR','POTCAR','INCAR']: 	# Copy VASP files
				sp.call(['cp',file,dest])	
			genSubScriptSimple(name,dest) 				# Make sub script	
			scripts += ['%s_submit'%name]
		else:
			print 'Length %d produces duplicate subdivisons, ignoring'%l
		print

	# Print summary table
	table = '\t'.join(['l','N1','N2','N3'])+'\n'
	for i in range(len(lList)):
		table += '\t'.join([str(lList[i])]+
			[str(n) for n in subList[i]])
		table += '\n'
	print table	

	# Make batch scripts
	genBatchScript(system,directories,scripts)

def genBatchScript(name,directories,scripts,dest='.'):
	""" 
	Creates script name_batch in dest to enter all
	directories and submit scripts to the Stampede queue
	"""
	s = '#!/bin/bash'
	for i in range(len(directories)):
		s += '\ncd %s'%directories[i]		# enter working directory
		s += '\nsbatch %s'%scripts[i]		# submit job
		s += '\ncd ..' 						# exit working directory
	
	filename = '%s/%s_batch'%(dest,name)
	f = open(filename,'w+')
	f.write(s)								# write file
	f.close()
	sp.call(['chmod','+x',filename])		# make executable

#===============================================================================
#  Main Program
#===============================================================================
print 'Add POSCAR, POTCAR and INCAR files to the working directory\n'

# Get user inputs
system = raw_input('System name: ')
print system, '\n'

center = raw_input('Center (Gamma or Monkhorst): ')
if center and center[0].upper() == 'M':
	center = 'Monkhorst'
else:
	center = 'Gamma'
print center, '\n'

minLength = int(raw_input('Minimum k length: '))
print minLength, '\n'

maxLength = int(raw_input('Maximum k length: '))
print maxLength, '\n'

nVal = int(raw_input('Number of values: '))
print nVal, '\n'

valid = False
while not valid:
	NCORES = raw_input('Number of cores per job (default 16): ')
	if not NCORES:
		NCORES = 16
		valid = True
	elif int(NCORES)%16 == 0:
		NCORES = int(NCORES)
		valid = True
	else:
		print 'Number of cores must be a multiple of 16'
print NCORES, '\n'

NPER = raw_input('Number of nodes per 16 cores (1 or 2, default 1): ')
if NPER and int(NPER) == 2:
	NNODES = NCORES*int(NPER)
else:
	NPER = 1
	NNODES = NCORES/16
print NPER, '\n'

RUNTIME = int(raw_input('Single job run time in minutes: '))
print RUNTIME, '\n'

# Make list of lengths
lengths = [int(l) for l in np.linspace(minLength, maxLength, nVal)]
print 'Length values: %s\n'%', '.join([str(l) for l in lengths])

# Create directories and batch script
print 'Making directories...\n'
prepareKPOINTSDirectories(system, lengths, center)