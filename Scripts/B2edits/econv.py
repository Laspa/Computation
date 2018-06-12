#!/usr/bin/env python
#===============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu 
#  June 3, 2016    
#  Script to run energy cutoff convergence calculations on Stampede
#===============================================================================
""" Add POSCAR, POTCAR, KPOINTS and INCAR files to the working directory. """
import subprocess as sp
import numpy as np
import collections
EMAIL = 'eyhwang@hmc.edu' # Change this to your own!
ALLOCATION = 'TG-DMR140093'

#===============================================================================
# INCAR Manipulation
#===============================================================================
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
	#s += '\nmodule load vasp'                  	# load vasp module
	s += '\nibrun vasp_std > vasp.output'	    # run vasp
	
	f = open('%s/%s_submit'%(dest,name),'w+')
	f.write(s)
	f.close()

def prepareENCUTDirectories(system, energies):
	""" Prepares directory for each length value """
	directories = []
	scripts = []
	tags = readINCAR()								# Read INCAR file

	# Go through energies
	for e in energies:
		dest = str(e)
		name = '%s_%s'%(system,e) 
		sp.call(['mkdir',dest])						# Make directory
		directories += [dest]
		tags['ENCUT'] = e
		makeINCAR(tags, dest)						# Make INCAR file
		for file in ['POSCAR','POTCAR','KPOINTS']: 	# Copy VASP files
			sp.call(['cp',file,dest])	
		genSubScriptSimple(name,dest) 				# Make sub script	
		scripts += ['%s_submit'%name]
	
	# Make batch scripts
        genBatchScript(system,directories,scripts)

def genBatchScript(name,directories,scripts,dest='.'):
	""" 
	Creates script name_batch in dest to enter all
	directories and submit scripts to the Stampede queue
	"""
	s = '#!/bin/bash'
        s += '\nmkdir enconv'
        for i in range(len(directories)):
                s += '\nmv %s enconv' % directories[i]
	for i in range(len(directories)):
		s += '\ncd enconv/%s'%directories[i]		# enter working directory
		s += '\nsbatch %s'%scripts[i]		# submit job
		s += '\ncd ../..' 						# exit working directory
	s += '\nmv %s_batch enconv' % name
        filename = '%s/%s_batch'%(dest,name)
	f = open(filename,'w+')
	f.write(s)								# write file
	f.close()
	sp.call(['chmod','+x',filename])		# make executable

#===============================================================================
#  Main Program
#===============================================================================
print 'Add POSCAR, POTCAR, KPOINTS and INCAR files to the working directory\n'

# Get user inputs
system = raw_input('System name: ')
print system, '\n'

minVal = int(raw_input('Minimum ENCUT: '))
print minVal, '\n'

maxVal = int(raw_input('Maximum ENCUT: '))
print maxVal, '\n'

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

NNODES = int(raw_input('Number of nodes? (default 1): '))
if NCORES / 64 > NNODES:
    NCORES = NNODES * 64
#NPER = raw_input('Number of nodes per 16 cores (1 or 2, default 1): ')
#if NPER and int(NPER) == 2:
#	NNODES = NCORES*int(NPER)
#else:
#	NPER = 1
#	NNODES = NCORES/16
print NNODES, '\n'

RUNTIME = int(raw_input('Single job run time in minutes: '))
print RUNTIME, '\n'

# Make list of energies
energies = [int(e) for e in np.linspace(minVal, maxVal, nVal)]
print 'ENCUT values: %s\n'%', '.join([str(e) for e in energies])

# Create directories and batch script
print 'Making directories...\n'
prepareENCUTDirectories(system, energies)
