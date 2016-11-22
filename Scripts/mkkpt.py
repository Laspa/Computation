#!/usr/bin/env python
#===============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu 
#  June 9, 2016
#  Script to make KPOINTS file from POSCAR automatically
#===============================================================================
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

#===============================================================================
#  Main Program
#===============================================================================
center = raw_input('Center (Gamma or Monkhorst): ')
gamma = True
if center and center[0].upper() == 'M':
	center = 'Monkhorst'
	gamma = False
else:
	center = 'Gamma'
print center, '\n'

length = int(raw_input('K length: '))
print length, '\n'

# Make file
subdivisions = autoSubdivisions(length)
makeKPOINTS(subdivisions,gamma=gamma)	
print
