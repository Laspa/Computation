#!/usr/bin/env python
#===============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  January 27, 2016
#  Script to make concatenated PBE POTCAR files
#===============================================================================
# PATHTOPBE = '/home1/03324/tg826232/PBE'
PATHTOPBE = '/Users/Jonas/Google Drive/LASPA/Potentials/PBE'

def makePOT(elements):
	""" generate POTCAR file for given elements """
	name = 'POTCAR_%s'%'-'.join(elements)
	print 'Making %s...\n'%name
	# Create blank POTCAR
	f = open(name,'w+')
	# Write new POTCAR
	for e in elements:
		p = open('%s/%s/POTCAR'%(PATHTOPBE,e),'r')
		f.write(p.read())
	f.close()

#===============================================================================
#  Main Program
#===============================================================================
# Get elements
elementSets = []
print 'Input sets of elements to make POTCAR files for. Press ENTER when done.'
while True:
	elements = raw_input('Element set (e.g. X,Y,...): ')
	if not elements: break
	else:
		elements = elements.split(',')
		elements.sort()
		elementSets += [elements]
		print ' '.join(elements)
print

# Make POTCAR files
for elements in elementSets:
	makePOT(elements)