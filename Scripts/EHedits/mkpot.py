#!/usr/bin/env python
#===============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu 
#  June 9, 2016
#  Script to make KPOINTS file from POSCAR automatically
# ============================================================================
# Update this path to point to your potentials
PATHTOPBE = '/work/05408/tg847075/stampede2/potentials'


def makePOT(elements):
    """ generate POTCAR file for given elements """
    name = 'POTCAR_%s' % '-'.join(elements)
    print 'Making %s...\n' % name
    # Create blank POTCAR
    f = open(name, 'w+')
    # Write new POTCAR
    for e in elements:
        p = open('%s/%s/POTCAR' % (PATHTOPBE, e), 'r')
        f.write(p.read())
    f.close()
    

def makeCatPOT(elementSets):
    """ generate a concatenated POTCAR file for given elements """
    allelements =''
    for x in elementSets:
	x = x[0]
        allelements += x
    name = 'POTCAR_%s' % allelements 
    print 'Making %s...\n' % name
    f = open(name, 'w+')
    for e in elementSets:
        e = e[0]
        p = open('%s/%s/POTCAR' % (PATHTOPBE, e), 'r')
        f.write(p.read())
    f.close()

# ============================================================================
#  Main Program
# ============================================================================
# Get elements
elementSets = []
print 'Input sets of elements for POTCAR files. Press ENTER when done.'
while True:
    elements = raw_input('Element set (e.g. X,Y,...): ')
    if not elements:
        break
    else:
        elements = elements.split(',')
        elementSets += [elements]
        print ' '.join(elements)
print

# Make POTCAR files
for elements in elementSets:
    makePOT(elements)
#makeCatPOT(elementSets)
