#!/usr/bin/env python
# ============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  January 27, 2016
#  Script to make concatenated POTCAR files
# ============================================================================
# Update this path to point to your potentials
PATHTOPBE = '/.../PBE'


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
        elements = sorted(elements.split(','))
        elementSets += [elements]
        print ' '.join(elements)
print

# Make POTCAR files
for elements in elementSets:
    makePOT(elements)
