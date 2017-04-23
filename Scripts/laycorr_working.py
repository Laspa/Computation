# =============================================================================
# Jonas Kaufman jlkaufman@hmc.edu
# June 28, 2016
# Script to analyze correlations between layers in SFE SQSs.
# =============================================================================
import subprocess as sp
import os.path

# =============================================================================
# ATAT Structure Manipulation
# =============================================================================

def getZPositions(inFile='lat.in'):
    """Read ATAT structure and return a list of distinct z positions."""
    f = open(inFile)                    # Read input file
    lines = f.readlines()
    f.close()
    zMax = float(lines[5].split()[2])   # Get maximum z
    zList = []
    for L in lines[6:]:                 # Get z values
        zValue = float(L.split()[2])
        if zValue >= zMax:              # Keep values under maximum
            zValue += -zMax
        if zValue not in zList:
            zList += [zValue]
    zList.sort()
    return zList

def selectLayers(layers=[], inFile='lat.in', outFile='llat.in'):
    """
    Read ATAT structure and write structure containing only specified layers 
    """
    CRIT = 0.000001                     # Criterion for equalitiy
    zList = getZPositions(inFile)       # Read z positions
    outString = ''
    f = open(inFile, 'r')               # Read input file
    lines = f.readlines()
    f.close()
    zMax = float(lines[5].split()[2])   # Get maximum z
    for L in lines[:6]:                 # Write coordinates, vectors
        outString += L
    for L in lines[6:]:                 # Select positions, write them
        zValue = float(L.split()[2])
        if zValue >= zMax:              # Keep values under maximum
            zValue += -zMax
        for i in layers:
            z = zList[i]
            if abs(zValue - z) < CRIT:
                outString += L
    f = open(outFile, 'w+')             # Write output file
    f.write(outString)
    f.close()

def randomStructure(inFile='lbestsqs.out', latFile='llat.in', outFile='lrnd.out' ):
    """ """
    f = open(latFile)                    # Read lat file
    lines = f.readlines()
    f.close()
    compo = lines[6].split()[3]          # Get composition
    compo = compo.split(',')
    elements = []
    concs = []
    for c in compo:
        c = c.split('=')
        elements += [c[0]]
        concs += [float(c[1])]
    f = open(inFile)                    # Read structure
    lines = f.readlines()
    f.close()
    nTotal = len(lines) - 6             # Get number of atoms
    atoms = []                          # Make list of atoms
    for i in range(len(elements)):
        n = int(round(concs[i]*float(nTotal)))
        atoms += [elements[i]]*n
    outString = ''       
    for L in lines[:6]:                 
        outString += L                 
    for i in range(nTotal):             # Add sites with new composition     
        L = lines[i+6]
        a = atoms[i]
        L = L.split()
        L = ' '.join(L[:3] + [a]) + '\n'
        outString += L
    f = open(outFile, 'w+')             # Write output file
    f.write(outString)
    f.close

# =============================================================================
# ATAT Commands
# =============================================================================

def getCorr(distances, struct='bestsqs.out', lat='lat.in', rnd=False):
    """Run corrdump on structure and return correlations"""
    if rnd:
        outFile = 'rndcorr.out'
    else:
        outFile = 'strcorr.out'
    command = 'corrdump -ro -crf=greg'
    command += ' -s=%s -l=%s'%(struct, lat)
    i = 2
    for dist in distances:
        if dist != 0:
            command += ' -%d=%f'%(i, dist)
        i += 1
    if rnd:
        command += ' -rnd'
    command += ' > %s'%outFile
    print command
    sp.call(command, shell=True)
    f = open(outFile, 'r')
    correlations = f.read().split()
    f.close()
    correlations = map(float, correlations)
    return correlations

def getClus():
    """Run getclus and return cluster nubmers, lengths, and multiplicities."""
    outFile = 'getclus.out'
    command = 'getclus > %s'%outFile
    print command
    sp.call(command, shell=True)
    numbers = []
    lengths = []
    multips = []
    f = open(outFile, 'r')
    clusters = f.readlines()
    f.close()
    for clus in clusters:
        clus = clus.split()
        numbers += [int(clus[0])]
        lengths += [float(clus[1])]
        multips += [int(clus[2])] 
    return (numbers, lengths, multips)

# =============================================================================
# ATAT Correlation Analysis
# =============================================================================

def correlationsReport(layers, distances, struct='bestsqs.out', lat='lat.in'):
    """Analyze correlations of specified layers and print a report."""
    if os.path.exists('clusters.out'):
        isClust = True
        move = 'mv clusters.out oclusters.out'  # Save original clusters 
        sp.call(move.split())
    else:
        isClust = False
    llat = 'l' + lat
    lstr = 'l' + struct
    # Select layers, get correlations
    selectLayers(layers, inFile=lat, outFile=llat)
    selectLayers(layers, inFile=struct, outFile=lstr)
    randomStructure(inFile=lstr, latFile=llat, outFile='lrnd.out')
    strCorr = getCorr(distances, lstr, llat, rnd=False)
    rndCorr = getCorr(distances, 'lrnd.out', llat, rnd=True)
    misCorr = []                         # Calculate mismatches
    for i in range(len(strCorr)):
        misCorr += [strCorr[i] - rndCorr[i]]
    numbers, lengths, multips = getClus()
    # absSum = sum(map(abs, misCorr))
    if isClust:
        move = 'mv oclusters.out clusters.out'  # Move clusters back
        sp.call(move.split())
    outString = ''
    outString += '\nnumber\tlength\tmultip\tstructure\trandom\tmismatch'                                  # Construct output
    for i in range(len(numbers)):
        # outString += '\n%d\t%.6f\t%.6f\t%.6f\t%.6f'%(numbers[i],
        #                                             lengths[i],
        #                                             strCorr[i],
        #                                             rndCorr[i],
        #                                             misCorr[i])
        outString += '\n%d\t%f\t%d\t%f\t%f\t%f'%(numbers[i],
                                                    lengths[i],
                                                    multips[i],
                                                    strCorr[i],
                                                    rndCorr[i], 
                                                    misCorr[i])
    # outString += '\nAbsolute mismatch sum = %f'%absSum
    outString += '\n'
    print outString                                 # Print, write output
    outFile = '_'.join(map(str, layers)) + '.out'
    f = open(outFile, 'w+')
    f.write(outString)
    f.close()
    # return absSum
   

# =============================================================================
# Main Program
# =============================================================================

def main():
    """Execute main functionality of the script."""

    inFile = raw_input("sqs file name: ")

    #Read layers, report results
    zList = getZPositions(inFile)
    print 'Found %d layers:'%len(zList)
    i = 0
    print 'L\tz'
    for z in zList:
        print '%d\t%.6f'%(i, z)
        i += 1
    print

    distances = [4, 0, 0]

    layersList = [[8,0]]

    for n in range(8):
        layersList += [[n,n+1]]  

    mismatches = []
    for layers in layersList:
        print 'Isolating layers ' + ' '.join(map(str, layers))
        correlationsReport(layers, distances, inFile)

    # summary table

if __name__ == '__main__':
    main()   
