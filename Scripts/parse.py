#!/usr/bin/env python
# ============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  June 6, 2016
#  Script to parse VASP output directories different calculation types
# ============================================================================
print 'Importing modules...\n'
import os
import matplotlib
matplotlib.use('Agg')  # Force mpl to not use any Xwindows backend.
import matplotlib.pyplot as plt
from pylab import *  # this includes numpy as np
import scipy.optimize as optimize
DPI = 300                       # image resolution for plots

# ============================================================================
# VASP File Parsing
# ============================================================================


def getNumberofIons(location='.', file='POSCAR'):
    """ parses a VASP 5 style POSCAR file and pulls out the number of ions """
    f = open('%s/%s' % (location, file), 'r')
    for i in range(7):
        nextLine = f.readline()
    number = sum([int(n) for n in nextLine.split()])
    f.close()
    return number


def getSubdivisions(location='.', file='KPOINTS'):
    """ parses a KPOINTS file and pulls out the subdivisions """
    f = open('%s/%s' % (location, file), 'r')
    for i in range(4):
        nextLine = f.readline()
    # Return list of subdivisions or length
    subdivisions = [int(n) for n in nextLine.split()[0:3]]
    f.close()
    return subdivisions


def getIrrKpoints(location='.', file='OUTCAR'):
    """ parses an OUTCAR file and pulls out number of irreducible k-points """
    f = open('%s/%s' % (location, file), 'r')
    kPoints = 0
    while True:
        nextLine = f.readline()
        if not nextLine:
            break
        if 'irreducible' in nextLine:
            kPointsLine = nextLine.split()
            kPoints = int(kPointsLine[1])
    f.close()
    return kPoints


def getTime(location='.', file='OUTCAR'):
    """ parses an OUTCAR file and pulls out the run time """
    f = open('%s/%s' % (location, file), 'r')
    time = 0
    while True:
        nextLine = f.readline()
        if not nextLine:
            break
        if 'Total CPU time used' in nextLine:
            time = float(nextLine.split()[5])
    f.close()
    return time


def getENCUT(location='.', file='OUTCAR'):
    """ parses an OUTCAR file and pulls out the ENCUT value """
    f = open('%s/%s' % (location, file), 'r')
    cutoff = 0
    while True:
        nextLine = f.readline()
        if not nextLine:
            break
        if 'energy-cutoff' in nextLine:
            cutoff = float(nextLine.split()[2])
    f.close()
    return cutoff


def getEnergies(location='.', file='OUTCAR'):
    """ parses an OUTCAR file and pulls out the energy of the system
    after each ionic step """
    f = open('%s/%s' % (location, file), 'r')
    energies = []
    while True:
        nextLine = f.readline()
        if not nextLine:
            if not energies:
                energies = [0]
            break
        # this function returns E(sigma->0), not TOTEN
        if 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)' in nextLine:
            f.readline()    # line of dashes
            f.readline()    # TOTEN line
            f.readline()    # blank line
            energyLine = f.readline().split()
            energies += [float(energyLine[6])]
    f.close()
    return energies


def getPressures(location='.', file='OUTCAR'):
    """ parses an OUTCAR file and pulls out the pressure
    after each ionic step """
    f = open('%s/%s' % (location, file), 'r')
    pressures = []
    while True:
        nextLine = f.readline()
        if not nextLine:
            if not pressures:
                pressures = [0]
            break
        if 'external pressure' in nextLine:
            pressureLine = nextLine.split()
            pressures += [float(pressureLine[3])]
    f.close()
    return pressures


def getPullay(location='.', file='OUTCAR'):
    """ parses an OUTCAR file and pulls out Pullay stress
    after each ionic step """
    f = open('%s/%s' % (location, file), 'r')
    stresses = []
    while True:
        nextLine = f.readline()
        if not nextLine:
            if not stresses:
                stresses = [0]
            break
        if 'external pressure' in nextLine:
            pressureLine = nextLine.split()
            stresses += [float(pressureLine[8])]
    f.close()
    return stresses


def getVolumes(location='.', file='OUTCAR'):
    """ parses an OUTCAR file and pulls out the volume
    after each ionic step """
    f = open('%s/%s' % (location, file), 'r')
    volumes = []
    while True:
        nextLine = f.readline()
        if not nextLine:
            if not volumes:
                volumes = [0]
            break
        if 'VOLUME and BASIS-vectors are now :' in nextLine:
            f.readline()    # dashed line
            f.readline()    # cutoff energy
            volumeLine = f.readline().split()
            volumes += [float(volumeLine[4])]
    f.close()
    return volumes


def getVectors(location='.', file='OUTCAR'):
    """ parses an OUTCAR file and pulls out the lattice
    vectors after each ionic step """
    f = open('%s/%s' % (location, file), 'r')
    vectors = []
    while True:
        nextLine = f.readline()
        if not nextLine:
            if not vectors:
                vectors = [[[0, 0, 0], [0, 0, 0], [0, 0, 0]]]
            break
        if 'VOLUME and BASIS-vectors are now :' in nextLine:
            f.readline()    # dashed line
            f.readline()    # cutoff energy
            f.readline()    # volume line
            f.readline()    # direct lattice vectors
            a = [float(x) for x in f.readline().split()[:3]]
            b = [float(x) for x in f.readline().split()[:3]]
            c = [float(x) for x in f.readline().split()[:3]]
            vectors.append([a, b, c])
    return vectors


def getStatic(location='.', file='OUTCAR'):
    """ parses OUTCAR file for the NSW tag """
    f = open('%s/%s' % (location, file), 'r')
    while True:
        nextLine = f.readline()
        if not nextLine:
            static = False
            break
        if 'NSW' in nextLine:
            nsw = nextLine.split()[2]
            if int(nsw) == 0:
                static = True
            else:
                static = False
            break
    return static

# ============================================================================
#  Analysis, Display and Organization
# ============================================================================


def getVASPDirs():
    """ Get all VASP directories under immediately under current directory"""
    dirs = []
    for d in next(os.walk('.'))[1]:
        # Check if directory contains OUTCAR
        if os.path.isfile('%s/OUTCAR' % d):
            print 'Found OUTCAR in directory %s' % d
            if d not in dirs:
                dirs += [d]
        # Check if directories under it contain OUTCAR
        for root, sdirs, files in os.walk(d, topdown=False):
            for s in sdirs:
                path = os.path.join(root, s)
                if os.path.isfile('%s/OUTCAR' % path):
                    print 'Found OUTCAR in directory %s' % path
                    if d not in dirs:
                        dirs += [d]
    return dirs


def makeTable(headings, data):
    """ makes and returns a table string out of headings and data columns """
    # Get correct widths
    nCol = len(headings)
    widths = []
    for j in range(nCol):
        w = len(str(headings[j]))
        for i in range(len(data[j])):
            w = max(w, len(str(data[j][i])))
        widths += [w]

    # Make table
    table = ''
    for j in range(nCol):
        table += str(headings[j]).ljust(widths[j])
        table += '\t'
    table += '\n'
    for i in range(len(data[0])):
        for j in range(nCol):
            table += str(data[j][i]).ljust(widths[j])
            table += '\t'
        table += '\n'
    return table


def getDeltas(values):
    """ calculate value-to-value changes of a list """
    deltas = [0]
    for i in range(len(values) - 1):
        deltas += [values[i + 1] - values[i]]
    return deltas


def checkConvergence(values, CRIT=0.0001):
    """ checks convergence of a quantity """
    converged = False
    c = 0
    # hard convergence, make sure change from converged value to
    # all subsequent values is < CRIT
    for i in range(len(values)):
        deltas = [(values[i] - v) for v in values[i:]]
        if all([abs(d) < CRIT for d in deltas]):
            if not converged:
                c = i
            converged = True
        else:
            converged = False
    return c

# def displayRun(run):
#   data = []
#   dirName = run[0]
#   EList = run[1]
#   VList = run[2]
#   aList = run[3]
#   PList = run[4]
#   sList = run[5]
#   kPoints = run[6]
#   time = run[7]
#   print dirName
#   headings = ['E0','Volume','|a|','|b|','|c|','Pressure','Pullay stress']
#   nSteps = len(EList)
#   for i in range(nSteps):
#       lats = aList[i]
#       a = lats[0]
#       b = lats[1]
#       c = lats[2]
#       data += [[EList[i],VList[i],a,b,c,PList[i],sList[i]]]
#   printTable(headings,data)
#   print '%d irr k-points'%kPoints
#   print '%d ionic steps'%nSteps # subtract 1?
#   print '%d seconds'%time
#   print

# def displayFinal(runList): # fix the column formatting
#   fins = finalValues(runList)
#   data = []
#   print 'Final values'
#   headings = ['Name','E0','Volume','|a|','|b|','|c|',
#       'Pressure','Pullay stress','Irr k-points','Time']
#   for i in range(len(runList)):
#       d = fins[0][i]
#       E = fins[1][i]
#       V = fins[2][i]
#       lats = fins[3][i]
#       a = lats[0]
#       b = lats[1]
#       c = lats[2]
#       P = fins[4][i]
#       s = fins[5][i]
#       k = fins[6][i]
#       t = fins[7][i]
#       data += [[d,E,V,a,b,c,P,s,k,t]]
#   printTable(headings,data)
#   print

# def finalValues(runList):
#   dList = []
#   EFins = []
#   VFins = []
#   aFins = []
#   PFins = []
#   sFins = []
#   kList = []
#   tList = []
#   for run in runList:
#       dList += [run[0]]
#       EFins += [run[1][-1]]
#       VFins += [run[2][-1]]
#       aFins += [run[3][-1]]
#       PFins += [run[4][-1]]
#       sFins += [run[5][-1]]
#       kList += [run[6]]
#       tList += [run[7]]
#   return (dList,EFins,VFins,aFins,PFins,sFins,kList,tList)

# # plot property over the course of calculation, option to have multiple(R,RR,S)
# # with line in between to indicate

# ============================================================================
#  Parsing for Different Calculation Types
# ============================================================================


def parseENCUTDirectories(dirs):
    """ parses ENCUT directories and output information """
    dirs.sort(key=float)
    # Extract data
    nAtoms = getNumberofIons(dirs[0])
    cutoffs = []
    energies = []
    times = []
    for d in dirs:
        d = str(d)
        cutoffs += [getENCUT(d)]
        energies += [getEnergies(d)[-1]]  # only final energy
        times += [getTime(d)]

    # Calculate energy changes, per atom, check convergence
    energiesPer = [e / nAtoms for e in energies]
    deltas = getDeltas(energiesPer)
    c = checkConvergence(energiesPer)
    if c != 0:
        conv = 'Converged to within 0.1 meV/atom at ENCUT = %d\n' % cutoffs[c]
        conv += 'Difference between last and converged is %f eV/atom' %\
            (energiesPer[-1] - energiesPer[c])
    else:
        conv = 'No converence to within 1 meV/atom'

    # Graph
    plt.plot(cutoffs, energiesPer, '-ro', mec='r')             # Plot data
    if c != 0:
        plt.plot(cutoffs[c], energiesPer[c], 'bx', ms=10)  # Plot converged
    plt.xlabel('Cutoff energy (eV)', fontsize='x-large')
    plt.ylabel('Energy per atom (eV)', fontsize='x-large')
    plt.tight_layout()
    plt.savefig('%s.png' % SYS, dpi=DPI)

    # Make table, print, write datafile
    print '\nResults:'
    headings = ['ENCUT', 'E0', 'E0/atom', 'delta', 'time']
    data = [cutoffs, energies, energiesPer, deltas, times]
    table = makeTable(headings, data)
    text = '%s\n%s\n' % (table, conv)
    print text
    f = open('%s.txt' % SYS, 'w+')
    f.write(text)
    f.close()


def parseKPOINTSDirectories(dirs):
    """ parses KPOINTS directories and output information """
    dirs.sort(key=float)
    # Extract data
    nAtoms = getNumberofIons(dirs[0])
    N1 = []
    N2 = []
    N3 = []
    irr = []
    energies = []
    times = []
    for d in dirs:
        subdivisions = getSubdivisions(d)
        N1 += [subdivisions[0]]
        if len(subdivisions) < 3:
            N2 += [0]
            N3 += [0]
        else:
            N2 += [subdivisions[1]]
            N3 += [subdivisions[2]]
        irr += [getIrrKpoints(d)]
        energies += [getEnergies(d)[-1]]  # only final energy
        times += [getTime(d)]

    # Calculate energy changes, per atom, check convergence
    energiesPer = [e / nAtoms for e in energies]
    deltas = getDeltas(energiesPer)
    c = checkConvergence(energiesPer)
    if c != 0:
        conv = 'Converged to within 0.1 meV/atom at ' + \
            'length %d, subdivisions %d %d %d, %d irr k-points\n'\
            % (int(dirs[c]), N1[c], N2[c], N3[c], irr[c])
        conv += 'Difference between last and converged is %f eV/atom' %\
            (energiesPer[-1] - energiesPer[c])
    else:
        conv = 'No converence to within 1 meV/atom'

    # Graph
    plt.plot(irr, energiesPer, '-ro', mec='r')             # Plot data
    if c != 0:
        plt.plot(irr[c], energiesPer[c], 'bx', ms=10)  # Plot converged
    plt.xlabel('Irreducible k-points', fontsize='x-large')
    plt.ylabel('Energy per atom (eV)', fontsize='x-large')
    plt.tight_layout()
    plt.savefig('%s.png' % SYS, dpi=DPI)

    # Make table, print, write datafile
    print '\nResults:'
    headings = ['l', 'N1', 'N2', 'N3', 'ikp', 'E0', 'E0/atom', 'delta', 'time']
    data = [dirs, N1, N2, N3, irr, energies, energiesPer, deltas, times]
    table = makeTable(headings, data)
    text = '%s\n%s\n' % (table, conv)
    print text
    f = open('%s.txt' % SYS, 'w+')
    f.write(text)
    f.close()


def parseBirchDirectories(dirs, residuals=True):
    """ parses Birch-Murnaghan directories and output information """
    # Figure out directory structure
    rList = [n for n in os.listdir(dirs[0])
             if os.path.isdir(os.path.join(dirs[0], n))]

    # Sort directories by lattice parameter
    dirs.sort(key=float)

    # Extract data
    finalEnergies = []
    volumes = []
    totalTimes = []
    for d in dirs:
        if not rList:
            enrgs = getEnergies(d)
            vol = getVolumes(d)[-1]
            t = getTime(d)
        else:
            enrgs = []
            t = 0
            for r in rList:
                path = '%s/%s' % (d, r)

                enrgs += getEnergies(path)
                vol = getVolumes(path)[-1]
                t += getTime(path)

        finalEnergies += [enrgs[-1]]
        volumes += [vol]
        totalTimes += [t]

    # Perform Birch-Murnaghan Fitting
    birchPars, cov, infodict, mesg, ier = fitBirch(finalEnergies, volumes)
    if ier == 1 or ier == 2 or ier == 3 or ier == 4:
        birch = 'Birch-Murnaghan fit successful'
        E0, B0, BP, V0 = birchPars

        # Compute R^2
        ss_err = (infodict['fvec']**2).sum()
        y = np.array(finalEnergies)
        ss_tot = ((y - y.mean())**2).sum()
        rSquared = 1 - (ss_err / ss_tot)
        birch += ' with R^2 = %f' % rSquared

        # Make Graph
        if residuals:  # Plot Residuals, if requested
            f, (ax1, ax2) = plt.subplots(2, 1, sharex=True,
                                         gridspec_kw={'height_ratios': [1, 3]})
            ax1.plot(volumes, infodict['fvec'], 'bo', mec='b')
            ax1.set_ylabel('Residuals', fontsize='x-large')
        else:
            f, ax2 = plt.subplots(1, 1)
        NPOINTS = 100
        vFit = np.linspace(min(volumes), max(volumes), NPOINTS)
        ax2.plot(vFit, Birch(birchPars, vFit), '-b')          # Plot fit
        ax2.plot(volumes, finalEnergies, 'ro', mec='r')    # Plot data
        ax2.plot(V0, E0, 'bx', ms=10)                  # Plot min
        ax2.set_xlabel('Volume ($\AA^3$)', fontsize='x-large')
        ax2.set_ylabel('Energy (eV)', fontsize='x-large')
        # Save Figure
        f.tight_layout()
        plt.savefig('%s.png' % SYS, dpi=DPI)
        # Write out fit parameters
        birch += '\nE0  = \t%f eV' % E0
        birch += '\nV0  = \t%f A^3' % V0
        birch += '\nB0  = \t%f GPa' % (B0 * 160.21766)
        birch += '\nB0\' = \t%f' % BP
    else:
        birch = 'Birch-Murnaghan fit unsuccessful'

    # Make table, print, write datafile
    print '\nResults:'
    headings = ['a', 'V', 'E0', 'total time']
    data = [dirs, volumes, finalEnergies, totalTimes]
    table = makeTable(headings, data)
    text = '%s\n%s\n' % (table, birch)
    print text
    f = open('%s.txt' % SYS, 'w+')
    f.write(text)
    f.close()


def parseConstDirectories(dirs, residuals=True):
    """ parses elastic constants (C'/C44) directories and output info """
    # Figure out directory structure
    rList = [n for n in os.listdir(dirs[0])
             if os.path.isdir(os.path.join(dirs[0], n))]
    # print rList

    # Sort directories by delta
    dirs.sort(key=float)
    deltas = [float(d) for d in dirs]

    if 0.0 not in deltas:
        # Get E0 from user
        E0 = raw_input('\nInput the ground state energy E0 or press ENTER: ')
        if E0:
            deltas = [0.0] + deltas
            finalEnergies = [float(E0)]
            totalTimes = [0]
            volumes = [0]

        else:
            finalEnergies = []
            totalTimes = []
            volumes = []

    else:
        finalEnergies = []
        totalTimes = []
        volumes = []

    # Extract data
    for d in dirs:
        if not rList:
            enrgs = getEnergies(d)
            vols = getVolumes(d)
            t = getTime(d)
        else:
            enrgs = []
            vols = []
            t = 0
            for r in rList:
                path = '%s/%s' % (d, r)
                enrgs += getEnergies(path)
                vols += getVolumes(path)
                t += getTime(path)

        finalEnergies += [enrgs[-1]]
        volumes += [vols[-1]]
        totalTimes += [t]

    # Perform Fitting
    V0 = volumes[-1]
    pars, cov, infodict, mesg, ier = fitConst(finalEnergies, deltas, V0)
    const = pars[0]
    E0 = pars[1]
    if ier == 1 or ier == 2 or ier == 3 or ier == 4:
        elastic = 'Elastic constant fit successful'

        # Compute R^2
        ss_err = (infodict['fvec']**2).sum()
        y = np.array(finalEnergies)
        ss_tot = ((y - y.mean())**2).sum()
        rSquared = 1 - (ss_err / ss_tot)
        elastic += ' with R^2 = %f' % rSquared

        # Make Graph
        if residuals:  # Plot Re`siduals, if requested
            f, (ax1, ax2) = plt.subplots(2, 1, sharex=True,
                                         gridspec_kw={'height_ratios': [1, 3]})
            ax1.plot(deltas, infodict['fvec'], 'bo', mec='b')
            ax1.set_ylabel('Residuals', fontsize='x-large')
        else:
            f, ax2 = plt.subplots(1, 1)
        NPOINTS = 100
        dFit = np.linspace(0.0, max(deltas), NPOINTS)
        print const
        ax2.plot(dFit, elasticEqn(pars, V0, dFit), '-b')    # Plot fit
        ax2.plot(deltas, finalEnergies, 'ro', mec='r')         # Plot data
        ax2.set_xlabel('Delta', fontsize='x-large')
        ax2.set_ylabel('Energy (eV)', fontsize='x-large')
        # Save Figure
        f.tight_layout()
        plt.savefig('%s.png' % SYS, dpi=DPI)
        # Write out fit parameters
        elastic += '\nC  = \t%f eV/A^3' % const
        elastic += '\nC  = \t%f GPa' % (const * 160.21766)
        elastic += '\nE0 = \t%f eV' % E0

    else:
        elastic = 'Elastic constant fit unsuccessful'

    # Make table, print, write datafile
    print '\nResults:'
    headings = ['delta', 'V', 'E0', 'total time']
    data = [deltas, volumes, finalEnergies, totalTimes]
    table = makeTable(headings, data)
    text = '%s\n%s\n' % (table, elastic)
    print text
    f = open('%s.txt' % SYS, 'w+')
    f.write(text)
    f.close()

# ============================================================================
#  Birch Murnaghan Fitting
# ============================================================================


def Birch(parameters, vol):
    """
    given a vector of parameters and volumes, return a vector of energies.
    equation From Wikipedia
    """
    E0 = parameters[0]
    B0 = parameters[1]
    BP = parameters[2]
    V0 = parameters[3]
    term12 = (((V0 / vol)**(2.0 / 3.0)) - 1.0)
    term3 = (6.0 - (4.0 * ((V0 / vol)**(2.0 / 3.0))))
    E = E0 + (9.0 * V0 * B0 / 16.0) * \
        ((term12**3.0) * BP + (term12**2.0) * term3)
    return E


def objectiveBirch(pars, y, x):
    """ function to be minimized """
    err = y - Birch(pars, x)
    return err


def fitBirch(energies, volumes):
    """ fit energy/volume data to BM EoS """
    # fit a parabola to the data to get guesses
    a, b, c = np.polyfit(volumes, energies, 2)
    V0 = -b / (2 * a)
    E0 = a * V0**2 + b * V0 + c
    B0 = 2 * a * V0
    BP = 4
    x0 = [E0, B0, BP, V0]
    # fit the data to Birch
    return optimize.leastsq(objectiveBirch, x0, (energies, volumes),
                            full_output=True)

# ============================================================================
#  C'/C44 Fitting
# ============================================================================


def elasticEqn(params, V, delta):
    """ E = E0 +2*V*C'*delta^2 """
    const = params[0]
    E0 = params[1]
    E = E0 + 2 * V * const * (delta**2.0)
    return E


def objectiveConst(params, V, y, x):
    """ function to be minimized """
    err = y - elasticEqn(params, V, x)
    return err


def fitConst(energies, deltas, V):
    """ fit energy/delta data to the quadratic from Zaddach 2013 """
    # fit a parabola to the data to get guess with nonfixed intercept
    a, b, c = np.polyfit(energies, deltas, 2)
    c0 = c / (2 * V)
    E0 = a
    x0 = [c0, E0]
    # fit the data
    # E0 = [E0]*len(energies)
    deltas = np.array(deltas)
    return optimize.leastsq(objectiveConst, x0, (V, energies, deltas),
                            full_output=True)

# ============================================================================
#  Main Program
# ============================================================================
# Identify CWD
SYS = os.path.basename(os.getcwd())
print 'Directory: %s\n' % SYS

# Get parsing option
print 'Parsing options:'
print '\t'.join(['Birch-Murn', '- b'])
print '\t'.join(['Elastic const', '- c'])
print '\t'.join(['Encut conv', '- e'])
print '\t'.join(['K-point conv', '- k'])

p = raw_input('What kind of calculations are you parsing?: ')
if p and p[0].lower() in 'kebc':
    p = p[0].lower()
else:
    p = 'b'
print p, '\n'

# Get VASP directories immediately under CWD
dirs = getVASPDirs()
print 'Directories are %s' % ', '.join(dirs)

# Parse directories
if p == 'k':
    parseKPOINTSDirectories(dirs)
elif p == 'e':
    parseENCUTDirectories(dirs)
elif p == 'c':
    parseConstDirectories(dirs)
else:  # p == 'b'
    parseBirchDirectories(dirs)
