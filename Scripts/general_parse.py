#!/usr/bin/env python

#===============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  June 29, 2015
#  Script to parse VASP output directories located in a given results directory
#===============================================================================
import os
import fnmatch

# image resolution for plots
DPI = 300
#===============================================================================
#  Parsing
#===============================================================================
def getTime(file):
    """ parses an OUTCAR file and pulls out the run time """
    f = open(file,'r')
    time = 0
    while True:
        nextLine = f.readline()
        if not nextLine:
            break
        if 'Total CPU time used' in nextLine:
            time = float(nextLine.split()[5])
    return time

def getEnergies(file):  # this function returns E(sigma->0), not TOTEN
    """ parses an OUTCAR file and pulls out the energy of the system
    after each ionic step """
    energies = []
    f = open(file,'r')
    while True:
        nextLine = f.readline()
        if not nextLine:
            if not energies: energies = [0]
            break
        if 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)' in nextLine:
            f.readline()    # line of dashes
            f.readline()    # TOTEN line
            f.readline()    # blank line
            energyLine = f.readline().split()
            energies += [float(energyLine[6])]
    return energies
    
def getPressures(file):
    """ parses an OUTCAR file and pulls out the pressure and Pullay stress 
    after each ionic step """
    f = open(file,'r')
    pressures = []
    stresses = []
    while True:
        nextLine = f.readline()
        if not nextLine:
            if not pressures: pressures = [0]
            if not stresses: stresses = [0]
            break
        if 'external pressure' in nextLine:
            pressureLine = nextLine.split()
            pressures += [float(pressureLine[3])]
            stresses += [float(pressureLine[8])]
    return (pressures,stresses)

def getKpoints(file):
    """ parses an OUTCAR file and pulls out number of irreducible k-points """
    f = open(file,'r')
    kPoints = 0
    while True:
        nextLine = f.readline()
        if not nextLine:
            break
        if 'irreducible' in nextLine:
            kPointsLine = nextLine.split()
            kPoints = int(kPointsLine[1])
    return kPoints

def getSizes(file):
    """ parses an OUTCAR file and pulls out the volume and lattice
    parameters after each ionic step """
    f = open(file,'r')
    volumes = []
    vectors = []
    while True:
        nextLine = f.readline()
        if not nextLine:
            if not volumes: volumes = [0]
            if not vectors: vectors = [[0,0,0]]
            break
        if 'VOLUME and BASIS-vectors are now :' in nextLine:
            f.readline()    # dashed line
            f.readline()    # cutoff energy
            volumeLine = f.readline().split()
            volumes += [float(volumeLine[4])]
            for i in range(6):
                f.readline()    # text
            aLine = f.readline().split()
            a = float(aLine[0])
            b = float(aLine[1])
            c = float(aLine[2])
            vectors.append([a,b,c])
    return (volumes, vectors)    

def getStatic(file):
    """ parses OUTCAR file for the NSW tag """
    f = open(file,'r')
    while True:
        nextLine = f.readline()
        if not nextLine:
            static = False
            break
        if 'NSW' in nextLine:
            nsw = nextLine.split()[2]
            if int(nsw) == 0: static = True
            else: static = False 
            break
    return static

def parseResults(directories):
    """ parses each subdirectory of the given results directory """
    runs = []   # list of runs to be returned
    dList = []  # list of directory names
    ELists = [] # list of energy lists
    VLists = [] # list of volume lists
    aLists = [] # list of lattice parameters lists
    PLists = [] # list of pressure lists
    sLists = [] # list of Pullay stress lists
    kList = []  # list of k-points
    tList = []  # list of runtimes
    stList = [] # list of static info
    for d in directories:
        for file in os.listdir('%s/%s'%(CWD,d)):   
            pathToFile = '%s/%s/%s'%(CWD,d,file)    
            if fnmatch.fnmatch(file,'OUTCAR'):
                print 'Found OUTCAR in %s\nParsing...'%d
                dList += [d] # add directory name
                ELists += [getEnergies(pathToFile)] # add energy list
                sizeData = getSizes(pathToFile) 
                VLists += [sizeData[0]] # add volume list
                aLists += [sizeData[1]] # add lattice parameters lists 
                pressureData = getPressures(pathToFile)
                PLists += [pressureData[0]] # add pressure lists
                sLists += [pressureData[1]] # add Pullay stress lists
                kList += [getKpoints(pathToFile)] # add k-points
                tList += [getTime(pathToFile)] # add time
                stList += [getStatic(pathToFile)] # add static
                print 'OUTCAR read\n'
    for i in range(len(dList)):
        runs.append([dList[i],ELists[i],VLists[i],aLists[i],PLists[i],
            sLists[i],kList[i],tList[i],stList[i]])
    return runs

#===============================================================================
#  Display and Organization
#===============================================================================
def printTable(headings, data):
    data.insert(0, headings)
    nCol = len(headings)
    colWidths = []
    for i in range(nCol):
        width = max(len(str(row[i])) for row in data)
        colWidths += [width]
    for row in data:
        string = ''
        for i in range(nCol):
            string += str(row[i]).ljust(colWidths[i])
            string += '\t'
        print string

def displayRun(run):
    data = []
    dirName = run[0]
    EList = run[1]
    VList = run[2]
    aList = run[3]
    PList = run[4]
    sList = run[5]
    kPoints = run[6]
    time = run[7]
    print dirName
    headings = ['E0','Volume','|a|','|b|','|c|','Pressure','Pullay stress']
    nSteps = len(EList)
    for i in range(nSteps):
        lats = aList[i]
        a = lats[0]
        b = lats[1]
        c = lats[2]
        data += [[EList[i],VList[i],a,b,c,PList[i],sList[i]]]
    printTable(headings,data)
    print '%d irr k-points'%kPoints
    print '%d ionic steps'%nSteps # subtract 1?
    print '%d seconds'%time
    print

def displayFinal(runList): # fix the column formatting
    fins = finalValues(runList)
    data = []
    print 'Final values'
    headings = ['Name','E0','Volume','|a|','|b|','|c|',
        'Pressure','Pullay stress','Irr k-points','Time']
    for i in range(len(runList)):
        d = fins[0][i]
        E = fins[1][i]
        V = fins[2][i]
        lats = fins[3][i]
        a = lats[0] 
        b = lats[1]
        c = lats[2]
        P = fins[4][i]
        s = fins[5][i]
        k = fins[6][i]
        t = fins[7][i]
        data += [[d,E,V,a,b,c,P,s,k,t]]
    printTable(headings,data)
    print

def finalValues(runList):
    dList = []
    EFins = []
    VFins = []
    aFins = []
    PFins = []
    sFins = []
    kList = []
    tList = []
    for run in runList:
        dList += [run[0]]
        EFins += [run[1][-1]]
        VFins += [run[2][-1]]
        aFins += [run[3][-1]]
        PFins += [run[4][-1]]
        sFins += [run[5][-1]]
        kList += [run[6]]
        tList += [run[7]]
    return (dList,EFins,VFins,aFins,PFins,sFins,kList,tList)

#===============================================================================
#  Birch Murnaghan Fitting
#===============================================================================
def Birch(parameters,vol):
    """
    given a vector of parameters and volumes, return a vector of energies.
    equation From Wikipedia
    """
    E0 = parameters[0]
    B0 = parameters[1]
    BP = parameters[2]
    V0 = parameters[3]
    term12 = (((V0/vol)**(2.0/3.0)) - 1.0)
    term3 = (6.0 - (4.0*((V0/vol)**(2.0/3.0))))
    E = E0 + (9.0*V0*B0/16.0)*((term12**3.0)*BP + (term12**2.0)*term3)
    return E

def objective(pars,y,x):
    """ function to be minimized """   
    err =  y - Birch(pars,x)
    return err

def fitBirch(EList,VList,jobName):
    """ fit energy/volume data to BM EoS, plot results """
    NPOINTS = 100
    VFit = np.linspace(min(VList),max(VList),NPOINTS)
    # fit a parabola to the data to get guesses
    a,b,c = polyfit(VList,EList,2)
    V0 = -b/(2*a)
    E0 = a*V0**2 + b*V0 + c
    B0 = 2*a*V0
    BP = 4
    x0 = [E0, B0, BP, V0]
    # fit the data to Birch
    print 'Fitting data...'
    BirchPars, ier = optimize.leastsq(objective, x0, args=(EList,VList))
    print 'Done\n'
    # make a plot of the data and the fit
    plot(VFit, Birch(BirchPars,VFit),c='b')
    plot(VList,EList,'ro')
    xlabel('Volume ($\AA^3$)',fontsize='x-large')
    ylabel('Energy (eV)',fontsize='x-large')
    #ax = gca()
    savefig('%s_birch.png'%jobName,dpi=DPI)
    E0,B0,BP,V0 = BirchPars
    print 'Minimum energy:\t%f'%E0
    print 'Minimum volume:\t%f'%V0
    print 'Bulk modulus:\t%f'%(B0*160.2177)

#===============================================================================
#  Main Program
#===============================================================================
# list directories, get parent directory
CWD = os.getcwd()
print 'Found these directories:'
dirs = next(os.walk('.'))[1]
print '\n'.join(dirs), '\n'
valid = False
while not valid:
    parent = raw_input('Directory to parse: ')
    for d in dirs:
        if parent in d:
            parent = d
            print parent, '\n'
            valid = True

# parse the chosen results directory
descendants = [x[0] for x in os.walk(parent)] 
runList = parseResults(descendants)
runList.sort(key=lambda x: x[0]) # sort job list by directory name

# check if runs are static and display data
staticList = []
relaxList = []
for run in runList:
    if not run[8]:
        displayRun(run) 
        relaxList += [run]
    else:
        displayRun(run) 
        staticList += [run]
print 'All runs:'
displayFinal(runList)
if relaxList:
    print 'Relaxation runs:'
    displayFinal(relaxList)
if staticList:
    print 'Static runs:'
    displayFinal(staticList)

    # run fitting on data
    fitting = raw_input('Birch fitting on static runs? (y/n): ')

    if 'y' in fitting or 'Y' in fitting:
        print 'Importing modules...'
        import matplotlib
        matplotlib.use('Agg') # Force mpl to not use any Xwindows backend.
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from pylab import * # this includes numpy as np
        import scipy.optimize as optimize
        print 'Done\n'

        fins = finalValues(staticList)
        energies = fins[1]
        volumes = fins[2]
        fitBirch(energies,volumes,parent)
