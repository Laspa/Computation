#===============================================================================
#  Kate Reed kireed@hmc.edu 
#  June 13, 2017
#  Script to check that all tags have been included for magnetic calculations in 
#  VASP. Also prints parts of POSCAR, POTCAR, and INCAR files so I can check
#  that the elements and their magnetic moments are correctly ordered.
#===============================================================================
def fileExists(fileName):
    try:
        searchfile = open(fileName, "r")
        searchfile.close()

    except IOError as e:
        print "File " + fileName + " not found."

def searchFileForUnfilled(fileName, keyword):
    try:
        searchfile = open(fileName, "r")
        i = 0
        for line in searchfile:
            if keyword in line: 
                print "File " + fileName + ", line " + str(i) + ": " + str(line)
        searchfile.close()

    except IOError as e:
        print "File " + fileName + " not found."

def searchFileForExists(fileName, keyword):
    try:
        searchfile = open(fileName, "r")
        exists = False
        for line in searchfile:
            if keyword in line: 
                exists = True
        searchfile.close()

        if not exists:
            print "keyword " + str(keyword) + " not found in file " + str(fileName) + "."
        else:
            print "keyword " + str(keyword) + " specified."

    except IOError as e:
        print "File " + fileName + " not found."

def checkAtomOrder(POSCARFileName, POTCARFileName, INCARFileName):
    try:
        poscarFile = open(POSCARFileName, "r")
        poscarLines = poscarFile.readlines()
        print poscarLines[5]
        poscarFile.close()

    except IOError as e:
        print "File " + POSCARFileName + " not found."

    try:
        potcarFile = open(POTCARFileName, "r")
        for line in potcarFile:
            if "PAW_PBE" in line and "TITEL" not in line: 
                print line
        potcarFile.close()

    except IOError as e:
        print "File " + POTCARFileName + " not found."

    try:
        incarFile = open(INCARFileName, "r")
        for line in incarFile:
            if "MAGMOM" in line: 
                print line
        incarFile.close()
        
    except IOError as e:
        print "File " + POTCARFileName + " not found."

print "Checking POSCAR and POTCAR orderings"
checkAtomOrder("POSCAR", "POTCAR", "INCAR")

print "Checking for unfilled tags---------------------------------------"
searchFileForUnfilled("INCAR", "FILLTAG")
searchFileForUnfilled("KPOINTS", "FILLTAG")

print "Checking that magnetic settings are specified--------------------"
searchFileForExists("INCAR", "ISPIN")
searchFileForExists("INCAR", "MAGMOM")
searchFileForExists("INCAR", "LORBIT")







