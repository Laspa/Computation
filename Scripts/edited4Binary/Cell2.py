# ============================================================================
# Josh Sanz jsanz@hmc.edu
# July 25, 2014
# edited by Jonas Kaufman jlkaufman@hmc.edu
# July 23, 2015
# Classes and methods for manipulating POSCARs and the cells they describe. :)
# ============================================================================
import numpy as np


class Site:
    """Hold a single atomic site for the Cell class."""

    def __init__(self, position=[0, 0, 0], index=0, xfree=True,
                 yfree=True, zfree=True):
        self.position = position    # should be a three element list
        self.index = index
        self.xfree = xfree
        self.yfree = yfree
        self.zfree = zfree

    def __repr__(self):
        s = '%d: ' % self.index
        for x in self.position:
            s += '%f ' % x
        if not (self.xfree and self.yfree and self.zfree):
            if self.xfree:
                s += 'T '
            else:
                s += 'F '
            if self.yfree:
                s += 'T '
            else:
                s += 'F '
            if self.zfree:
                s += 'T '
            else:
                s += 'F '
        return s

    def toString(self):
        s = '%d: ' % self.index
        for x in self.position:
            s += '%f ' % x
        return s[:-1]

    def toStringSelectiveDynamics(self):
        s = '%d: ' % self.index
        for x in self.position:
            s += '%f ' % x
        if self.xfree:
            s += 'T '
        else:
            s += 'F '
        if self.yfree:
            s += 'T '
        else:
            s += 'F '
        if self.zfree:
            s += 'T '
        else:
            s += 'F '
        return s[:-1]

    def move(self, newPos):
        self.position = newPos

    def setFree(self, xFree, yFree, zFree):
        """Set the mobility of site."""
        self.xfree = xFree
        self.yfree = yFree
        self.zfree = zFree

    def copySite(self, otherSite):
        self.position = otherSite.position
        self.index = otherSite.index
        self.xfree = otherSite.xfree
        self.yfree = otherSite.yfree
        self.zfree = otherSite.zfree

    def equals(self, otherSite):
        if (self.index == otherSite.index and
                self.position == otherSite.position):
            return True
        return False


class Cell:
    """Hold lattice vectors and atomic sites of a simulation cell."""

    def __init__(self):
        self.CorD = 'Direct'    # direct coordinates by default
        self.latticeVectors = []
        self.sites = []
        self.elements = []
        self.elementCounts = []
        self.a0 = 1.0
        self.header = ''
        self.SelectiveDynamics = False

    def __repr__(self):
        string = self.header + '\n'
        string += str(self.a0) + '\n'
        for j in self.latticeVectors:
            for x in j:
                string += str(x) + ' '
            string += '\n'
        string += ''.join(self.elements) + '\n'
        string += ' '.join(map(str, self.elementCounts)) + '\n'
        if self.SelectiveDynamics:
            string += 'Selective Dynamics\n'
        string += self.CorD + '\n'
        for el in self.sites:
            for s in el:
                if self.SelectiveDynamics:
                    string += s.toStringSelectiveDynamics() + '\n'
                else:
                    string += s.toString() + '\n'
        return string[:-1]

    def setHeader(self, newHeader):
        """Set string at start of POSCAR."""
        if isinstance(newHeader, str):
            self.header = newHeader
        else:
            print "New header must be a string!"

    def setA0(self, newA0):
        """Set scaling for the cell."""
        if isinstance(newA0, float) or isinstance(newA0, int):
            self.a0 = newA0
        else:
            print "a0 must be a number!"

    def setCoordinateSystem(self, newCorD):
        """Set the coordinate system to Cartesian or Direct."""
        if newCorD[0] in ['C', 'c', 'D', 'd']:
            self.CorD = newCorD
        else:
            print "Invalid coordinate system!"

    def setLatticeVectors(self, newLVs):
        """Set the lattice vectors of the cell."""
        valid = True
        if not isinstance(newLVs[0], list):
            valid = False
        if len(newLVs) != 3:
            valid = False
        for v in newLVs:
            if len(v) != 3:
                valid = False
        if valid:
            self.latticeVectors = newLVs
        else:
            print "Incorrect lattice vectors!"

    def setElements(self, elements):
        """Set the elements in the cell."""
        self.elements = elements

    def sortElements(self):
        """Sort the elements alphabetically"""
        d = []
        for i in range(len(self.elements)):
            d += [[self.elements[i], self.elementCounts[i], self.sites[i]]]
        d.sort(key=lambda x: x[0])
        self.elements = [x[0] for x in d]
        self.elementCounts = [x[1] for x in d]
        self.sites = [x[2] for x in d]

    def setSiteMobilities(self, xFree, yFree, zFree):
        """Set the mobility of each site."""
        for i in range(self.numberOfElements()):
            for j in range(self.numberOfAtomsOfElement(i)):
                self.sites[i][j].xfree = xFree
                self.sites[i][j].yfree = yFree
                self.sites[i][j].zfree = zFree
        if (xFree and yFree and zFree):
            self.SelectiveDynamics = False
        else:
            self.SelectiveDynamics = True

    def addSite(self, element, site):
        """Add a site to the element-th sublist in self.sites."""
        notNew = True
        if element >= self.numberOfElements():
            self.sites.append([])
            self.elementCounts.append(0)
            notNew = False
        repeat = False
        if notNew:
            for s in self.sites[element]:
                if s.equals(site):
                    repeat = True
        if not repeat:
            self.sites[element].append(site)
            self.elementCounts[element] += 1

    def moveSite(self, element, nthSite, position):
        """Move the nth site of the element-th element to a new position."""
        self.sites[element][nthSite].move(position)

    def removeSite(self, element, nthSite):
        """Remove the nth site of the element-th element."""
        self.sites[element] = self.sites[element][0:nthSite] + \
            self.sites[element][nthSite + 1:]
        self.elementCounts[element] += -1
        if self.elementCounts[element] == 0:
            self.elementCounts = (self.elementCounts[0:element] +
                                  self.elementCounts[element + 1:])
            self.sites = self.sites[0:element] + self.sites[element + 1:]

    def newSiteList(self, newSites):
        """Get rid of old sites, replace with a new set."""
        if self.validSiteList(newSites):
            self.sites = newSites
            self.elementCounts = [len(x) for x in self.sites]
        else:
            print "Not a valid list of sites!"

    def depth(self, L):
        """Get maximum depth of list."""
        return 1 + max(map(self.depth, L)) if L and isinstance(L, list) else 0

    def validSiteList(self, siteList):
        """Check whether the list of sites is the right format."""
        if self.depth(siteList) == 2:
            for element in siteList:
                for s in element:
                    if not isinstance(s, Site):
                        return False
            return True
        else:
            return False

    def numberOfAtoms(self):
        """Return total number of atoms in the cell."""
        return sum(self.elementCounts)

    def numberOfElements(self):
        """Return the number of unique elements in the cell."""
        return len(self.sites)

    def numberOfAtomsOfElement(self, element):
        """Return the number of atoms of one element in the cell."""
        return len(self.sites[element])

    def readFile(self, fileName='POSCAR'):
        """Read in the lines from a POSCAR file and return a list of them."""
        p = open(fileName, 'r')
        lines = []
        # read file into list of lines
        while True:
            newline = p.readline()
            if newline == '':
                p.close()
                break
            lines += [newline]
        return lines

    def loadFromSQS(self, fileName='bestsqs.out', scale=False):
        """Read in SQS, convert to direct POSCAR."""
        # Get scaling vectors
        lines = self.readFile(fileName)
        scalings = (k.split()[0:3] for k in lines[0:3])
        scalings = map(lambda x: map(float, x), scalings)
        scaleMatrix = np.matrix(scalings)
        # Set a0, divide out of scaling by norm of first lattice vector
        if scale:
            a0 = float(np.linalg.norm(scaleMatrix[0]))
            scaleMatrix = scaleMatrix / a0
            self.setA0(a0)
        else:
            self.setA0(1.0)
        # Get lattice vectors
        lats = (k.split()[0:3] for k in lines[3:6])
        lats = map(lambda x: map(float, x), lats)
        latMatrix = np.matrix(lats)
        latMatrixInv = latMatrix.transpose().getI()  # transpose, invert
        # Scale lattice vectors
        latMatrix = latMatrix * scaleMatrix
        self.latticeVectors = latMatrix.tolist()
        # Add the sites
        for line in lines[6:]:
            xc, yc, zc, e = line.split()[0:4]
            cartCoords = [xc, yc, zc]
            cartCoords = map(float, cartCoords)
            cartCoordsTrans = np.matrix(cartCoords).transpose()
            directCoords = latMatrixInv * cartCoordsTrans
            directCoords = directCoords.flatten().tolist()[0]
            if e in self.elements:
                i = self.elements.index(e)
                self.elementCounts[i] += 1
            else:
                self.elements.append(e)
                i = self.elements.index(e)
                self.elementCounts.append(1)
                self.sites.append([])
            s = Site(directCoords, i)
            self.sites[i].append(s)
        self.CorD = 'Direct'
        # Pack data, sort by elements and re-extract
        d = []
        for i in range(len(self.elements)):
            d += [[self.elements[i], self.elementCounts[i], self.sites[i]]]
        d.sort(key=lambda x: x[0])
        self.elements = [x[0] for x in d]
        self.elementCounts = [x[1] for x in d]
        self.sites = [x[2] for x in d]
        # Set header
        header = ''
        for i in range(len(self.elements)):
            header += '%s%d' % (self.elements[i], self.elementCounts[i])
        self.header = header
        return self

    def isNumber(self, s):
        """Check if a string represents a number or not."""
        try:
            float(s)
            return True
        except ValueError:
            return False

    def loadFromPOSCAR(self, fileName='POSCAR'):
        """Read in 4/5 style POSCAR."""
        lines = self.readFile(fileName)
        self.sites = []
        # lattice vectors, scaling, comment at top of file
        latVecs = [k.split()[0:3] for k in lines[2:5]]
        self.latticeVectors = map(lambda x: map(float, x), latVecs)
        self.a0 = float(lines[1])
        self.header = lines[0].strip()
        # Check if POSCAR is VASP 4 or 5 style
        pointer = 5
        q = lines[pointer].split()[0]
        if self.isNumber(q):        # VASP 4 Style
            self.elementCounts = map(int, lines[pointer].split())
            self.elements = [chr(x + 97)
                             for x in range(len(self.elementCounts))]
        else:                   # VASP 5 Style
            self.elements = lines[pointer].split()
            pointer += 1
            self.elementCounts = map(int, lines[6].split())
        pointer += 1
        # Check for selective dynamics
        if lines[pointer][0] in ['s', 'S']:
            self.SelectiveDynamics = True
            pointer += 1
        else:
            self.SelectiveDynamics = False
        # Check for cartesian or direct
        self.CorD = lines[pointer].strip()
        pointer += 1
        # Add all the atom sites
        newSites = []
        index = 0
        for i in range(len(self.elementCounts)):
            for j in range(self.elementCounts[i]):
                L = lines[pointer].split()
                position = map(float, L[0:3])
                if len(L) == 6:  # selective dynamics tags present
                    s = Site(position, index, L[3], L[4], L[5])
                else:
                    s = Site(position, index)
                newSites.append(s)
                index += 1
                pointer += 1
            self.sites.append(newSites)
            newSites = []
        return self

    def sendToPOSCAR(self, fileName='POSCAR', style=5):
        """Make a POSCAR from the current Cell data."""
        # preamble
        string = (self.header + '\n' +
                  '%f\n' % self.a0 +
                  '%f %f %f\n' % (self.latticeVectors[0][0],
                                  self.latticeVectors[0][1],
                                  self.latticeVectors[0][2]) +
                  '%f %f %f\n' % (self.latticeVectors[1][0],
                                  self.latticeVectors[1][1],
                                  self.latticeVectors[1][2]) +
                  '%f %f %f\n' % (self.latticeVectors[2][0],
                                  self.latticeVectors[2][1],
                                  self.latticeVectors[2][2]))
        if style == 5:
            temp = (' %s' * len(self.elements)) % tuple(self.elements)
            temp = temp.lstrip()
            string += temp + '\n'
        string += ((' %d' * len(self.elementCounts)) % tuple(self.elementCounts)).lstrip()
        string += '\n'
        if self.SelectiveDynamics:
            string += 'Selective Dynamics\n'
        string += self.CorD + '\n'
        # add atom sites
        for element in self.sites:
            for k in element:
                if self.SelectiveDynamics:
                    pos = ' '.join(k.toStringSelectiveDynamics().split()[1:])
                    string += pos + '\n'
                else:
                    pos = ' '.join(k.toString().split()[1:])
                    string += pos + '\n'
        f = open(fileName, 'w+')
        f.write(string)
        f.close()

    def copyCell(self, otherCell):
        """Create a deep copy of another Cell object.
        Must be called from a new Cell object which is updated to be equal to
        the old Cell.
        """
        self.CorD = otherCell.CorD
        self.latticeVectors = otherCell.latticeVectors
        self.sites = []
        for i in range(len(otherCell.sites)):
            self.sites.append([])
            for j in range(len(otherCell.sites[i])):
                s = Site()
                s.copySite(otherCell.sites[i][j])
                self.sites[i].append(s)
        self.elementCounts = otherCell.elementCounts
        self.a0 = otherCell.a0
        self.header = otherCell.header
        self.SelectiveDynamics = otherCell.SelectiveDynamics

    def returnCopyOfCell(self):
        """Return a deep copy of self, create its own new Cell."""
        new = Cell()
        new.copyCell(self)
        return new

    def applyStrain(self, strainTensor):
        """Apply strain tensor in Voigt notation to cell."""
        e1, e2, e3, e4, e5, e6 = strainTensor
        newVectors = []
        for i in range(len(self.latticeVectors)):
            x, y, z = self.latticeVectors[i]
            a = (1 + e1) * x + (e6 / 2) * y + (e5 / 2) * z
            b = (e6 / 2) * x + (1 + e2) * y + (e4 / 2) * z
            c = (e5 / 2) * x + (e4 / 2) * y + (1 + e3) * z
            newVectors += [[a, b, c]]
        self.setLatticeVectors(newVectors)
