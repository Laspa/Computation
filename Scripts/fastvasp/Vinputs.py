# =============================================================================
# Jonas Kaufman jlkaufman@hmc.edu
# June 19, 2016
# Classes and methods for making and manipulating VASP input files.
# =============================================================================
import collections


class Kpoints:
    """Hold information for VASP KPOINTS file."""

    def __init__(self):
        self.header = 'Automatic mesh'
        self.scheme = 'Gamma'
        self.auto = False
        self.length = 0
        self.subdivisions = [1, 1, 1]
        self.offsets = [0, 0, 0]

    def __repr__(self):
        return self.toString()

    def toString(self):
        """Return string representation of KPOINTS."""
        s = self.header
        s += '\n0'
        s += '\n%s' % self.scheme
        if self.auto:
            s += '\n%d' % self.length
        else:
            s += '\n%d %d %d' % tuple(self.subdivisions)
            s += '\n%d %d %d' % tuple(self.offsets)
        return s

    def setHeader(self, newHeader):
        """Set string at start of INCAR."""
        if isinstance(newHeader, str):
            self.header = newHeader
        else:
            print 'New header must be a string!'

    def setScheme(self, newScheme):
        """Set center/scheme of k-mesh."""
        if isinstance(newScheme, str):
            if newScheme:
                newScheme = newScheme[0].upper()
                if newScheme == 'A':
                    self.scheme = 'Auto'
                    self.auto = True
                elif newScheme == 'G':
                    self.scheme = 'Gamma'
                    self.auto = False
                elif newScheme == 'M':
                    self.scheme = 'Monkhorst'
                    self.auto = False
                else:
                    print 'Invalid scheme!'
        else:
            print 'New scheme must be a string!'

    def setSubdivisions(self, newSubs):
        """Set subdivisions of k-mesh."""
        # list of three elements
        if isinstance(newSubs, list):
            if len(newSubs) == 3:
                self.subdivisions = newSubs
            else:
                print 'Subdivisons must be three elements!'
        else:
            print 'New subdivisions must be list!'

    def setLength(self, newLength):
        """Set length of fully automatic k-mesh."""
        if isinstance(newLength, int):
            self.length = newLength
        else:
            print 'New length must be an integer!'

    def sendToKPOINTS(self, outFile='KPOINTS'):
        """Make a KPOINTS file from the current data."""
        s = self.toString()
        f = open(fileName, 'w+')
        f.write(s)
        f.close()

    def loadFromKPOINTS(self, inFile='KPOINTS'):
        """Read from KPOINTS file."""
        f = open(inFile, 'r')
        
        
        return self


def getSubdivisions(location='.', file='KPOINTS'):
    """ parses a KPOINTS file and pulls out the subdivisions """
    f = open('%s/%s'%(location,file),'r')
    for i in range(4):
        nextLine = f.readline()
    # Return list of subdivisions or length
    subdivisions = [int(n) for n in nextLine.split()[0:3]]
    f.close()
    return subdivisions

class Incar:
    """Hold information for VASP INCAR file."""

    def __init__(self, static=True, continuation=False):
        self.tags = collections.OrderedDict()
        self.tags['SYSTEM'] = 'INCAR'
        self.tags['PREC'] = 'Accurate'
        self.tags['ENCUT'] = 450
        self.tags['LREAL'] = '.FALSE.'
        self.tags['ISMEAR'] = -5
        self.static = static
        if static:

        self.continuation = continuation

        if continuation:

    def __repr__(self):

    def toString(self):
        """Return string representation of INCAR."""


    def loadFromINCAR(inFile='INCAR'):
    """Read INCAR file and return an ordered dictionary of tags."""

    return self

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
     

    def sendToINCAR(tags, dest='.', filename='INCAR'):
        """ Takes in ordered dictionary of tags and makes INCAR file in dest """
        s = ''
        for k, v in tags.items():
            s += '%s = %s\n'%(k,str(v))
        f = open('%s/%s'%(dest,filename),'w+')
        f.write(s)
        f.close()


    def changeTags(self, newTags):
        """"""
        return

    def inputTags(self):
        """Add/change tags from user input."""
        return

