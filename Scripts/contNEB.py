#!/usr/bin/env python
#===============================================================================
#  Simone Griffith      sgriffith@hmc.edu
#  November 10, 2018
#  Script to quickly continue NEB runs on Stampede
#===============================================================================
""" Add POSCAR, POTCAR, KPOINTS and INCAR files to
the working directory. """
import subprocess as sp
import collections

#===============================================================================
# INCAR Manipulation 
#===============================================================================
def readINCAR(filename='INCAR'):
        """ Reads INCAR file and returns dictionary
 of tags """
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

def makeINCAR(tags, dest='.', filename='INCAR'):
        """ Takes in ordered dictionary of tags and
 makes INCAR file in dest """
        s = ''
        for k, v in tags.items():
                s += '%s = %s\n'%(k,str(v))
        f = open('%s/%s'%(dest,filename),'w+')
        f.write(s)
        f.close()

#===============================================================================
# File Preparation
#===============================================================================

def prepareINCARcont():
        """ Prepares the INCAR for a continued run"""
""
        tags = readINCAR() #Read INCAR file

        # change INCAR values for continuation
        tags['ISTART'] = 1
        tags['ICHARG'] = 0
        makeINCAR(tags,'continued')

def genContScript(name,copy,dest='.'):
        """ 
        Creates script name_batch in dest to move files over 
        to the continued directory, remove the POSCAR and fe.dat,
        and rename CONTCAR as POSCAR
        """
        s = '#!/bin/bash'

        s += '\ncp KPOINTS continued'
        s += '\ncp submit continued'

        for i in range(len(copy)):
                s += '\ncp -r %s' %copy[i]
                s_len = len(s)+1
                s = s.ljust(s_len)
                s += 'continued'

        copy = copy[1:-1]
        s += '\ncd continued'

        for i in range(len(copy)):
                s += '\ncd %s' %copy[i]  #enter working directory
                s += '\nrm POSCAR fe.dat' #remove old files
                s += '\nmv CONTCAR POSCAR' #rename CONTCAR as POSCAR
                s += '\ncd ..' #exit working directory

        filename = '%s/%s_batch'%(dest,name)
        f = open(filename,'w+')
        f.write(s)                                                              # write file
        f.close()
        sp.call(['chmod','+x',filename])                # make executable


#===============================================================================
#  Main Program
#===============================================================================

images=raw_input('How many images is this run?')
print images, '\n'

sp.call(['mkdir','continued'])
# print 'mkdir'

copy = []
dirNum = 0
images = int(images)

while dirNum <= images+1:
        dirN = str(dirNum)
        dirName = '0'
        dirName += dirN
        copy += [dirName]
        dirNum = dirNum + 1

#create continued script
genContScript('cont',copy)

#edit INCAR
prepareINCARcont()
