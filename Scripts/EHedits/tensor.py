#!/usr/bin/env python
f = open("OUTCAR")
data = f.read()
data = data.split("TOTAL ELASTIC MODULI")[1]
#print data
data = data.split()
#print data
C11 = float(data[10])
C12 = float(data[11])
C22 = float(data[18])
C33 = float(data[26])
C44 = float(data[34])
C55 = float(data[42])
C66 = float(data[50])
trace = C11 + C22 + C33 + C44 + C55 + C66
print "The trace of the elastic tensor in GPa is: " + str(trace/10)
B = C11/3 + 2*C12/3
G = (3*C44+C11-C12)/5
A = 2*C44 / (C11-C12)
nu = C12 / (C11 + C12)
CP = C12-C44
print "The bulk modulus (B) in GPa is: "+ str(B/10)
print "The shear modulus (G) in GPa is: "+ str(G/10)
print "The ratio of bulk to shear modulus (B/G) is: "+ str(B/G)
print "The Poisson ratio (nu) is: " + str(nu)
print "The shear anisotropy factor (A) is: " + str(A)
print "The Cauchy pressure in GPa is: " + str(CP/10) 
