#!/usr/bin/env python
f = open("OUTCAR")
data = f.read()
data = data.split("ELASTIC MODULI CONTR FROM IONIC RELAXATION (kBar)")[1]
#print data
data = data.split()
#print data
C11 = float(data[10])
C12 = float(data[11])
C44 = float(data[34])
print "C11 is "+data[10]
print "C12 is "+data[11]
print "C44 is "+data[34]
B = C11/3 + 2*C12/3
G = (3*C44+C11-C12)/5
print "The bulk modulus (B) is: "+ str(B)
print "The shear modulus (G) is: "+ str(G)
print "The ratio of bulk to shear modulus (B/G) is: "+ str(B/G) 
