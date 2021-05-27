#!/usr/bin/python
# -*- coding: utf-8 -*-

# Import system module
import sys
import math

from numpy import zeros
import numpy as NP


nr_input_arg = len(sys.argv)

try:
  infilename = sys.argv[1]
  outfilename = sys.argv[2]
except:
  print "\033[1;38mWrong number of input arguments (%d)\033[1;m" % (nr_input_arg - 1)
  print "\033[1;38mUse as %s inputfile outputfile\033[1;m" % sys.argv[0]

# Array which remembers the number of nodes for each element type defined in gmsh:
gmsh_elem_nodes = NP.array([2,3,4,4,8,6,5,3,6,9,10,27,18,14,1,8,20,15,13,9,10,12,15,15,21,4,5,6,20,35,56])

print "Starting to read the input file ..."

infile = open(infilename, 'r')

buffer = "";

## START READING THE FILE

while buffer.find('Nodes') != 1 :
  buffer = infile.readline();

NNodes = int(infile.readline())
print "\tReading \033[1;31m%d\033[1;m nodes" % NNodes


x = NP.zeros(NNodes,float)
y = NP.zeros(NNodes,float)
z = NP.zeros(NNodes,float)

for inode in range(NNodes):
  nodeinfo = (infile.readline()).split() # Line containing all coordinates
  x[inode] = float(nodeinfo[1])
  y[inode] = float(nodeinfo[2])
  z[inode] = float(nodeinfo[3])

print "\t...done"

#for inode in range(NNodes):
#  print "x[%d] = %f" % (inode,x[inode])

buffer = "";

medit_idx = NP.zeros(NNodes,int)

while buffer.find('Elements') != 1 :
  buffer = infile.readline()

NElements = int(infile.readline())

print "\tReading \033[1;31m%d\033[1;m elements - part 1" % NElements


############## FIRST READING - WE USE IT TO COUNT THE ELEMENT TYPES ############## 
# (HOW MANY LINES, TRIANGLES ETC.)

NrP1Line = 0
NrP1Tri = 0
NrP1Tet = 0

for ielem in range(NElements):
  eleminfo = (infile.readline()).split()
  elem_nr = int(eleminfo[0])
  elem_type = int(eleminfo[1])

  if elem_type == 1: # Simple line segment (P1 line)
    NrP1Line = NrP1Line + 1
 
  elif elem_type == 2: # 3-node triangle
    NrP1Tri = NrP1Tri + 1

  elif elem_type == 4:
    NrP1Tet = NrP1Tet+1

infile.close() # The first pass finished, start the second:

print "\tReading \033[1;31m%d\033[1;m elements - part 2" % NElements
infile = open(infilename, 'r')

buffer = "";

## START READING THE FILE

while buffer.find('Elements') != 1 :
  buffer = infile.readline();

buffer = infile.readline(); #This reads the number of elements again
                            #we have this information from the first reading

# Create matrices which hold connectivity:
LineConnectivity = zeros( (NrP1Line,3), int )
TriagConnectivity = zeros( (NrP1Tri,4), int )
TetraConnectivity = zeros( (NrP1Tet,5), int )

iline = 0
itriag = 0
itetra = 0

for ielem in range(NElements):
  eleminfo = (infile.readline()).split()
  elem_nr = int(eleminfo[0])
  elem_type = int(eleminfo[1])
  nr_tags = int(eleminfo[2])
  phys_tag = int(eleminfo[3])


  if elem_type == 1: # 2-node line
    LineConnectivity[iline,0] = int(eleminfo[2+nr_tags+1]) # 2 nodes of the line
    LineConnectivity[iline,1] = int(eleminfo[2+nr_tags+2])
    LineConnectivity[iline,2] = phys_tag # physical tag of triag needed for medit file format
    iline = iline + 1

  elif elem_type == 2: # 3-node triangle
    # We save the information about the nodes of the triangle - they will have
    # the same physical entity index as the triangle
    TriagConnectivity[itriag,0] = int(eleminfo[2+nr_tags+1]) # 3 nodes of the triag
    TriagConnectivity[itriag,1] = int(eleminfo[2+nr_tags+2])
    TriagConnectivity[itriag,2] = int(eleminfo[2+nr_tags+3])
    TriagConnectivity[itriag,3] = phys_tag # physical tag of triag needed for medit file format
    itriag = itriag + 1

  elif elem_type == 4: # 4-node tetrahedron
    TetraConnectivity[itetra,0] = int(eleminfo[2+nr_tags+1]) # 4 nodes of the tetra
    TetraConnectivity[itetra,1] = int(eleminfo[2+nr_tags+2])
    TetraConnectivity[itetra,2] = int(eleminfo[2+nr_tags+3])
    TetraConnectivity[itetra,3] = int(eleminfo[2+nr_tags+4])
    TetraConnectivity[itetra,4] = phys_tag # physical tag of tetra needed for medit file format
    itetra = itetra + 1

infile.close() # Second pass finished


########################## SET THE PHYSICAL INDEX OF NODES NEEDED FOR MEDIT ##########################
# FIRST, SET THE PHYSICAL INDEX OF ALL VERTICES OF ALL 3D ELEMENTS (TETRAHEDRA, HEXAHEDRA,...)
# THEN SET THE PHYSICAL INDEX OF ALL 2D ELEMENTS (TRIANGLES, QUADS) 
# FINALLY SET THE PHYSICAL INDEX OF ALL 1D ELEMENTS (LINES)

for itet in range(NrP1Tet):
  phys_tag = TetraConnectivity[itet,4]
  for node in range(4):
    node_nr = TetraConnectivity[itet,node]
    # Because python indexes an array of length N from 0 to N-1, we store the physical index for node i (in gmsh: i goes from 1 to N)
    # in location medit_idx[i-1]
    medit_idx[node_nr-1] = phys_tag 


for itri in range(NrP1Tri):
  phys_tag = TriagConnectivity[itri,3]
  for node in range(3):
    node_nr = TriagConnectivity[itri,node]
    medit_idx[node_nr-1] = phys_tag 


for iline in range(NrP1Line):
  phys_tag = LineConnectivity[iline,2]
  for node in range(2):
    node_nr = LineConnectivity[iline,node]
    medit_idx[node_nr-1] = phys_tag 



################################### ALL READING FINISHED. START WRITING. ###################################

#Open file for writing:
outfile = open(outfilename, 'w')

outstring = "MeshVersionFormatted 1\n\nDimension 3\n\nVertices\n%d\n" % (NNodes)
outfile.write(outstring)

for inode in range(0,NNodes):
  outputstr = "%0.16f %0.16f %0.16f %d\n" % (x[inode],y[inode],z[inode],medit_idx[inode])
  outfile.write(outputstr)

outputstr = "\nTriangles\n%d\n" % (NrP1Tri)
outfile.write(outputstr)

for itri in range(0,NrP1Tri):
  outputstr = "%d %d %d %d\n" % (TriagConnectivity[itri,0],TriagConnectivity[itri,1],TriagConnectivity[itri,2],
                                 TriagConnectivity[itri,3])
  outfile.write(outputstr)


outputstr = "\nTetrahedra\n%d\n" % (NrP1Tet)
outfile.write(outputstr)

for itet in range(0,NrP1Tet):
  outputstr = "%d %d %d %d %d\n" % (TetraConnectivity[itet,0],TetraConnectivity[itet,1],
                                 TetraConnectivity[itet,2],TetraConnectivity[itet,3],TetraConnectivity[itet,4])
  outfile.write(outputstr)


outfile.close()
