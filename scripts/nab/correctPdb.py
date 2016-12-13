#!/usr/bin/python

import sys

"""
correct a PDB file from buildGJ.nab by adding the correct TER lines
"""

f = open(sys.argv[1])
contents = f.readlines()
f.close()

# search for CONECT
connect = []
for line in contents:
  if line.startswith("CONECT"):
    connect.append( line[6:12] )
#   print "|%s|" % (line[6:12])
#print connect

#curRes = False
#ter = False
for line in contents:
# print "|%s|" % (line[6:12])
  if line.startswith("ATOM "):
#   if line[20:26] != curRes:
#     curRes = line[20:26]
#     if ter:
#       ter = False
    if line[6:12] in connect:
        print "TER"
#     ter = True
  print line,

