#!/usr/bin/python

# Copyright (C) 2016  Jean-Patrick Francoia, Jean-Christophe Rossi,
# Gerald Monard, and Laurent Vial

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""
correct a PDB file from buildGJ.nab by adding the correct TER lines
"""


import sys

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

