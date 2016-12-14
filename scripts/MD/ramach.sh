#!/bin/bash

#Create a cpptraj script that will be used to generate Ramachandran maps

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

echo "parm prmtop-mass" > cpp_ramach.in
echo "trajin nvt*.crd" >> cpp_ramach.in

#phi
echo "printDihedrals @C @CA @N @C" | parmed.py -p prmtop-mass | grep CA | awk '{print "dihedral @"$1 " @" $5 " @" $9 " @" $13 " out dihedral/" $5}' >> cpp_ramach.in

#psi
echo "printDihedrals @N @C @CA @N" | parmed.py -p prmtop-mass | grep CA | awk '{print "dihedral @"$1 " @" $5 " @" $9 " @" $13 " out dihedral/" $9}' >> cpp_ramach.in

#psi branching
echo "printDihedrals @NZ @C @CA @N" | parmed.py -p prmtop-mass | grep CA | awk '{print "dihedral @"$1 " @" $5 " @" $9 " @" $13 " out dihedral/" $9}' >> cpp_ramach.in
