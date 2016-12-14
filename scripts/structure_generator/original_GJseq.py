#!/usr/bin/env python
# coding: utf-8

"""
Original script to build DGL structures.
Adapted from Gerald's script GJseq.py, for Python 3 compatibility.
This script will reproduce the experimental synthesis of the DGLs and will
output random structures of DGLs (but respecting their experimental
parameters).
"""

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


import random


class Lysine(object):

    """define a lysine residue"""

    def __init__(self, state="protected", father=None):
        self.father = None
        self.N = None   # lysine linked by peptidic N
        self.NZ = None  # lysine linked by sidechain N
        self.state = state  # default is "protected" (see chemistry)


    def isProtected(self):

        return self.state == "protected"


    def unProtect(self):

        self.state = "unprotected"


    def isNlinked(self):

        return self.N is not None


    def isNZlinked(self):

        return self.NZ is not None


    def linkN(self, lysine=None):

        if lysine is None:
            lysine = Lysine()
        self.N = lysine
        lysine.father = self

        return lysine


    def linkNZ(self, lysine=None):

        if self.isProtected is True:
            raise RuntimeError("Cannot link protected residue by NZ \
                                side chain")
        if lysine is None:
            lysine = Lysine()

        self.NZ = lysine
        lysine.father = self

        return lysine


    def __repr__(self):

        # nomenclature:
        # "Z": "LYS" # CO2-NH3-NH3 init w/o  branch (first residue)
        # "z": "LYS" # CO2-NH3-NH3 init w/o  branch (first residue) protected
        # "A": "LPC" # CO2-NH3-NH  init w/o  branch
        # "a": "LPC" # CO2-NH3-NH  init w/o  branch protected
        # "B": "LNC" # CO2-NH -NH  init with branch
        # "C": "LNP" # CO -NH3-NH  add  w/o  branch
        # "c": "LNP" # CO -NH3-NH  add  w/o  branch protected
        # "D": "LNN" # CO -NH -NH  add  with branch
        # "E": "LPN" # CO -NH -NH3 end  with branch
        # "F": "LPP" # CO -NH3-NH3 end  w/o branch
        # "f": "LPP" # CO -NH3-NH3 end  w/o branch protected
        F = False   # no father
        N = False   # no peptide bond
        NZ = False  # no branching
        P = True    # protected
        if self.father is not None: F = True
        if self.N is not None: N = True
        if self.NZ is not None: NZ = True
        if self.state != "protected": P = False
        #     (  F  ,   N  ,   NZ ,   P  )
        L = {(False, False, False, False): "Z",
             (False, False, False, True ): "z",
             (False, True , False, False): "A",
             (False, True , False, True ): "a",
             (False, True , True , False): "B",
             (True , True , False, False): "C",
             (True , True , False, True ): "c",
             (True , True , True , False): "D",
             (True , False, True , False): "E",
             (True , False, False, False): "F",
             (True , False, False, True ): "f",
             }
        return L[(F, N, NZ, P)]


class GJ():

    def __init__(self):
        self.root = Lysine()  # only one lysine as a seed

    def __repr__(self):
        return GJ.strwalk(self.root, 0, "")[0]

    def __len__(self):
        return GJ.nwalk(self.root, 0)


    @staticmethod
    def nwalk(lysine, resid):

        resid += 1
        if lysine.N is not None:
            resid = GJ.nwalk(lysine.N, resid)
        if lysine.NZ is not None:
            resid = GJ.nwalk(lysine.NZ, resid)

        return resid


    @staticmethod
    def strwalk(lysine, resid, strRep):
        resid += 1
        lysine.number = resid
        strRep += repr(lysine)
        if lysine.N is not None:
            strRep, resid = GJ.strwalk(lysine.N, resid, strRep)
        if lysine.NZ is not None:
            strRep += "%d" % (lysine.number)
            strRep, resid = GJ.strwalk(lysine.NZ, resid, strRep)
        return strRep, resid


    def Nlinkable(self):

        """return all lysine residues that can be linked through their peptidic
        N atom"""

        return GJ.linkNwalk(self.root, [])


    @staticmethod
    def linkNwalk(lysine, linkable):

        if lysine.N is None:
            linkable.append(lysine)  # nothing on N -> it is linkable
        else:
            linkable = GJ.linkNwalk(lysine.N, linkable)

        if lysine.NZ is not None:
            linkable = GJ.linkNwalk(lysine.NZ, linkable)

        return linkable


    def NZlinkable(self):

        """return all lysine residues that can be linked through their
        sidechain NZ atom"""

        return GJ.linkNZwalk(self.root, [])


    @staticmethod
    def linkNZwalk(lysine, linkable):

        if not lysine.isProtected() and lysine.NZ is None:
            linkable.append(lysine)
        if lysine.N is not None:
            linkable = GJ.linkNZwalk(lysine.N, linkable)
        if lysine.NZ is not None:
            linkable = GJ.linkNZwalk(lysine.NZ, linkable)

        return linkable


    def unProtect(self):

        """unprotect the lysine sidechain at the end of the synthesis"""

        GJ.unProtectWalk(self.root)


    @staticmethod
    def unProtectWalk(lysine):
        lysine.unProtect()
        if lysine.N is not None:
            GJ.unProtectWalk(lysine.N)
        if lysine.NZ is not None:
            GJ.unProtectWalk(lysine.NZ)


DGL = GJ()
# 1) we build G1:
L = DGL.root

# Coupling of 8 residues
while (len(DGL) < 8):
    L = L.linkN()

# Deprotection
DGL.unProtect()

# ---------------------- G2 -----------------------

# List branchable residues
branchable = DGL.NZlinkable()

# experimentally: 6 branches
# pick 6 residues to branch
branches = random.sample(branchable, 6)


# Branching !
for b in branches:
    b.linkNZ()

# Add new residues to complete the new generation
while len(DGL) < 48:

    # Random elongation
    expandables = DGL.Nlinkable()
    l = random.choice(expandables)
    l.linkN()

# Deprotection
DGL.unProtect()

# # ---------------------- G3 -----------------------

branchable = DGL.NZlinkable()
branches = random.sample(branchable, 24)

for b in branches:
    b.linkNZ()

while len(DGL) < 123:
    expandables = DGL.Nlinkable()
    l = random.choice(expandables)
    l.linkN()

DGL.unProtect()

# Replace some residues that don't physically exist
str_DGL = str(DGL)
str_DGL = str_DGL.replace('F', 'h')
str_DGL = str_DGL.replace('E', 'e')

print(str_DGL)

# ---------------------- G4 -----------------------

branchable = DGL.NZlinkable()
branches = random.sample(branchable, 54)

for b in branches:
    b.linkNZ()

while len(DGL) < 365:
    expandables = DGL.Nlinkable()
    l = random.choice(expandables)
    l.linkN()

DGL.unProtect()

# Replace some residues that don't physically exist
str_DGL = str(DGL)
str_DGL = str_DGL.replace('F', 'h')
str_DGL = str_DGL.replace('E', 'e')

print(str_DGL)
