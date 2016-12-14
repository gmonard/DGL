#!/usr/bin/env python
# coding: utf-8

"""
Imports G2 structures, and builds a G3 structure for each
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
import re
from original_GJseq import Lysine, GJ


# ---------------------- Import G2 -----------------------

# List of the G2 structures we used
structures = [
              "BCDCDDDDCCh8CCCCCh7CCCCCCCCCh6CCCCCh5CCCCh3CCCCh1CCCCh",
              "BDDDCDDCCCCCCCCh7CCCh6CCh4CCh3CCCCCCCCh2CCCCCCh1CCCCCh",
              "ACDDDDDDCCCCh8CCCCCCCCh7CCCCCCCCh6CCCh5CCh4CCCCCCh3CCh",
              "BDCDDDDCCCCh7CCCCCCCCh6CCCCh5CCCCCh4CCCCh2CCCCCCh1CCCh",
              "BDDDCCDDCCCCCCh8CCCCh7CCCCCCCh4CCCh3CCCCCCCh2CCCh1CCCh",
              "ADDDDCDDCCCCCCh8CCCCh7CCCCCCh5CCCCCCCh4CCCCh3CCCCh2CCh",
              "BDCDCDDDCCh8CCCCCCh7CCCCCCCh6CCCCCCh4CCCCh2CCCh1CCCCCh",
              "BDDCDDDCCCCh7CCCCCh6CCCCCCCh5CCh3CCCCCCCCCh2CCCCCh1CCh",
              ]

for inp in structures:

    nbr_res = 0
    dico_res = {}

    chains = re.findall("([a-zA-Z]+)", inp)
    bp = re.findall("([0-9]+)", inp)

# chains: ['ACDDDDDDCCCCh', 'CCCCCCCCh', 'CCCCCCCCh', 'CCCh', 'CCh', 'CCCCCCh', 'CCh']
# bp: ['8', '7', '6', '5', '4', '3']

    dgl = GJ()

    # Create first residue
    L = dgl.root
    nbr_res += 1
    dico_res[nbr_res] = L

    # Build the core
    while len(dgl) < len(chains[0]):
        L = L.linkN()
        nbr_res += 1
        dico_res[nbr_res] = L
    dgl.unProtect()

    for number, sequence in zip(bp, chains[1:]):
        # Branch the residue at position 'number'
        L = dico_res[int(number)].linkNZ()
        nbr_res += 1
        dico_res[nbr_res] = L

        # Elongate after branching. -1 residue bc already branched
        for res in sequence[:-1]:
            L = L.linkN()
            nbr_res += 1
            dico_res[nbr_res] = L

    # Deprotects the residues
    dgl.unProtect()

    str_DGL = str(dgl)
    str_DGL = str_DGL.replace('F', 'h')
    str_DGL = str_DGL.replace('E', 'e')

    # Check if the importation succeeded
    if str_DGL != inp:
        print("ERROR import structure")

    # Check that the input string and the new tree built match
    print(inp)
    print(str_DGL)


# ---------------------- G3 -----------------------

    # Pick random residues to branch
    branchable = dgl.NZlinkable()
    branches = random.sample(branchable, 24)

    for b in branches:
        b.linkNZ()

    while len(dgl) < 123:
        # Random elongation
        expandables = dgl.Nlinkable()
        l = random.choice(expandables)
        l.linkN()

    dgl.unProtect()

    str_DGL = str(dgl)
    str_DGL = str_DGL.replace('F', 'h')
    str_DGL = str_DGL.replace('E', 'e')

    print(str_DGL)
