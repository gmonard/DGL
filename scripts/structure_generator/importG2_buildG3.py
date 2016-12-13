#!/usr/bin/env python
# coding: utf-8

"""
Imports G2 structures, and builds a G3 structure for each
"""


import random
import re



class Lysine():

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

        """return all lysine residues (as objects) that can be linked
        through their peptidic N atom"""

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

        """return all lysine residues (as objects) that can be linked through their
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
