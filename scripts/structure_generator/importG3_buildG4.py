#!/usr/bin/env python
# coding: utf-8


"""
Import G3 structures, and build a G4 structure for each
See importG2_buildG3 for full comments, this script basically does the
same thing.
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


# ---------------------- Import G3 -----------------------

structures = [
              "BCDCDDDDDDCCCh10Ch9CCh8DCDDDCCCh24CCh23CCh22CCCCh20CCh7CDDDCDCDCe52CCh50h48CCh46CCh45h44Ch6CCDCDe71h70CCCCCh68CCCh5CDDCh85CCh84h3DDDCDh96Ch94Ch93CCh92h1DCCDDh110CCCh109CCCCh106CCh",
              "BDDDDDDCCCDDCDCDh16h14CCCCCCh12h11h7DCDCh30CCh28h6CDDh39Ch38h5CCh4DCDh49Ch47CCh3CCCCCDDDDCCh64Ch63CCh62Ch61CCh2DDDDCDDCh84Ch83CCCh81h80CCCh79CCCh78Ch1DCDCDCCh108Ch106CCCh104CCCCCh",
              "ADDDDDDDDCDDDCCCh13Ch12CCh11h9CCh8CCCDDDCCDCh35h32CCCh31CCCh30CCCh7CDDDDCDCDCh59h57h55CCh54Ch53CCh52Ch6DDDDh77Ch76CCCh75h74Ch5CDCCh89CCh4CCDCCCCCh98CCh3DDDh110CCh109CCh108CCh2CCh",
              "BDCDDDDDCDDCh11CCCh10CCh8CCCh7DDCCCCDCCh31Ch26CCCh25Ch6CDDDe47Ch46Ch45CCh44Ch5DDDCCe62CCCh59Ch58Ch57CCCCCh4CDDDCCh80h79CCh78CCCh2DCDCDDh97CCCh96CCh94Ch92CCh1DDDCh113h112CCCCh111Ch",
              "BDDDCCDDDCCDDCDh15Ch13CCh12CCh9Ch8CCCCh7DDCDDCDe39Ch38Ch36CCh35Ch33h32CCh4DDDDh56h55Ch54Ch53CCh3DDDCDDDCh72CCCCh71CCCCh70Ch68CCh67CCCCCh66Ch2DDCCCh99Ch98CCCCh1CCDDh114CCCh113CCCh",
              "BDDDDDDDCCDCCCCh11CCh8DDCDCh23CCh21CCh20Ch7DDCCDDDCCCCh40Ch39CCh38CCh35CCCh34CCCh6CCCh5DCCDDDDDh73Ch72CCh71CCh70Ch69Ch66CCCh4DCDDDh95h94CCh93CCh91CCh3CDCCDh111CCh108Ch2CCe120Ch1h",
              "BDCDDDDDDCDh11CCCCCh9h8CCDDDCDh26CCh24Ch23h22CCh7CDDDDDDDh44h43h42CCCh41CCh40Ch39Ch38CCCCh6CDDDDCCCh68h67h66h65CCCh5CCCh4DCCDh87CCCh84CCCCCh2DCCCh99Ch1CCCDDDCh111CCCh110Ch109CCCh",
              "BDDCDDDCDCDDCh12CCCh11CCh9Ch7DDCDCDCh29CCh27CCh25CCCh24Ch6DCDDDDCe51Ch49CCh48CCCh47h46h44CCCh5CDh68Ch3DDCCCDCCDDCCh81Ch80CCh77h73h72CCCCh2CCDDCCCCh100CCh99CCCh1DDDh115CCh114Ch113Ch",
              ]

for inp in structures:

    nbr_res = 0
    dico_res = {}


    chains = re.findall("([a-zA-Z]+)", inp)
    bp = re.findall("([0-9]+)", inp)
    # bp = [int(value) for value in bp]


# ['ACDDDDDDCCCCh', 'CCCCCCCCh', 'CCCCCCCCh', 'CCCh', 'CCh', 'CCCCCCh', 'CCh']
# ['8', '7', '6', '5', '4', '3']

    dgl = GJ()
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

    dgl.unProtect()

    str_DGL = str(dgl)
    str_DGL = str_DGL.replace('F', 'h')
    str_DGL = str_DGL.replace('E', 'e')

    if str_DGL != inp:
        print("ERROR import structure")

    # break
    # print(inp)
    # print(str_DGL)


# ---------------------- G4 -----------------------

    branchable = dgl.NZlinkable()
    branches = random.sample(branchable, 54)

    for b in branches:
        b.linkNZ()

    while len(dgl) < 365:
        expandables = dgl.Nlinkable()
        l = random.choice(expandables)
        l.linkN()

    dgl.unProtect()

    str_DGL = str(dgl)
    str_DGL = str_DGL.replace('F', 'h')
    str_DGL = str_DGL.replace('E', 'e')

    print(str_DGL)
