#!/usr/bin/python
# coding: utf-8

"""This script will generate the 2D structure a DGL from its 1D sequence.
See fig 2 of the paper."""

# https://networkx.readthedocs.io/en/stable/examples/drawing/circular_tree.html


import networkx as nx
import matplotlib.pyplot as plt
import re
import pygraphviz
from networkx.drawing.nx_agraph import graphviz_layout



# List of 1D structures
inp = [
       "BDDDDDDDDDCCDDCh14CCCh13Ch10Ce24CCh9DCCCCh28CCh8DDDDDDDCCCCCh43Ch42CCCCh41DDDh59Ch58CCCh57CCh40CCDCCh72CCCCCCh39CDDCCh85Ch84h38Ch37DDDCh96Ch95h94CCCh7CDDDDDCDCDh115CDCCh118CCh113Dh125Ch111CDe131CCCh130h110Ch109DDCCCCh140CCh139CCCCh108CCCh107DDh159CCCh158CCCh6DDDDDDCh174CCCh173DCDCDDCCCh186CCh185CCh183CCh181Ch172Ch171DDDCCCh206CCh205CCh204Ch170CCCCh169CCCh5DDDCCCh230DDe237CCCh236CCCh235Ch229CCCCh228CCCCCCCh4CCCh3DDDCDDCCh270Ch269DDh277h276CCh267CCCh266DDCCh288CCCh287CCCCCh265DCCCCh302CCh2CCCCCh1DCCDDDh322Ch321DDCe329CCCCCCCCh327CCh326CCCCCh320CCDCCCh350CCCh317DCe361Ch359Ch",
       "BDDDDDDDDCDDCDCDh16DCh18CCCCh14DDCDDDDh32CCh31h30CCh29CCh27Ch26CCh12Ch11Dh51CCCCh9CCCCh8CCCCh7DDDCDCh72Ch70DDDCh79CCCh78CCh77CCh69CCh68CCCCh6CDDCCh102CDh107CCh101DCCh112Ch5CCDCCCh120CCCh4DDDDCCh132CCh131CCCCCh130Ch129DDCCh148h147Ch3CDCCDDDDDDDCCCh165CCCCCCCh164CCCh163CDCh183Ch162CCe190CCCh161DCh195Ch160CCCh159Ch156CCCCh2DDDDDDDDDh219CCCh218Ch217DDCCh228CCh227CCCCCh216CCDDCCCCh244CCCh243Ch215CCCh214e260CCh213CDDCCCCh266CCCh265CCCCh212CCDDCh284CCCh283CCh211DCCCCh294CCh1DCDCDDDDCCh310Ch309Ch308CCh307DDh322h321CCh305DDCe331Ch329CCCCh328CCCh303DDCDDCh347CCh346CCCh344CCCh343CCCCh",
       "ADDDDDDDDDDDDDCDDh17h16CCCCh14CCh13DDCCCCh29CCCh28Ch12CDe43CCh42Ch11CCh10Ch9DCDCh56CCh54Ch8CCCDDDDCDCh72DCCCCh75CCh70CCh69DCDCCCCh89Ch87CCCh68CCDDCCCh104CCh103Ch67CDDDCh117Ch116CCh115Ch7DDDDDDDDDCh135CCh134CCCh133CCCh132CCh131CDDCCCCh154CCh153CCh130Dh166CCCh129CDDCh174CCCh173CCCh128CDCCh186CCh127CCCCCCCCh6DDDDDCh206CCCh205DDCCh214CCh213Ch204DCDDh226CCCCh225CCh223CCCh203DCCh240CCCCh202DCh249CCh5DDCDDCh259CCh258CCCCh256DDDCCCh272CCh271CCCCh270Ch255CCh4CCDDCDDCCCCCCh296CCh295CCh293CCCh292DDe316h315CCh314CCCCCCh3DDDDh331CCCh330DCCCh337Ch329DCh344h328CCCCh2DDDh355CCh354CCh353CCh",
       "BDDDDDDDDDDCDCh13CCCCh11CDDh23CCCCCh22CCCh10CCDCh37Ch9Ch8CDCDCCCh47h45Ch7DDDDDDDDDDh64Ch63Ch62CCh61CDCh74CCh60Ch59CCCh58CCh57CCCCh56CDDDh97Ch96CCh95CCCCh55CDCh110CCh6CDDDDh120DDCCh123Ch122Ch119CCCh118CDDCCh137CCCh136CCh117CCCCCCh5DDDDDDh160DDDDh165CCh164CCCh163CCCh162CCCh159CCh158CCCCCCh157DCCCh192CCCh156CCCh155CCDDCDCCCh210CCCCCCh208h207Ch4DDDDDCDCh231CCCCCCCCh229h228h227CCe247CCCh226CDDDh255CCCh254Ch253Ch225h3CCh2DDDDDDDh275CCh274DDCDh283CCh281CCCCh280CCCh273DDDCh299CCh298CCh297CCh272CCCh271CCh270CCh269DCCh321Ch1DDDCCCCCh329CCCh328CCDDDCCCCh344CCh343CCCCh342Ch327CCCCCh",
       "BDDDDCDDDCCDDCDCCh15CDCCCCCCCh20h13DDDCh32Ch31CCh30CCh12CDDCCCh45Ch44CCCh9DCCCCh56CCCCh8DCDDh70CCh69CCh67CCCCh7DDCDDDDDCh90CDCh94CCCCh89De103CCCCCh102CCh88Ch87CCDh117CCCh86DCh123Ch84CCh83DDDCh133Ch132CCCh131h5Ch4DDDDCCh148Dh152Ch147DDh157CCCCCh156CCh146De169CCh168Ch145CCDCCCCh177CCCh3DDDDDDDDDCh195Ch194Ch193CCDDDh206CCCh205CCCCCh204CCh192CCCDDCCh225CCh224CCCh191Dh236CCCCh190CCCCCh189DDh250CCh249CCCCh188CCDDDCCCh264CCCh263Ch262CCh187DDCCh279CCh278CCCCCh2DDCDCCh295CCh293CCCh292DDCDCCCCh309CCh307CCCCh306CCCCh1DDDDDCh332h331DDCCh337Ch336h330CDDDh347CCh346CCCh345CCCh329CCCh328Ch",
       "BDDDDDDDDDDDDCCCCCCh13h12CCh11DDDCh27CCh26CCh25Ch10CCCh9CCh8DDDDDCCCh49h48DDDCh57CCh56CCCCh55h47CCCh46DCDCCh75CCh73CCCh45DCCCCCh86h7DDDCDDDDDCCDh105h102CCCh101CCCCh100DCh117CCCCh99DDCh126CCh125CCCh98DDDCCCCCCh138Ch137CCCCh136CCCCCCh96CCCh95DDCDCCCh167Ch165Ch164CCCh94DCDCh182Ch180CCCh6CCDDh194Ch193CCCh5DCDDDDDDDh210Ch209DCCCh214CCCh208DCDCCh225Ch223Ch207CCCCCh206DDh240Ch239CCCh205CDh249CCh204h202DDCDCCh258CCh256CCh255CCCCCh4DDDDDCCCCh278CCCh277CCDCCh290Ch276CDDCCh298CCCh297CCh275CCCh274DCh313CCCCCCh3DDDCDh327CCDCh331CCh325CCCh324DDh342CCh341CCh323Ch2DCDh354CCh352CCCCh1Ch",
       "BDCDDDDDDCDCh11CCDCDCCCCh18CCCCh16CCCCh9CCCh8DDDDDDDCh44DDDCCh49CCh48CCh47CCh43CCh42CDCh66Ch41DCh71CCh40CDh78CCh39CCh38Ch7DDDDDDDDDCh96CCCh95Ch94Dh105CCCCCCh93DDDDCCh117CCCCCCh116CCCh115CCCh114Ch92CCDCCh140CCCh91CDCh149CCCCh90CCCh89CCCDCCCh164Ch88Ch6DDDDDDDCCCCh179CCh178Ch177Ch176Dh192Ch175DCCCh196Ch174DDCDCh206CCCh204h203Ch173CCh5CDDCCh221h220h4DDDDDh231CCCCh230CDDDCCh241Ch240Ch239CCCh229CCCh228CCCh227DDCCDCh265CCh262CCCh261CCCCh2DCDCDh284CCCh282CCh280DDh294h293CCCh1DDCDDDDDh308CCCh307CCCh306CDCDCCh321CCCCCh319h305CDCh333Ch304DDDCCCCh340CCh339CCCh338CCCCh302CCCCCh301Ch",
       "BDDCDDDDDCDDDCh13CCh12DDCDh22CCCh20CCCh19CCCCCh11CDDCh40CCh39CCCh9Ch8CCh7DDDDDDCCCh60DCCCCCCh65CCCh59CCCh58DDDCh83h82CCCCh81CCh57Ch56DDDDh100CCCh99Ch98CCh97Ch55CDCCh114CCh6DDDDDDDDh128Ce131h127CCCh126CCCCCCh125DCCCCCCh144CCh124DCCCCh155CCCCh123e166CCCCh122CCh121DDDCh177CCh176Ch175Ch5DDDh189h188CCCCh187CCCh3DDDCDDDCDDDCDCh213CCCh211CCh210CCCh209DCDCCCh229CCh227Ch207CCh206CCh205CCCCCCh203Ch202Dh254CCCh201DDDDDh264Ch263Ch262Ch261CCCh260CCh2CDDDCCDDDCCh287CCCh286CCCh285CCCh282CDCCCh304CCCCh281CCCDh317CCCCh280CCh1DDDDCh330Ch329CDh336CCCh328DDCCh343CCCCCCh342CCCh327DDCh359Ch358Ch",
      ]

for struct in inp:

    G = nx.Graph()
    G.add_node(1, G=1)

    bp = [int(value) for value in re.findall("([0-9]+)", struct)]
    chains = re.findall("([a-zA-Z]+)", struct)

    # Build the "core" first
    while G.number_of_nodes() < len(chains[0]):

        nbr_nodes = G.number_of_nodes() + 1
        G.add_node(nbr_nodes)
        G.add_edge(G.number_of_nodes() - 1, G.number_of_nodes())


    # Build the branches
    for number, sequence in zip(bp, chains[1:]):

        G.add_node(G.number_of_nodes() + 1)
        G.add_edge(number, G.number_of_nodes())
        seq = sequence[:-1]

        for i, res in enumerate(seq):
            G.add_node(G.number_of_nodes() + 1)
            G.add_edge(G.number_of_nodes() - 1, G.number_of_nodes())

    print(G.number_of_nodes())

    # Name the residues with their type
    dico_res = {'A': 'CBG',
                'B': 'CHG',
                'C': 'IBG',
                'D': 'IHG',
                'e': 'IHD',
                'h': 'IBD'
                }

    labels = {}
    for i, res in enumerate(''.join(chains)):
        labels[i + 1] = "{}:{}\n{}".format(i + 1, res, dico_res[res])

    pos = nx.nx_agraph.graphviz_layout(G)
    fig = plt.figure(1)
    nodes = nx.draw_networkx_nodes(G, pos, node_size=2000, node_color='w')

    # Set edge color to black
    nodes.set_edgecolor('k')
    # nodes.set_linewidth(20)

    nx.draw_networkx_edges(G, pos)
    nx.draw_networkx_labels(G, pos, labels, font_size=11)

    fig.set_size_inches(50, 50)
    plt.axis('off')
    plt.savefig('{}.png'.format(inp.index(struct) + 1), bbox_inches='tight')
    plt.clf()