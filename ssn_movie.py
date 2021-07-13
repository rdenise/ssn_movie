#! /usr/bin/env python3
# -*- coding: utf-8 -*-

##########################################################################################
##########################################################################################
##
##                                Library
##
##########################################################################################
##########################################################################################

import argparse
from textwrap import dedent
import sys, os
import re 
import networkxgmml
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.cm as cm
import matplotlib.colors as colors
from itertools import cycle
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

##########################################################################################
##########################################################################################
##
##                                Function
##
##########################################################################################
##########################################################################################

def create_folder(mypath):

    """
    Created the folder that I need to store my result if it doesn't exist
    :param mypath: path where I want the folder (write at the end of the path)
    :type: string
    :return: Nothing
    """

    try:
        os.makedirs(mypath)
    except OSError:
        pass

    return

##########################################################################################
##########################################################################################

def get_color_cmap(name, n_colors=6):

    """
    Return discrete colors from a matplotlib palette.

    :param name: Name of the palette. This should be a named matplotlib colormap.
    :type: str
    :param n_colors: Number of discrete colors in the palette.
    :type: int
    :return: List-like object of colors as hexadecimal tuples
    :type: list
    """

    brewer_qual_pals = {"Accent": 8, "Dark2": 8, "Paired": 12,
                        "Pastel1": 9, "Pastel2": 8,
                        "Set1": 9, "Set2": 8, "Set3": 12, 'tab20':20, 'tab20b':20}


    if name == 'tab20' and n_colors > 19:
        second = 'tab20b'
        ncolor2 = n_colors - 19
        n_colors = 19
    else :
        second = False

    cmap = getattr(cm, name)
    
    if name in brewer_qual_pals:
        bins = np.linspace(0, 1, brewer_qual_pals[name])
        if 'tab20' == name :
            len_bins = len(bins)
            bins = [bins[i] for i in range(len_bins) if i != 14][:n_colors]
        else :
            bins = bins[:n_colors]
    else:
        bins = np.linspace(0, 1, n_colors + 2)[1:-1]

    palette = list(map(tuple, cmap(bins)[:, :3]))

    if second :
        cmap = getattr(cm, second)
        bins = np.linspace(0, 1, brewer_qual_pals[second])[:ncolor2]
        palette += list(map(tuple, cmap(bins)[:, :3]))

        pal_cycle = cycle(palette)
        palette = [next(pal_cycle) for _ in range(n_colors+ncolor2)]
    else :
        pal_cycle = cycle(palette)
        palette = [next(pal_cycle) for _ in range(n_colors)]

    return [colors.rgb2hex(rgb) for rgb in palette]

##########################################################################################
##########################################################################################

def visu_graph(graph, output, threshold, all_hit_id, annot_df, hit_id2node) :
    
    """
    Visualisation of the graph using the weight for the color. If weight >0.5
    put the edge red. And write the graph in graphMl format.

    :params adjacency_mtrix: adjacency matrix calculate with get_matrix_interaction_system()
    :type: pandas.DataFrame
    :params output: Name of the graphml file
    :type: str
    :return: Nothing
    
    """ 
    
    annot_df = annot_df[annot_df.Hit_Id.isin(all_hit_id)]
    all_gene = annot_df.Gene.unique()
    num_gene = all_gene.shape[0]
    all_gene = sorted(all_gene)

    palette = get_color_cmap("tab20", n_colors = num_gene)
    tmp_color = {all_gene[i]:palette[i] for i in range(num_gene)}
    
    hit_id2gene = annot_df.set_index('Hit_Id').Gene.to_dict()
    
    dict_color = {}
    for hit_id in all_hit_id :
        if hit_id in hit_id2gene :
            dict_color[hit_id2node[hit_id]] = tmp_color[hit_id2gene[hit_id]]
        else :
            dict_color[hit_id2node[hit_id]] = 'grey'

    edge2remove = []
    
    # Parsing the edges
    for n1, n2, edge_dict in graph.edges.data() :
        if edge_dict['alignment_score'] >= threshold :
            graph.edges[(n1, n2)]["color"] = "lightgrey"
        else :
            edge2remove.append((n1, n2))
                
    for e in edge2remove :
        graph.remove_edge(*e)

    graph = graph.to_undirected()
    # print("Calculating layout...")    

    # Color edges
    edges,edge_colors = zip(*nx.get_edge_attributes(graph,'color').items())   

    # Choose between : dot, neato, fdp, sfdp, twopi, circo
    pos=graphviz_layout(graph, prog="neato")

    # Put the color of the node
    nx.set_node_attributes(graph, dict_color, "color")    

    # Color nodes
    nodes,nodes_colors = zip(*nx.get_node_attributes(graph,'color').items())   
    
    # Write the graph
    # nx.write_graphml(graph, output)
 
    plt.figure(figsize=(12,10))

    # If you have too much node maybe reduce the node_size to 50
    nx.draw_networkx_nodes(graph, pos, node_color=nodes_colors, node_size=50, edgecolors="black")

    nx.draw_networkx_edges(graph,pos,edgelist=edges, edge_color=edge_colors)

    custom_lines = []
    custom_text = []
    
    for gene, color in tmp_color.items() :
        custom_text.append(gene)
        custom_lines.append(Line2D(range(1), range(1), color="white", marker='o', markerfacecolor=color, linestyle='none'))
    
    plt.legend(custom_lines, custom_text, bbox_to_anchor=(1.05, 1), loc='upper left', prop={"size":'xx-small'})
    
    #Label drawing as well
    # nx.draw_networkx_labels(graph,pos,font_size=8)

    plt.axis('off')
    plt.tight_layout()
    
    plt.title(f"SNN organisation {threshold}")
    
    if output!=None: plt.savefig(output.replace("graphml", "pdf"), dpi=300, bbox_inches='tight')

    plt.close('all')

    return

##########################################################################################
##########################################################################################
##
##                                Main
##
##########################################################################################
##########################################################################################


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
     description=dedent("""See the effect of alignement threshold based on annotation color""") )

general_option = parser.add_argument_group(title = "General input dataset options")
general_option.add_argument("-g",'--xgmml',
                            metavar="<XGMML>",
                            dest="xgmml",
                            help="XGMML file from the analysis of EFI-EST : https://efi.igb.illinois.edu/efi-est/",
                            required=True)
general_option.add_argument("-k",'--kofam',
                            metavar="<KOFAM>",
                            dest="kofam",
                            help="KOFAM file from the analysis of kofamscan")
general_option.add_argument("-e",'--eggnog',
                            metavar="<EGGNOG>",
                            dest="eggnog",
                            help="EGGNOG file from the analysis of eggnog-mapper")
general_option.add_argument("-a",'--annotation',
                            metavar="<annotation>",
                            dest="annotation",
                            help="Tabulated file with your own annotation, need to have a header with the columns 'Hit_Id' and 'Gene' where hit_id id the name of the sequence in your network and gene is the annotation on the sequence")
general_option.add_argument("-o",'--output',
                            default=None,
                            dest="output",
                            metavar='<OUTPUT>',
                            help="Name of the output file (default: [NAME_OF_XGMML])")

##########################################################################################

args = parser.parse_args()

##########################################################################################

if args.output :
    OUTPUT = args.output
else :
    OUTPUT = os.path.join(os.getcwd(), os.path.basename(args.xgmml))

##########################################################################################
print()

create_folder(OUTPUT)

XGMML = args.xgmml
KOFAM = args.kofam
EGGNOG = args.eggnog
ANNOT = args.annotation

##########################################################################################

with open(XGMML, 'rb') as r_file :
    G = networkxgmml.XGMMLReader(r_file)

all_alignment_score = set([edge_dict['alignment_score'] for n1, n2, edge_dict in G.edges.data()])
all_hit_id = [node_dict['Description'][0].split()[0] for n1, node_dict in G.nodes.data()]
hit_id2node = {node_dict['Description'][0].split()[0]:n1 for n1, node_dict in G.nodes.data()}

all_score = sorted(all_alignment_score, reverse=True)
num_score = len(all_score)
max_str_len_score = len(str(int(max(all_score))))


print(f'Network range of alignment score :: Min = {min(all_score)}, Max = {max(all_score)}')
##########################################################################################

if KOFAM :

    create_folder(os.path.join(OUTPUT, 'KOFAM'))

    kofam = pd.read_table(KOFAM)
    kofam = kofam[kofam.Hit_Id.isin(all_hit_id)].reset_index(drop=True)



    for score_index in range(num_score) :
        score = all_score[score_index]

        score_txt = str(int(score)).zfill(max_str_len_score)

        visu_graph(graph = G.copy(),
                   output = os.path.join(OUTPUT, 'KOFAM', os.path.basename(args.xgmml).replace('.xgmml', f'.{score_txt}.png')),
                   threshold = score,
                   all_hit_id = all_hit_id,
                   annot_df = kofam, 
                   hit_id2node = hit_id2node)

        print(f'Score done for kofam ::: {score_index + 1}/{num_score} : {round(((score_index + 1)/num_score)*100, 2)}%', end='\r', flush=True)

##########################################################################################

if EGGNOG :

    create_folder(os.path.join(OUTPUT, 'EGGNOG'))

    eggnog = pd.read_table(EGGNOG, comment='#')
    eggnog = eggnog.rename(columns={'query':'Hit_Id', 'Preferred_name':'Gene'})
    eggnog = eggnog[~(eggnog.Gene == '-')]

    for score_index in range(num_score) :
        score = all_score[score_index]

        score_txt = str(int(score)).zfill(max_str_len_score)

        visu_graph(graph = G.copy(),
                   output = os.path.join(OUTPUT, 'EGGNOG', os.path.basename(args.xgmml).replace('.xgmml', f'.{score_txt}.png')),
                   threshold = score,
                   all_hit_id = all_hit_id,
                   annot_df = eggnog, 
                   hit_id2node = hit_id2node)    

        print(f'Score done for eggnog ::: {score_index + 1}/{num_score} : {round(((score_index + 1)/num_score)*100, 2)}%', end='\r', flush=True)

##########################################################################################

if ANNOT :

    create_folder(os.path.join(OUTPUT, 'ANNOTATION'))

    annot = pd.read_table(ANNOT, comment='#')

    for score_index in range(num_score) :
        score = all_score[score_index]

        score_txt = str(int(score)).zfill(max_str_len_score)

        visu_graph(graph = G.copy(),
                   output = os.path.join(OUTPUT, 'ANNOTATION', os.path.basename(args.xgmml).replace('.xgmml', f'.{score_txt}.png')),
                   threshold = score,
                   all_hit_id = all_hit_id,
                   annot_df = annot, 
                   hit_id2node = hit_id2node)    

        print(f'Score done for annotation ::: {score_index + 1}/{num_score} : {round(((score_index + 1)/num_score)*100, 2)}%', end='\r', flush=True)
