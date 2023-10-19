#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "FAUCHOIS Antoine"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["FAUCHOIS Antoine"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "FAUCHOIS Antoine"
__email__ = "antoine.fauchois92@gmail.com"
__status__ = "Developpement"

def isfile(path): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file
    
    :raises ArgumentTypeError: If file doesn't exist
    
    :return: (str) Path 
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Extract reads from fastq files.

    :param fastq_file: (str) Path to the fastq file.
    :return: A generator object that iterate the read sequences. 
    """
    if not os.path.exists(fastq_file):
        sys.exit(f"ERROR: {fastq_file} doesn't exist")
    with open(fastq_file, "r") as input_file:
        for i, line in enumerate(input_file):
            #Sequence is the second line among a block of 4 lines
            if i % 4 == 1:
                yield line.strip()


def cut_kmer(read, kmer_size):
    """Cut read into kmers of size kmer_size.
    
    :param read: (str) Sequence of a read.
    :return: A generator object that iterate the kmers of of size kmer_size.
    """
    for i in range(len(read) - kmer_size + 1):
        kmer = read[i:i + kmer_size]
        yield kmer


def build_kmer_dict(fastq_file, kmer_size):
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    kmer_dict = {}
    #Get reads from fastq_file
    reads = read_fastq(fastq_file)
    for read in reads:
        kmers = cut_kmer(read, kmer_size)
        for kmer in kmers:
            if kmer not in kmer_dict.keys():
                kmer_dict[kmer] = 1
            else:
                kmer_dict[kmer] += 1
    return kmer_dict


def build_graph(kmer_dict):
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    #Create empty oriented graph
    graph = nx.DiGraph()
    for kmer in kmer_dict.keys():
        #Extract prefix, suffix and weight
        prefix = kmer[:-1]
        suffix = kmer[1:]
        weight = kmer_dict[kmer]
        #Create edge
        graph.add_edge(prefix, suffix, weight = weight)
    return graph



def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    for path in path_list:
        #Create a path to remove according to boolean arguments
        if delete_entry_node and delete_sink_node:
            to_remove = path
        elif not delete_entry_node and delete_sink_node:
            to_remove = path[1:]
        elif delete_entry_node and not delete_sink_node:
            to_remove = path[:-1]
        else:
            if len(path) >= 3:
                to_remove = path[1:-1]
            else:
                to_remove = path   
        graph.remove_nodes_from(to_remove)
    return graph


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    #Create empty list of path to remove
    path_to_remove = []
    #Select paths to remove
    ##Selection based on variation in path weight
    if statistics.stdev(weight_avg_list) > 0:
        for i in range(len(path_list)):
            if weight_avg_list[i] != max(weight_avg_list):
                path_to_remove.append(path_list[i])
    ##Selection on variation in length path            
    elif statistics.stdev(path_length) > 0:
        for i in range(len(path_list)):
            if path_length[i] != max(path_length):
                path_to_remove.append(path_list[i])
    ##Random choice
    else:
        index_to_store = random.randint(0, len(path_list) - 1)
        for i in range(len(path_list)):
            if i != index_to_store:
                path_to_remove.append(path_list[i])
    #Delete all path in path to remove
    graph = remove_paths(graph, path_to_remove, delete_entry_node, delete_sink_node)
    return graph

def path_average_weight(graph, path):
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph 
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    #Create enpty weight_avg_list and path_length list
    weight_avg_list = []
    path_length = []

    #Call all simple path between start and end node of a given bubble
    path_list = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    #Compute mean weigth and length for each one
    for path in path_list:
        weight_avg = path_average_weight(graph, path)
        length = len(path)
        weight_avg_list.append(weight_avg)
        path_length.append(length)
    #Call select_best_path function
    graph = select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False)
    return graph

def simplify_bubbles(graph):
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    bubble = False
    for node in graph.nodes():
        #Extract predecessors nodes for analysed node
        previous_nodes = list(graph.predecessors(node))
        if len(previous_nodes) > 0:
            #Looking for a common ancestror node for each possible couple 
            for i in range(len(previous_nodes)):
                for j in range(i + 1, len(previous_nodes)):
                    ancestor_node = nx.lowest_common_ancestor(graph, previous_nodes[i], previous_nodes[j])
                    if ancestor_node is not None:
                        bubble = True
                        break
        if bubble:
            break
    if bubble:
        graph = simplify_bubbles(solve_bubble(graph, ancestor_node, node))
    return graph


def solve_entry_tips(graph, starting_nodes):
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    tip = False
    #Looking for a ending node couple which have the same ancestor
    for i in range(len(starting_nodes)):
        for j in range(i + 1, len(starting_nodes)):
            ancestor_node = nx.lowest_common_ancestor(graph.reverse(), starting_nodes[i], starting_nodes[j])
            if ancestor_node is not None:
                tip = True
                break
        if tip:
            break
    if tip:
        #Extract path list
        path_list = [list(nx.all_simple_paths(graph, starting_nodes[i], ancestor_node))[0]]
        path_list.append(list(nx.all_simple_paths(graph, starting_nodes[j], ancestor_node))[0])
        #Compute weight average and length
        weight_avg_list = []
        path_length = []        
        for path in path_list:
            weight_avg = path_average_weight(graph, path)
            length = len(path)
            weight_avg_list.append(weight_avg)
            path_length.append(length)
        #Update graph and ending_nodes       
        graph = select_best_path(graph, path_list, path_length, weight_avg_list,
        delete_entry_node = True, delete_sink_node = False)
        starting_nodes = get_starting_nodes(graph)
        #Recall solve_out_tips
        graph = solve_entry_tips(graph, starting_nodes)
    return graph


def solve_out_tips(graph, ending_nodes):
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    tip = False
    #Looking for a ending node couple which have the same ancestor
    for i in range(len(ending_nodes)):
        for j in range(i + 1, len(ending_nodes)):
            ancestor_node = nx.lowest_common_ancestor(graph, ending_nodes[i], ending_nodes[j])
            if ancestor_node is not None:
                tip = True
                break
        if tip:
            break
    if tip:
        #Extract path list
        path_list = [list(nx.all_simple_paths(graph, ancestor_node, ending_nodes[i]))[0]]
        path_list.append(list(nx.all_simple_paths(graph, ancestor_node, ending_nodes[j]))[0])
        #Compute weight average and length
        weight_avg_list = []
        path_length = []        
        for path in path_list:
            weight_avg = path_average_weight(graph, path)
            length = len(path)
            weight_avg_list.append(weight_avg)
            path_length.append(length)
        #Update graph and ending_nodes       
        graph = select_best_path(graph, path_list, path_length, weight_avg_list,
        delete_entry_node=False, delete_sink_node = True)
        ending_nodes = get_sink_nodes(graph)
        #Recall solve_out_tips
        graph = solve_out_tips(graph, ending_nodes)
    return graph

def get_starting_nodes(graph):
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    #Create empty list of starting nodes
    starting_nodes = []
    #Looking for predecessrs for each node
    for node in graph.nodes():
        if len(list(graph.predecessors(node))) == 0:
            starting_nodes.append(node)
    return starting_nodes

def get_sink_nodes(graph):
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    #Create empty list of ending nodes
    ending_nodes = []
    #Looking for sucessors for each node
    for node in graph.nodes():
        if len(list(graph.successors(node))) == 0:
            ending_nodes.append(node)
    return ending_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object 
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    #Create empty contig list
    contigs = []
    #Test all possible couple of starting and ending nodes 
    for starting_node in starting_nodes:
        for ending_node in ending_nodes:
            #extract list of simple paths
            simple_paths = list(nx.all_simple_paths(graph, starting_node, ending_node))
            #Build contig for each simple paths
            for simple_path in simple_paths:
                #Add the first node
                contig = simple_path[0]
                #Add the last base for the other node
                for node in simple_path[1:]:
                    contig += node[-1]
                to_add = (contig, len(contig))
                contigs.append(to_add)
    return contigs

def save_contigs(contigs_list, output_file):
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (str) Path to the output file
    """
    with open(output_file, "w") as output_fasta:
        for i, contig in enumerate(contigs_list):
            output_fasta.write(f">contig_{i} len={contig[1]}")
            output_fasta.write("\n")
            output_fasta.write(textwrap.fill(contig[0], width=80))
            output_fasta.write("\n")
    


def draw_graph(graph, graphimg_file): # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (str) Path to the output file
    """                                   
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    #Build kmer dict from fastq file
    print("Build kmer dictionary...")
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    
    #Build graph
    print("Build graph...")
    graph = build_graph(kmer_dict)

    #Remove bubble form graph
    print("Remove bubbles...")
    graph = simplify_bubbles(graph)

    print("Remove tips...")

    #Remove entry and out tips
    starting_nodes = get_starting_nodes(graph)
    graph = solve_entry_tips(graph, starting_nodes)

    ending_nodes = get_sink_nodes(graph)
    graph = solve_out_tips(graph, ending_nodes)

    print("Get contig from clean graph...")
    #Re extract starting and ending nodes after cleaning
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)

    #Get contig from graph
    contigs = get_contigs(graph, starting_nodes, ending_nodes)

    #Save created contigs
    print("Saving contig...")
    save_contigs(contigs, args.output_file)
    #Draw graph
    #if args.graphimg_file:
    #    draw_graph(graph, args.graphimg_file)


if __name__ == '__main__': # pragma: no cover
    main()
