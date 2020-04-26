#!/usr/bin/python3
# -*- coding: utf-8 -*-

import io
import copy
import numpy as np
import networkx as nx
from Bio import Phylo


def load_distance_matrix(file):
    """ Load a distance matrix from a file """
    D = []
    taxa = []

    with open(file, "r") as f:
        lines = f.readlines()

        count = 0
        for line in lines:
            line_arr = line.rstrip().split()
            if count != 0:
                row = []
                for entry in line_arr:
                    if entry == line_arr[0]: # add taxon names
                        taxa.append(entry)
                    else:                    # add distance
                        row.append(float(entry))
                D.append(row)
            count += 1

        return np.array(D), taxa


def calculate_matrix_N(D, S):
    """ Calculate N as defined by Saitou and Nei """
    N = np.zeros(D.shape)
    n, m = D.shape

    for i in range(n):
        for j in range(m):
            if i != j:
                N[i, j] = D[i, j] - (1 / (len(S) - 2) * sum(D[i, k] for k in range(n)) + 1 / (len(S) - 2) * sum(D[j, k] for k in range(n)))

    return N


def find_neighbors(N):
    """ Find a pair of neighbors which minimized Nij as defined by Saitou and Nei """
    n, m = N.shape

    min = float('inf')
    pair = None
    for i in range(n):
        for j in range(m):
            if N[i, j] < min and i != j:
                min = N[i, j]
                pair = (i, j)

    return pair


def create_newick_tree(T):
    """ Create a tree in newick format from the initial list of taxa """
    tree = "("
    for taxon in T:
        if taxon != T[-1]:
            tree += taxon + ", "
        else:
            tree += taxon + ")"

    return tree
    # return Phylo.read(io.StringIO(tree), "newick")


def create_tree_graph(taxa):
    """ Create a graph that represents our tree """
    tree = nx.Graph()
    tree.add_node("root") # insert a dummy root node which will be deleted later
    tree.add_nodes_from(taxa)

    for taxon in taxa:
        tree.add_edge("root", taxon) # connect each taxon to the root

    return tree


def update_dissimilarity_matrix(D, i, j):
    """ Update the dissamilarity matrix """
    shift = -1
    if i > j:
        shift = 0

    D_ = np.delete(D, i, 0)
    D_ = np.delete(D_, i, 1)
    D_ = np.delete(D_, j + shift, 0) # shift j one to the left if i was already deleted
    D_ = np.delete(D_, j + shift, 1)

    n, _ = D_.shape
    new = []
    for m in range(n):
        new.append(1 / 2 * (D[i, m] + D[j, m] - D[i, j]))

    new.append(0.0)

    D = np.zeros((n + 1, n + 1))
    D[:n, :n] = D_
    D[-1, :] = np.array(new)
    D[:, -1] = np.array(new)

    return D


def update_S(S, i, j):
    """ Update new list of taxa """
    taxon_i = S[i]
    taxon_j = S[j]

    # remove old nodes
    S.remove(taxon_i)
    S.remove(taxon_j)

    # add new taxon "k"
    S.append(taxon_i + taxon_j)

    return S


def update_tree(tree, taxa, i, j, D):
    """ Delete old nodes  and insert new nodes into the tree model """
    print(tree.nodes())
    tree.remove_node(taxa[i])
    tree.remove_node(taxa[j])

    tree.add_node(taxa[-1])
    tree.add_node(taxa[i])
    tree.add_node(taxa[j])

    # weight_i = 1 / 2 * (D[i, j] + )
    # weight_j =



def neighbor_joining(D, taxa):
    """ Creates an unrooted binary tree such that similar species are grouped in the
        same sub tree and the tree distances correspond to the distance matrix if possible """
    S = copy.deepcopy(taxa)
    tree = create_tree_graph(taxa) # create a newick format from the given taxa
    # Phylo.draw(tree)
    # Phylo.draw_ascii(tree)

    while len(S) > 3:
        # 1. a) Compute matrix N
        N = calculate_matrix_N(D, S)

        # 1. b) Select i, j such that it is a minimum entry in N
        i, j = find_neighbors(N)

        # 2. Add node k to the tree

        # 3. Add edges with weights to the tree

        # 4. Update the dissimilarity matrix
        D = update_dissimilarity_matrix(D, i, j)

        # 5. Delete i and j from S and append the new taxon k
        S = update_S(S, i, j)

    # Add new node v


    # Add edges to the tree

    return 0


def main():
    D, taxa  = load_distance_matrix("example_slide4.phy")

    nj = neighbor_joining(D, taxa)



if __name__ == '__main__':
    main()
