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


def tree_object(taxa):
    """ Save name, parent and edge_weight for each taxon """
    tree = {}
    for taxon in taxa:
        tree[taxon] = ["root", 0.0]

    return tree


def update_tree_object(tree_obj, taxa, i, j, S, D):
    """ Insert new nodes and update properties of the old nodes """
    n, _ = D.shape

    print(taxa[i], taxa[j])

    # add new taxon
    tree_obj[taxa[i] + taxa[j]] = ["root", 0.0]

    # update old taxons with new parent and edge weights
    edge_weight = 1/2 * (D[i, j] + 1 / (len(S) - 2) * (1 / (len(S) - 2)) * sum(D[i, k] for k in range(n)))
    tree_obj[taxa[i]] = [taxa[i] + taxa[j], edge_weight]
    edge_weight = D[i, j] - edge_weight
    tree_obj[taxa[j]] = [taxa[i] + taxa[j], edge_weight]

    return tree_obj


def calculate_edge_weight(i, j, m, D):
    """ Calculate the weight after the termination of the while loop """
    return (D[i, j] + D[i, m] - D[j, m]) / 2


def last_update_tree_obj(tree_obj, S, D):
    """ Update tree after while loop has terminated """
    # Add new node v
    v = ""
    for taxon in S:
        v += taxon
    tree_obj[v] = ["root", 0.0]

    # Add new weights
    edge_weight = D[0, 1] + D[0, 2] - D[1, 2]
    tree_obj[S[0]] = [v, edge_weight]

    edge_weight = D[0, 1] + D[1, 2] - D[0, 2]
    tree_obj[S[1]] = [v, edge_weight]

    edge_weight = D[0, 2] + D[1, 2] - D[0, 1]
    tree_obj[S[2]] = [v, edge_weight]

    return tree_obj


def convert_to_newick_format(tree):
    """ Convert the tree object to newick format """
    return 0


def neighbor_joining(D, taxa):
    """ Creates an unrooted binary tree such that similar species are grouped in the
        same sub tree and the tree distances correspond to the distance matrix if possible """
    S = copy.deepcopy(taxa)
    tree_obj = tree_object(taxa)
    # Phylo.draw(tree)
    # Phylo.draw_ascii(tree)

    while len(S) > 3:
        # 1. a) Compute matrix N
        N = calculate_matrix_N(D, S)

        # 1. b) Select i, j such that it is a minimum entry in N
        i, j = find_neighbors(N)

        # 2. Add new node k to the tree
        # 3. Add edges with weights to the tree
        tree_obj = update_tree_object(tree_obj, taxa, i, j, S, D)

        # 4. Update the dissimilarity matrix
        D = update_dissimilarity_matrix(D, i, j)

        # 5. Delete i and j from S and append the new taxon k
        S = update_S(S, i, j)

    # Add new node v and edges to the tree
    tree_obj = last_update_tree_obj(tree_obj, S, D)
    return convert_to_newick_format(tree_obj)


def main():
    D, taxa  = load_distance_matrix("example_slide4.phy")

    nj = neighbor_joining(D, taxa)



if __name__ == '__main__':
    main()
