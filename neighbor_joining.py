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


def calculate_matrix_N(D):
    """ Calculate N as defined by Saitou and Nei """
    N = np.zeros(D.shape)
    n, m = D.shape

    for i in range(n):
        for j in range(m):
            if i != j:
                N[i, j] = D[i, j] - (1 / (n - 2) * sum(D[i, k] for k in range(n)) + 1 / (n - 2) * sum(D[j, k] for k in range(n)))

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


def create_tree(T):
    """ Create a tree in newick format using a list of taxa """
    tree = "("
    for taxon in T:
        if taxon != T[-1]:
            tree += taxon + ", "
        else:
            tree += taxon + ")"

    return Phylo.read(io.StringIO(tree), "newick")


def update_dissimilarity_matrix(D, i, j):
    """ Update the dissamilarity matrix """
    D_ = np.delete(D, i, 0)
    D_ = np.delete(D_, i, 1)
    D_ = np.delete(D_, j - 1, 0) # -1 since we already deleted i from the matrix
    D_ = np.delete(D_, j - 1, 1)

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


def update_taxa(S, i, j):
    """ Update new list of taxa """
    taxon_i = S[i]
    taxon_j = S[j]

    # remove old nodes
    S.remove(taxon_i)
    S.remove(taxon_j)

    # add new taxon "k"
    S.append(taxon_i + taxon_j)

    return S


def neighbor_joining(D, taxa):
    """ Creates an unrooted binary tree such that similar species are grouped in the
        same sub tree and the tree distances correspond to the distance matrix if possible """
    S = copy.deepcopy(taxa)
    # tree = create_tree(taxa)
    # Phylo.draw(tree)
    # Phylo.draw_ascii(tree)

    while len(S) > 3:
        N = calculate_matrix_N(D)
        i, j = find_neighbors(N)

        S = update_taxa(S, i, j)
        print(S)
        D = update_dissimilarity_matrix(D, i, j)

    print(D)

    return 0


def main():
    D, taxa  = load_distance_matrix("example_slide4.phy")

    nj = neighbor_joining(D, taxa)



if __name__ == '__main__':
    main()
