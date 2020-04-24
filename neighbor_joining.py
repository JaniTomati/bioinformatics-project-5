#!/usr/bin/python3
# -*- coding: utf-8 -*-

import io
import numpy as np
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


def calculate_Q(D):
    """ Calculate the Q-matrix based on our distance matrix D """
    Q = np.zeros(D.shape)
    n, m = D.shape

    for i in range(n):
        for j in range(m):
            Q[i, j] = (n - 2) * D[i, j] - sum(D[i, k] for k in range(n)) - sum(D[j, k] for k in range(n))

    return Q


def find_minimum_Q(Q):
    """ Find the pair of nodes for which Q is minimized """
    n, m = Q.shape

    min = float('inf')
    pair = None
    for i in range(n):
        for j in range(m):
            if Q[i, j] < min:
                min = Q[i, j]
                pair = (i, j)

    return pair


def select_nearest_pair(D):
    n, m = D.shape

    min = float('inf')
    pair = None
    for i in range(n):
        for j in range(m):
            if D[i, j] < min and i != j:
                min = D[i, j]
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


def generic_neighbor_joining(D, taxa):
    """ Creates an unrooted binary tree such that similar species are grouped in the
        same sub tree and the tree distances correspond to the distance matrix if possible """
    S = taxa
    Q = calculate_Q(D)
    print(select_nearest_pair(D))

    print(create_tree(taxa))
    Phylo.draw(create_tree(taxa))

    while len(S) > 3:
        break

    return 0


def main():
    D, taxa  = load_distance_matrix("example_slide4.phy")

    nj = generic_neighbor_joining(D, taxa)



if __name__ == '__main__':
    main()
