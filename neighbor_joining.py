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


def calculate_matrix_N(S, D):
    """ Calculate N as defined by Saitou and Nei """
    N = np.zeros(D.shape)
    n, _ = D.shape

    for i in range(n):
        for j in range(n):
            if i != j:
                N[i, j] = D[i, j] - (r(i, S, D) + r(j, S, D))

    return N


def find_neighbors(N):
    """ Find a pair of neighbors which minimized Nij as defined by Saitou and Nei """
    n, m = N.shape

    min = float('inf')
    pair = None
    for i in range(n):
        for j in range(n):
            if N[i, j] < min:
                min = N[i, j]
                pair = (i, j)

    return pair


def update_S(S, i, j):
    """ Update new list of taxa """
    taxon_i = S[i]
    taxon_j = S[j]

    # remove old nodes
    S.remove(taxon_i)
    S.remove(taxon_j)

    # add new taxon "k" = str(i) + str(j)
    S.append(taxon_i + taxon_j)
    return S


def update_dissimilarity_matrix(S, D, i, j):
    """ Update the dissamilarity matrix """
    # delete rows and columns for i and j
    D_reduced = np.delete(D, i, 0)
    D_reduced = np.delete(D_reduced, i, 1)
    D_reduced = np.delete(D_reduced, j - 1, 0) # reduce j by if i was already deleted in an index before j
    D_reduced = np.delete(D_reduced, j - 1, 1)

    n, _ = D_reduced.shape # dimensions without i and j
    D_new = np.zeros((n + 1, n + 1)) # create D_new that with one more row and column
    D_new[:n, :n] = D_reduced # fill D_new with D_reduced

    l, _ = D.shape # original dimensions
    new_row_column = []
    for m in range(l):
        if m != i and m != j:
            new_row_column.append(0.5 * (D[i, m] + D[j, m] - D[i, j]))

    new_row_column.append(0.0) # append 0.0 for new taxon in D
    D_new[-1, :], D_new[:, -1]= new_row_column, new_row_column # set new row and column

    return D_new


def r(i, S, D):
    """ Calculate r_i = 1 / (|S| - 2) Σ D_i,m """
    n, _ = D.shape
    return sum(D[i, m] for m in range(n))  / (len(S) - 2)


def tree_object(taxa):
    """ Save name, parent and edge_weight for each taxon """
    tree = {}
    for taxon in taxa:
        tree[taxon] = ["root", 0.0]

    return tree


def update_tree_object(tree_obj, i, j, S, D):
    """ Insert new nodes and update properties of the old nodes """
    n, _ = D.shape

    # add new taxon
    tree_obj[S[i] + S[j]] = ["root", 0.0]

    # update old taxons with new parent and edge weights
    edge_weight = 0.5 * (D[i, j] + r(i, S, D) - r(j, S, D))
    tree_obj[S[i]] = [S[i] + S[j], edge_weight]
    edge_weight = D[i, j] - edge_weight
    tree_obj[S[j]] = [S[i] + S[j], edge_weight]

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
    # return Phylo.read(io.StringIO(tree), "newick")
    return 0


def neighbor_joining(D, taxa):
    """ Creates an unrooted binary tree such that similar species are grouped in the
        same sub tree and the tree distances correspond to the distance matrix if possible """
    S = copy.deepcopy(taxa)
    tree_obj = tree_object(taxa)

    while len(S) > 3:
        # print information about the dissimilarity matrix
        print("--------------------------")
        print("Dissimilarity matrix")
        print(S)
        print(D)
        print("--------------------------\n")

        # 1. a) Compute matrix N
        N = calculate_matrix_N(S, D)
        print("--------------------------")
        print("Calculated N matrix")
        print(S)
        print(N)

        # 1. b) Select i, j such that it is a minimum entry in N
        i, j = find_neighbors(N)

        print("\nDeleting indices " + str(i) + " (" + S[i] + ") " + "and " + str(j)  + " (" + S[j] + ") " + "from D.")
        print("Creating new taxon " + "'" + S[i] + S[j] + "'.")
        print("--------------------------\n")

        # 2. Add new node k to the tree
        # 3. Add edges with weights to the tree
        tree_obj = update_tree_object(tree_obj, i, j, S, D)

        # 4. Update the dissimilarity matrix
        D = update_dissimilarity_matrix(S, D, i, j)

        # 5. Delete i and j from S and append the new taxon k
        S = update_S(S, i, j)

    # Add new node v and edges to the tree
    tree_obj = last_update_tree_obj(tree_obj, S, D)
    # print(tree_obj)
    return convert_to_newick_format(tree_obj)


def main():
    # read in the dissimilarity matrix
    file = "example_slide4.phy"
    print("Reading from " + file + ".\n")
    D, taxa  = load_distance_matrix(file)

    # run the neighbor joining algorithm
    nj_tree = neighbor_joining(D, taxa)

    # draw the tree
    # Phylo.draw(tree)
    # Phylo.draw_ascii(tree)


if __name__ == '__main__':
    main()
