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


def update_leaves(leaves, i, j, S, D):
    """ Insert new nodes and update properties of the old nodes """
    # add new nodes and children
    edge_weight_i = 0.5 * (D[i, j] + r(i, S, D) - r(j, S, D))
    edge_weight_j= D[i, j] - edge_weight_i

    if S[i] not in leaves:
        leaves[S[i]] = Phylo.BaseTree.Clade(branch_length=edge_weight_i, name=S[i])
    if S[j] not in leaves:
        leaves[S[j]] = Phylo.BaseTree.Clade(branch_length=edge_weight_j, name=S[j])

    combined = S[i] + S[j]
    leaves[combined] = Phylo.BaseTree.Clade(branch_length=None, name=combined, clades=[leaves[S[i]], leaves[S[j]]])
    return leaves


def calculate_edge_weight(i, j, m, D):
    """ Calculate the weight after the termination of the while loop """
    return (D[i, j] + D[i, m] - D[j, m]) / 2


def create_tree(leaves, S, D):
    """ Update tree after while loop has terminated """
    # Add nodes and children
    edge_weights = []
    edge_weights.append((D[0, 1] + D[0, 2] - D[1, 2]) / 2)
    edge_weights.append((D[0, 1] + D[1, 2] - D[0, 2]) / 2)
    edge_weights.append((D[0, 2] + D[1, 2] - D[0, 1]) / 2)

    for i in range(len(S)):
        if S[i] not in leaves:
            leaves[S[i]] = Phylo.BaseTree.Clade(branch_length=edge_weights[i], name=S[i])
        else:
            leaves[S[i]] = Phylo.BaseTree.Clade(branch_length=edge_weights[i],
                                                name=S[i],
                                                clades=leaves[S[i]].clades)

    combined = S[0]+S[1]+S[2]
    leaves[combined] = Phylo.BaseTree.Clade(branch_length=None,
                                                name=combined,
                                                clades=[leaves[S[0]], leaves[S[1]], leaves[S[2]]])
    tree = Phylo.BaseTree.Tree(root=leaves[combined])
    return tree


def neighbor_joining(D, taxa):
    """ Creates an unrooted binary tree such that similar species are grouped in the
        same sub tree and the tree distances correspond to the distance matrix if possible """
    S = copy.deepcopy(taxa)
    leaves = {}

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
        leaves = update_leaves(leaves, i, j, S, D)

        # 4. Update the dissimilarity matrix
        D = update_dissimilarity_matrix(S, D, i, j)

        # 5. Delete i and j from S and append the new taxon k
        S = update_S(S, i, j)

    # Add new node v and edges to the tree
    tree = create_tree(leaves, S, D)
    return tree


def main():
    # read in the dissimilarity matrix
    file = "example_slide4.phy"
    print("Reading from " + file + ".\n")
    D, taxa  = load_distance_matrix(file)

    # run the neighbor joining algorithm
    nj_tree = neighbor_joining(D, taxa)

    # draw the tree
    # Phylo.draw(nj_tree)
    print("--------------------------")
    print("ASCII")
    Phylo.draw_ascii(nj_tree)
    print("--------------------------\n")

    # print tree in newick format
    print("--------------------------")
    print("Newick")
    print(nj_tree.format("newick").strip())
    print("--------------------------\n")


if __name__ == '__main__':
    main()
