#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np


def load_distance_matrix(file):
    """ Load a distance matrix from a file """
    D = []
    nodes = []

    with open(file, "r") as f:
        lines = f.readlines()

        count = 0
        for line in lines:
            line_arr = line.rstrip().split()
            if count != 0:
                row = []
                for entry in line_arr:
                    if entry == line_arr[0]: # add node names
                        nodes.append(entry)
                    else:                    # add distance
                        row.append(float(entry))
                D.append(row)
            count += 1

        return np.array(D), nodes


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


def neighbor_joining(D):
    """ Creates an unrooted binary tree such that similar species are grouped in the
        same sub tree and the tree distances correspond to the distance matrix if possible """

    Q = calculate_Q(D)
    pair = find_minimum_Q(Q)
    
    return 0


def main():
    D, nodes  = load_distance_matrix("example_slide4.phy")

    nj = neighbor_joining(D)



if __name__ == '__main__':
    main()
