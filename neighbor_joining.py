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
                        row.append(entry)
                D.append(row)
            count += 1

        return np.array(D), nodes


def calculate_Q(D):
    """ Calculate the Q-matrix based on our distance matrix D """
    pass


def neighbor_joining(D):
    """ Creates an unrooted binary tree such that similar species are grouped in the
        same sub tree and the tree distances correspond to the distance matrix if possible """
    pass


def main():
    D, nodes  = load_distance_matrix("example_slide4.phy")



if __name__ == '__main__':
    main()
