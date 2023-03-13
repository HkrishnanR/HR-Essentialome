"""
Created April 2021
@author: Harikrishnan Ramadasan <harikrishnan@students.iisertirupati.ac.in>
- HR Essentialome analysis
"""

import random
from collections import Counter
import pandas as pd
import networkx as nx
import numpy as np
import os
import sys
import time

current = time.time()


def get_degree_preserving_randomization(edges):
    """

    :param edges: list with each tuple stored as an edge
    :return: a list of new edges
    """
    edges = set([tuple(e) for e in edges]) # create a tuple of all edges

    degrees = []
    [degrees.extend(e) for e in edges]  # Stores degrees in a list

    degree_counter = Counter(degrees) #counts th

    new_edges = set()

    nodes = np.array(
        [degree for degree, count in degree_counter.items() if count != 0])

    while len(nodes) > 0:

        first, second = -1, -1

        while first == second and len(nodes) > 1:
            first, second = np.random.choice(nodes, size=(2,), replace=False)  # randomly select two nodes

        if first != second and \
                (first, second) not in new_edges and \
                (second, first) not in new_edges and \
                len(nodes) > 1:
            new_edges.add(
                (first, second))  # Create an edge with the random node pair and subtract one from degrees of each node
            degree_counter[first] -= 1
            degree_counter[second] -= 1
        else:
            edge = random.sample(new_edges, 1)[0]  # remove a random edge and update the degrees
            new_edges.remove(edge)
            degree_counter[edge[0]] += 1
            degree_counter[edge[1]] += 1

        nodes = np.array(
            [degree for degree, count in degree_counter.items() if count != 0])

    return list(new_edges)


df = pd.read_csv(sys.argv[1], sep="\t")  # First system argument is the file

os.makedirs("Randomized", exist_ok=True)

for i in range(int(sys.argv[2])):  # Second system argument is the number of randomisations
    G = nx.from_pandas_edgelist(df, 'UniprotA', 'UniprotB')

    edges = [i for i in G.edges()]
    print("Reading Edgelist..........")
    print("Done \n")

    print("Randomizing....")

    li = get_degree_preserving_randomization(edges)

    print("Randomized {} time".format(i + 1))

    dff = pd.DataFrame(li, columns=["UniprotA", "UniprotB"])
    dff.to_csv("Randomized/Human_ppi_ran_" + str(i + 1) + ".txt", sep="\t")

    print("Done \n \n")

print(time.time() - current)
# print(dff)
