import networkx as nx
from canonization import test_canon
import os

for (root, dirs, files) in os.walk("cfi"):
    print(root, dirs, files)

    for f in files:

        with open(os.path.join(root, f), "r") as file:
            print(f"Opening file {f}")
            g = nx.Graph()
            edges = []
            file.readline()
            for string in file:
                temp = string[2:]
                nodes = temp.split(" ")
                g.add_edge(int(nodes[0]) - 1, int(nodes[1]) - 1)
            
            # print(g.edges)
            if len(g.nodes) < 300:
                print(test_canon(g))