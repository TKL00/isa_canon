import os
from canonization import permute_edges, graph_canon, first_non_trivial, rightmost_first_non_trivial, create_adjacency
import networkx as nx
import numpy as np 
import random
import time
import timeout_decorator

CANON_TIMEOUT = 600

def permuted_graph(G):
    n_nodes = len(G.nodes)
    permutation = np.random.permutation(n_nodes)
    permutation_mapping = {i: permutation[i] for i in range(n_nodes)}
    permuted_edges = permute_edges(permutation_mapping, G.edges)
    permuted_graph = nx.Graph()
    permuted_graph.add_nodes_from(G.nodes)
    permuted_graph.add_edges_from(permuted_edges)

    return permuted_graph


@timeout_decorator.timeout(CANON_TIMEOUT)
def call_graph_canon(G, Q, T):
    return graph_canon(G, Q, T)


def table1(target_cell_selector, traces):
    """
    n nodes | m edges | time G_1 | time G_2 | automorphisms found G_1 | automorphisms found G_2 | trace impact 
    """

    for (root, dirs, files) in os.walk("random_graphs"):
        
        for file_name in sorted(files):
            full_path = os.path.join(root, file_name)

            test_graph = nx.read_adjlist(full_path, nodetype=int)
            n_nodes = len(test_graph.nodes)
            edges = len(test_graph.edges)
            
            permuted_test_graph = permuted_graph(test_graph)

            time_before = time.time()
            labeling_1, automorphisms_1, trace_impact_1 = graph_canon(test_graph, target_cell_selector, traces)
            time_after = time.time()
            time_1 = time_after - time_before
            
            time_before = time.time()
            labeling_2, automorphisms_2, trace_impact_2 = graph_canon(permuted_test_graph, target_cell_selector, traces)
            time_after = time.time()
            time_2 = time_after - time_before
            
            if traces:
                print(f"{n_nodes}\t{edges}\t{time_1}\t{time_2}\t{len(automorphisms_1)}\t{len(automorphisms_2)}\t{trace_impact_1}\t{trace_impact_2}")
            else:
                print(f"{n_nodes}\t{edges}\t{time_1}\t{time_2}\t{len(automorphisms_1)}\t{len(automorphisms_2)}")

            ## CHECK FOR EQUALITY
            # G_canon_edge = permute_edges(labeling_1, test_graph.edges)
            # perm_G_canon_edge = permute_edges(labeling_2, permuted_test_graph.edges.edges)
            # G_canon_adj = create_adjacency(node_class, G_canon_edge)
            # perm_G_canon_adj = create_adjacency(node_class, perm_G_canon_edge)
            # np.array_equal(G_canon_adj, perm_G_canon_adj)


def table2(dirpath, target_cell_selector, traces):
    """
        filename | n nodes | m edges | time(s) | automorphism found | trace impact
    """
    for (root, dirs, files) in os.walk(dirpath):

        for f in files:

            with open(os.path.join(root, f), "r") as file:
                g = nx.Graph()
                file.readline()
                for string in file:
                    temp = string[2:]
                    nodes = temp.split(" ")
                    g.add_edge(int(nodes[0]) - 1, int(nodes[1]) - 1)
                
                # print(g.edges)
                if len(g.nodes) < 250:
                    try:
                        time_before = time.time()
                        labeling, automorphisms, trace_impact = call_graph_canon(g, target_cell_selector, traces)
                        time_after = time.time()
                        time_spent = time_after - time_before

                        if traces:
                            print(f"{f}\t{len(g.nodes)}\t{len(g.edges)}\t{time_spent}\t{len(automorphisms)}\t{trace_impact}")
                        else:
                            print(f"{f}\t{len(g.nodes)}\t{len(g.edges)}\t{time_spent}\t{len(automorphisms)}")

                    except: 
                        print(f"{f}\t{len(g.nodes)}\t{len(g.edges)}\ttimeout")

def table_folder(dir_name):
    print(f"The target cell: first_non_trivial and traces: false")
    print(f"graph_name\tnode_count\tedge_count\ttime(s)\taut")
    table2(dir_name, first_non_trivial, False)

    print(f"The target cell: first_non_trivial and traces: true")
    print(f"graph_name\tnode_count\tedge_count\ttime(s)\taut\ttrace_impact")
    table2(dir_name, first_non_trivial, True)


    print(f"The target cell: rightmost_first_non_trivial and traces: false")
    print(f"graph_name\tnode_count\tedge_count\ttime(s)\taut")
    table2(dir_name, rightmost_first_non_trivial, False)

    print(f"The target cell: rightmost_first_non_trivial and traces: true")
    print(f"graph_name\tnode_count\tedge_count\ttime(s)\taut\ttrace_impact")
    table2(dir_name, rightmost_first_non_trivial, True)


if __name__ == "__main__":
    print(f"The target cell: first_non_trivial and traces: false")
    print(f"node_count\tedge_count\ttime_1(s)\ttime_2(s)\taut_1\taut_2")
    table1(first_non_trivial, False)
    
    print(f"The target cell: first_non_trivial and traces: true")
    print(f"node_count\tedge_count\ttime_1(s)\ttime_2(s)\taut_1\taut_2\ttrace_imp_1\ttrace_imp_2")
    table1(first_non_trivial, True)

    print(f"The target cell: rightmost_first_non_trivial and traces: false")
    print(f"node_count\tedge_count\ttime_1(s)\ttime_2(s)\taut_1\taut_2")
    table1(rightmost_first_non_trivial, False)

    print(f"The target cell: rightmost_first_non_trivial and traces: true")
    print(f"node_count\tedge_count\ttime_1(s)\ttime_2(s)\taut_1\taut_2\ttrace_imp_1\ttrace_imp_2")
    table1(rightmost_first_non_trivial, True)

    # table_folder("mz-aug2")

    # table_folder("cfi")

    # table_folder("ag")


    

    