import networkx as nx
from TreeNode import TreeNode
import copy
import numpy as np
from itertools import permutations 


def permute_edges(mapping, edge_set):
    new_edge_list = []
    for edges in edge_set:
        (u, v) = edges
        new_edge_list.append((mapping[u], mapping[v]))

    return new_edge_list

def create_adjacency(number_of_nodes, edge_set):
    adj = np.zeros((number_of_nodes, number_of_nodes))

    for (u, v) in edge_set:
        adj[u][v] = adj[v][u] = 1
    return adj


def graph_canon(G):
    """
        Produces the canonical labeling of the Graph G using the Graph Canonicalization 
        algorithm suggested by Brendan McKay.

        `Paramters`:
            G (Graph): NetworkX Graph

        `Returns`:
            Canonical_labeling (dict: node -> node): Relabeling mapping representing G in its canonical representation
            Automorphisms (list( dict: node -> node )): The generating set for the automorphism group of G.
    """

    G_NODE_AMT = len(G.nodes)
    
    def equitable_refinement(p):
        """
            Produces an equitable ordered partition based on the input partition 'p' of the graph 'G'
            using the algorithm suggested by Brendan McKay.

            `Parameters`
                p (list( list(int) )): List of partitioned nodes (in lists).
        """
        def shatter_set(p, adj):
            """
                Computes the list of tuples (i, j) where partition 'j' in tau shatters partition 'i'.

                `Paramters`:
                    p (list( list(int) )): List of partitioned nodes (in lists).
                    
                    adj: The adjacency matrix.

                `Returns`:
                    B (list( (int, int) )): A list of tuples, each (i, j) indicating that set V_i can be shattered by V_j.
            """
            B = []
            ## All parts in 'p'
            for i in range(len(p)):
                V_i = p[i]
                ## Ignore trivial parts (containing only one vertex)
                if len(V_i) > 1:

                    for j in range(len(p)):
                        V_j = p[j]

                        ## Reset limit for each chosen V_j
                        target_degree = 0

                        ## Counting degrees for all u in V_i to nodes v in V_j
                        for nodes in range(len(V_i)):
                            current_degree = 0
                            for target in range(len(V_j)):
                                u = V_i[nodes]
                                v = V_j[target]
                                if adj[u][v] == 1:
                                    current_degree += 1
                                
                            ## First vertex stores the target degree allowed for the verticies in V_i towards V_j
                            if nodes == 0:
                                target_degree = current_degree
                            elif current_degree != target_degree:
                                B.append( (i, j) )
                                ## Ignore rest of nodes in V_j and continue to V_(j+1) if one exists
                                break
            return B
        
        def shatter_partition(p, adj, to_shatter, shatter_with):
            """
                Computes a new partition from 'p' resulting from shattering the 
                'to_shatter' part with the 'shatter_with' part.

                `Parameters`:
                    p (list (list(int))): Current list of partitions (node lists).

                    adj: The adjacency matrix.

                    to_shatter: Index of part in `p` to shatter.

                    shatter_with: Index of part in `p` to shatter with.

                `Returns`:
                    A new list corresponding to the adjusted partition of 'p', 
                    having replaced `to_shatter` with its shattered parts.
            """
            V_i = p[to_shatter]
            V_j = p[shatter_with]
            degree_map = {}

            for shatter_node in V_i:
                degree = 0
                for node in V_j:
                    if adj[shatter_node][node] == 1:
                        degree += 1
                ## Node in shatter set must have 'degree' edges into the shattered with set                
                degree_map[shatter_node] = degree
            
            ## get list of (nodes, degree) and sort it on degree (ascending order)
            node_degree = sorted([ (v, degree_map[v]) for v in degree_map], key=lambda tuple: tuple[1])
            
            ## track the first node degree in the list
            track_degree = node_degree[0][1]

            ## list of shattered partitions of V_i
            shatter_partitions = [[]]
            current_partition = 0
            for values in node_degree:
                (node, degree) = values
                ## Current part is done, track new degree and begin new part.
                if degree != track_degree:
                    track_degree = degree
                    shatter_partitions.append([])
                    current_partition += 1
                shatter_partitions[current_partition].append(node)
            
            ## Replace partition V_i in 'partition' with all partitions in 'shatter_partitions'
            return_partition = []
            for i in range(len(p)):
                if i != to_shatter:
                    return_partition.append(p[i])
                else:
                    for j in range(len(shatter_partitions)):
                        return_partition.append(shatter_partitions[j])
            
            return return_partition

        adj_matrix = create_adjacency(G_NODE_AMT, G.edges)

        tau = p
        
        B = shatter_set(p, adj_matrix)

        ## NOTE: Abort refinement if current partition is lexicographically larger than the global minimum
        while len(B) != 0:
            ## Choose individualizing the lexicographically smallest
            (i, j) = sorted(B)[0]
            tau = shatter_partition(tau, adj_matrix, i, j)
            ## Update B to see if partitions remains to be shattered
            B = shatter_set(tau, adj_matrix)

        return tau
   
    def generate_tree(root):

        """
            Generates the search tree for the given root partition. Returns the global minimum partition
            corresponding to the labeling of the graph G.
        """

        def individualize(p, v):
            """Individualizes the node 'v' in partition 'p' such that if 'v' in V_i:
            
            p = (V_1, V_2, ..., V_i, ..., V_n)
            => p indv. v = (V_1, V_2, ..., {v}, V_i \ {v}, ..., V_n)

            `Returns`:
                Returns a new list corresponding to p indv. v
            
            """
            new_p = []
            for part in p:
                if v in part:
                    ## Create and insert {v}
                    new_p.append([v])
                    new_part = []
                    ## Create V_i \ {v}
                    for nodes in part:
                        if nodes != v:
                            new_part.append(nodes)
                    ## Insert V_i \ {v}
                    new_p.append(new_part)
                ## All parts not containing v stays the same in the same relative position in the partition
                else:
                    new_p.append(part)
            
            return new_p

        def update_automorphisms(new_leaf_partition, new_leaf_adj, automorphisms):            
            ## If partitions are isomorphic, compute automorphism between them
            if np.array_equal(new_leaf_adj, global_minimum[1]):
                best_partition = global_minimum[0]
                ## pi_one goes from the canonical labelling to the input graph, the inverse should go from the input graph to the canonical labelling
                pi_one_inverse = {new_leaf_partition[i][0]: i for i in range(G_NODE_AMT)}
                pi_two = {i: best_partition[i][0] for i in range(G_NODE_AMT)}
                new_automorphism_math = {i: pi_two[pi_one_inverse[i]] for i in range(G_NODE_AMT)}

                ## NOTE: incorrect automorphism right now. Delete below appending.
                new_automorphism = {best_partition[i][0]: new_leaf_partition[i][0] for i in range(G_NODE_AMT)}
                print(f"New automorphism discovered between")
                print(new_leaf_partition)
                print(best_partition)
                print("Automorphism")
                print(new_automorphism_math)
                automorphisms.append(new_automorphism_math)

        def all_fixed(traverse_sequence, automorphism):
                """
                    Given a sequence of individualized nodes and an automorphism, returns true if all
                    nodes in the sequence are mapped to themselves (fixed) in the automorphism. Otherwise,
                    return false.
                """
                all_fixed = True
                for node in traverse_sequence:
                    if automorphism[node] != node: return False

                return all_fixed

        def calculate_orbit(children_list, automorphisms):
            """
                Calculates the orbits for each node child in the children list under closure of the automorphism,
                thereby partitioning the children list, where each orbit is a part.

                `Paramters`:

                    children_list (list(int)): List of children to possibly travers

                    automorphism (list( dict(int -> int) )): List of automorphisms represented as dictionaries (mapping one node to another)

                `Returns`:
                    all_orbits (list( list(int))): Partition of children_list corresponding to the orbits of the children.
            """

            all_orbits = []
            
            for child in children_list:
                member_of_orbit = False
                ## check for occurrence in other orbits. If so, no new orbit should be calculated.
                for orbit in all_orbits:
                    if child in orbit:
                        member_of_orbit = True
                        break
                ## calculate non-orbitted child's orbit
                if not member_of_orbit:
                    ## orbit is a set, accounting for duplicates
                    child_orbit = {child}
                    for automorphism in automorphisms:
                        child_orbit.add(automorphism[child])
                    ## Sort orbit as set-usage may interfere with ordering
                    all_orbits.append(sorted(list(child_orbit)))
            
            return all_orbits

        def generate_subtree(parent_node, partition, current_seq, to_indiv, automorphisms, global_minimum):
            """
                Subroutine to generate the rest of the tree spanning from this node

                `Parameters`:

                    parent_node (TreeNode): A TreeNode representing the parent node in the tree of the node being generated.

                    partition (list (list(int))): The current partition of the nodes in the graph.
                    
                    current_seq (list(int)): The current sequence of individualization that led to (including) this TreeNode's partition.

                    to_indiv (int): The element in the partition being individualized.

            """
            def array_less_than(x, y):
                """Returns true if x < y, false if x >= y"""
                for i in range(G_NODE_AMT):
                    for j in range(G_NODE_AMT):
                        if x[i][j] != y[i][j]:
                            if x[i][j] < y[i][j]: return True
                            else: return False
                ## equal
                return False

            indiv_partition = individualize(partition, to_indiv)
            refinement = equitable_refinement(indiv_partition)
            
            ## If parent is root node, the current_seq is empty and will only contain this individualized node
            if not current_seq:
                new_sequence = [to_indiv]
            else:
                new_sequence = copy.deepcopy(current_seq)
                new_sequence.append(to_indiv)

            this_node = TreeNode(refinement, parent_node, new_sequence)

            ## Find first non-trivial part of the refined partition
            ## NOTE: Input parameter-based target cell selector function here instead.
            children_list = []
            for i in range(len(refinement)):
                if len(refinement[i]) != 1:
                    children_list = copy.deepcopy(refinement[i])
                    break
            ## If no choices available, this is a leaf node
            if len(children_list) == 0:
                
                ## create mapping from current discrete partition. Note that partition[i] says that the node (= partition[i]) of the input graph 
                ## goes to color "i" i.e. node "i" in the canonical labelling.
                relabeling_map = {this_node.get_partition()[i][0]: i for i in range(G_NODE_AMT)}
                relabeled_edges = permute_edges(relabeling_map, G.edges)
                leaf_adj = create_adjacency(G_NODE_AMT, relabeled_edges)

                ## First encountered leaf is the global minimum
                if not global_minimum[0]:
                    global_minimum[0] = this_node.get_partition()
                    global_minimum[1] = leaf_adj
                
                new_leaf_partition = this_node.get_partition()
                print(f"\t \t \t Found leaf: {new_leaf_partition}")
                print(f"\t \t \t Permuting edges to {relabeled_edges}")
                ## Check if automorph with current best leaf. 
                update_automorphisms(new_leaf_partition, leaf_adj, automorphisms)
                
                ## Check if current relabeling is better than previous relabeling. Resulting array contains true if one elemen
                if array_less_than(leaf_adj, global_minimum[1]):
                    global_minimum[0] = this_node.get_partition()
                    global_minimum[1] = leaf_adj
                ## Check for automorphism between other leaf partitions
            else:
                this_node.set_children(children_list)
                travers_seq = this_node.get_travers_seq()

                ## for each child visited, calculate new orbits and decide whether remaining children
                ## should be visisted
                for child in this_node.get_children():
                    ## orbit containing this child
                    A_prime = list(filter(lambda auto: all_fixed(travers_seq, auto), automorphisms))
                    orbits = calculate_orbit(children_list, A_prime)
                    ## Filter the partitioned children list to retrieve the child's orbit (list containing one orbit)
                    child_orbit = list(filter(lambda orbit: child in orbit, orbits))[0]
                    ## Only the first element of the orbit should be discovered. Since everything is lexicographically sorted, a child
                    ## that appears late in an orbit should not be traversed.
                    if child_orbit[0] == child: generate_subtree(this_node, refinement, new_sequence, child, automorphisms, global_minimum)

                                             ## GENERATE TREE ROOT
        
        ## list of automorphism (dicts)
        automorphisms = []
        ## the global minimum consist of 1) a partition and 2) the resulting adjacency matrix when permuting G.
        global_minimum = [[], []]

        ## Generate the root node's canonical

        print(f"Root children {root.get_children()}")
        for child in root.get_children():
            ## NOTE: Check for pruning here using all automorphisms
            orbits = calculate_orbit(root.get_children(), automorphisms)
            print(f"All orbits: {orbits}")
            child_orbit = list(filter(lambda orbit: child in orbit, orbits))[0]
            print(f"Orbits of {child}: {child_orbit}")
            if child_orbit[0] == child: generate_subtree(root, root.get_partition(), [], child, automorphisms, global_minimum)

        return global_minimum[0], automorphisms
    
    
    ##                                  CANONICAL
    root_partition = [sorted(list(G.nodes))]
    init_refinement = equitable_refinement(root_partition)
    print(f"INIT REFINEMENT {init_refinement}")   
    root_node = TreeNode(init_refinement, None, [])

    ## Find first non-trivial part of the refined partition
    children_list = []
    for i in range(len(init_refinement)):
        if len(init_refinement[i]) != 1:
            children_list = copy.deepcopy(init_refinement[i])
            break

    root_node.set_children(children_list)

    canonical_partition, automorphisms = generate_tree(root_node)
    canonical_labeling = {i:canonical_partition[i][0] for i in range(G_NODE_AMT)}

    return canonical_labeling, automorphisms



if __name__ == "__main__":

    graph = nx.Graph()
    graph.add_nodes_from([i for i in range(9)])
    graph.add_edges_from([(0, 1), (0, 3), (1, 2), (1, 4), (2, 5), (3, 4), (3, 6), (4, 5), (4, 7), (5, 8), (6, 7), (7, 8)])

    print("FIRST RUN:")
    labeling, automorphisms = graph_canon(graph)
    new_edges = permute_edges(labeling, graph.edges)
    print()
    print()

    

    canonical_graph = nx.Graph()
    canonical_graph.add_edges_from(new_edges)

    # all_mappings = []
    # l = list(permutations(range(9)))
    # for perm in l:
    #     new_mapping = {}
    #     for i in range(9):
    #         new_mapping[i] = perm[i]
    #     all_mappings.append(new_mapping)

    # canon_adj = create_adjacency(9, new_edges)
    # # print("Original canon adjacency:")
    # # print(canon_adj)

    # iterator = 0
    # for mapping in all_mappings:
    #     print(f"RUN {iterator}")
    #     print(f"Permutation: {mapping}")
    #     new_edge_set = permute_edges(mapping, graph.edges)
    #     relabeled_graph = nx.Graph()
    #     relabeled_graph.add_edges_from(new_edge_set)

    #     new_labeling, new_automorphisms = graph_canon(relabeled_graph)
    #     print(f"Resulting in labeling {new_labeling}")

    #     canonical_edge_set = permute_edges(new_labeling, relabeled_graph.edges)
    #     new_canonical_adj = create_adjacency(9, canonical_edge_set)

    #     if not np.array_equal(canon_adj, new_canonical_adj):
    #         print("Graphs inequal")
    #         break
    #     iterator += 1


## NOTE, ADJACENCY MATRIX DOES NOT SEEM TO BE UPDATED UPON RELABELLING