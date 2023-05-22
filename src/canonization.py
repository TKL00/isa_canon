import networkx as nx
from TreeNode import TreeNode
import copy
import numpy as np
from itertools import permutations
import time


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

def graph_canon(G, Q, traces=False):
    """
        Produces the canonical labeling of the Graph G using the Graph Canonicalization 
        algorithm suggested by Brendan McKay.

        `Paramters`:
            G (Graph): NetworkX Graph
            Q (function ( list( list(int) ) ) -> list): The target cell selector that, given a non-discrete partition, 
                                                        returns a reference to a non-trivial part in the partition.
            traces (bool): Declares whether the node invariant on cell spilt positions should be used.

        `Returns`:
            Canonical_labeling (dict: node -> node): Relabeling mapping representing G in its canonical representation
            Automorphisms (list( dict: node -> node )): The generating set for the automorphism group of G.
    """

    trace_had_impact = False

    G_NODE_AMT = len(G.nodes)

    ## list of automorphism (dicts)
    automorphisms = []
    ## the global minimum consist of 1) a partition and 2) the resulting adjacency matrix when permuting G.
    global_invariants = {
        "least_partition": [],
        "least_adjacency": [],
        "max_trace": [[]]
    }

    def less_than_max(current_trace):

            max_trace = global_invariants["max_trace"]
            
            ## Determine what level the current trace is at in the search tree
            tree_depth = min(len(max_trace), len(current_trace))
            
            ## Gets the specific trace introduced at that level in the max trace,
            ## and in the current level in the current trace
            node_max_trace = max_trace[tree_depth - 1]
            node_current_trace = current_trace[-1]
            
            ## Determine how many elements should be compared
            prefix_length = min(len(node_max_trace), len(node_current_trace))

            ## Slicing of the nodes

            compare_num_max_trace = node_max_trace[0:prefix_length]
            compare_num_current_trace = node_current_trace[0:prefix_length]

            return compare_num_current_trace < compare_num_max_trace 

    def greater_than_max(current_trace):

            max_trace = global_invariants["max_trace"]
            
            ## Determine what level the current trace is at in the search tree
            tree_depth = min(len(max_trace), len(current_trace))
            
            ## Gets the specific trace introduced at that level in the max trace,
            ## and in the current level in the current trace
            node_max_trace = max_trace[tree_depth - 1]
            node_current_trace = current_trace[-1]
            
            ## Determine how many elements should be compared
            prefix_length = min(len(node_max_trace), len(node_current_trace))

            ## Slicing of the nodes

            compare_num_max_trace = node_max_trace[0:prefix_length]
            compare_num_current_trace = node_current_trace[0:prefix_length]

            return compare_num_current_trace > compare_num_max_trace 
    
    def equitable_refinement(p, v):
        """
            Produces an equitable ordered partition based on the input partition 'p' of the graph 'G'
            using the algorithm suggested by Brendan McKay.

            `Parameters`
                p (list( list(int) )): List of partitioned nodes (in lists).

                v (TreeNode): The current node in the search tree
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

            ## When using traces, we add the trace of this shatter to the node's current trace
            if(traces):
                split_position = 0
                for i in range(to_shatter):
                    split_position += len(p[i])
                for i in range(len(shatter_partitions) -1):
                    split_position += len(shatter_partitions[i])
                    v.get_trace()[-1].append(split_position)
                                                       
            ## Simply replacing the partition V_i with the partitions in shatter_partition
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

        ## NOTE: Abort refinement if current partition is lexicographically smaller than the global minimum
        while len(B) != 0:
            if traces and less_than_max(v.get_trace()):
                # print("Skipped branch due to trace")
                trace_had_impact = True
                return False

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
                and returns the position of the split
            
            """
            new_p = []
            split_position = 0
            has_individualized = False
            for part in p:
                if not has_individualized:
                    if v in part:
                        split_position += 1
                        ## Create and insert {v}
                        new_p.append([v])
                        new_part = []
                        ## Create V_i \ {v}
                        for nodes in part:
                            if nodes != v:
                                new_part.append(nodes)
                        ## Insert V_i \ {v}
                        new_p.append(new_part)
                        has_individualized = True
                    else: 
                        split_position += len(part)
                        new_p.append(part)
                ## All parts not containing v stays the same in the same relative position in the partition
                else:
                    new_p.append(part)
            
            return new_p, split_position

        def update_automorphisms(new_leaf_partition, new_leaf_adj):            
            ## If partitions are isomorphic, compute automorphism between them
            if np.array_equal(new_leaf_adj, global_invariants["least_adjacency"]):
                best_partition = global_invariants["least_partition"]
                ## pi_one goes from the canonical labelling to the input graph, the inverse should go from the input graph to the canonical labelling
                ## partition[i][0] = j --> whatever node at position "i" in the partition 
                ## corresponds to a node in the original graph and is mapped to node (color) 'j' in canonical graph
                pi_one_inverse = {new_leaf_partition[i][0]: i for i in range(G_NODE_AMT)}
                pi_two = {i: best_partition[i][0] for i in range(G_NODE_AMT)}
                new_automorphism_math = {i: pi_two[pi_one_inverse[i]] for i in range(G_NODE_AMT)}

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

        def calculate_orbit(child, fixed_automorphisms):
            """
                Calculates the orbits for each node child in the children list under closure of the automorphism,
                thereby partitioning the children list, where each orbit is a part.

                `Paramters`:

                    children ((int)): The child to calculate the orbit of

                    automorphism (list( dict(int -> int) )): List of automorphisms represented as dictionaries (mapping one node to another)

                `Returns`:
                    all_orbits (list( list(int))): Partition of children_list corresponding to the orbits of the children.
            """

            
            L_i = set([child])
            found_members = set([child])

            ## as long as new members of orbits are encountered, we continue
            while L_i:
                L_i_s = set([])
                for vals in L_i:
                    for automorphism in fixed_automorphisms:
                        L_i_s.add(automorphism[vals])
                L_j = L_i_s - found_members
                for new_discovery in L_j:
                    found_members.add(new_discovery)
                
                L_i = L_j
            
            return sorted(list(found_members))

        def generate_subtree(parent_node, partition, current_seq, to_indiv):
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

            ## If parent is root node, the current_seq is empty and will only contain this individualized node
            if not current_seq:
                new_sequence = [to_indiv]
            else:
                new_sequence = copy.deepcopy(current_seq)
                new_sequence.append(to_indiv)

            indiv_partition, split_position = individualize(partition, to_indiv)
            
            ## create dummy node with parent's partition and append individualized node to the split
            this_node = TreeNode(partition, parent_node, new_sequence)
            if traces: 
                this_node.set_trace(copy.deepcopy(parent_node.get_trace()))
                this_node.get_trace().append([split_position])
            
            refinement = equitable_refinement(indiv_partition, this_node)
            ## Only applicable for when traces is true
            if traces and not refinement:
                return
            else: 
                ## update node with refined partition
                this_node.set_partition(refinement)
            
            ## Target cell selector
            children_list = copy.deepcopy(Q(refinement))

            ## If no choices available, this is a leaf node
            if len(children_list) == 0:
                
                ## create mapping from current discrete partition. Note that partition[i] says that the node (= partition[i]) of the input graph 
                ## goes to color "i" i.e. node "i" in the canonical labelling.
                relabeling_map = {this_node.get_partition()[i][0]: i for i in range(G_NODE_AMT)}
                relabeled_edges = permute_edges(relabeling_map, G.edges)
                leaf_adj = create_adjacency(G_NODE_AMT, relabeled_edges)

                ## First encountered leaf is the global minimum
                if not global_invariants["least_partition"]:
                    global_invariants["least_partition"] = this_node.get_partition()
                    global_invariants["least_adjacency"] = leaf_adj
                    if traces: ## initial trace
                        #print("first trace")
                        global_invariants["max_trace"] = this_node.get_trace()                    

                ## Overwrite max trace if leaf's trace is larger
                if traces and greater_than_max(this_node.get_trace()):
                    #print(f"overwrite trace")
                    global_invariants["max_trace"] = this_node.get_trace()

                
                new_leaf_partition = this_node.get_partition()
                ## Check if automorph with current best leaf. 
                update_automorphisms(new_leaf_partition, leaf_adj)
                
                ## Check if current relabeling is better than previous relabeling. Resulting array contains true if one elemen
                if array_less_than(leaf_adj, global_invariants["least_adjacency"]):
                    global_invariants["least_partition"] = this_node.get_partition()
                    global_invariants["least_adjacency"] = leaf_adj
                ## Check for automorphism between other leaf partitions
            else:
                this_node.set_children(children_list)
                travers_seq = this_node.get_travers_seq()

                ## for each child visited, calculate new orbits and decide whether remaining children
                ## should be visisted
                for child in this_node.get_children():
                    ## orbit containing this child
                    A_prime = list(filter(lambda auto: all_fixed(travers_seq, auto), automorphisms))
                    orbit = calculate_orbit(child, A_prime)
                    ## Only the first element of the orbit should be discovered. Since everything is lexicographically sorted, a child
                    ## that appears late in an orbit should not be traversed.
                    if orbit[0] == child: generate_subtree(this_node, refinement, new_sequence, child)

                                             ## GENERATE TREE ROOT  

        ## Generate the root node's canonical
        for child in root.get_children():
            orbit = calculate_orbit(child, automorphisms)
            if orbit[0] == child: generate_subtree(root, root.get_partition(), [], child)

        return
    
    ##                                  CANONICAL

    root_partition = [sorted(list(G.nodes))]
    root_node = TreeNode([], None, [])
    if traces: root_node.set_trace([[]])
    init_refinement = equitable_refinement(root_partition, root_node)
    root_node.set_partition(init_refinement)
    # if traces: global_invariants["max_trace"] = root_node.get_trace()

    ## Target cell selector
    children_list = copy.deepcopy(Q(init_refinement))
    
    ## If initial equitable partition is discrete, the canonical representation has been achieved.
    if not children_list:
        return {init_refinement[i][0]:i for i in range(G_NODE_AMT)}, [], trace_had_impact

    root_node.set_children(children_list)

    ## generate_tree has side effects on the global invariants dict
    generate_tree(root_node)
    canonical_partition = global_invariants["least_partition"]
    canonical_labeling = {canonical_partition[i][0]:i for i in range(G_NODE_AMT)}

    return canonical_labeling, automorphisms, trace_had_impact

def test_canon(G, Q, traces=False):
    """
        Runs the graph_canon function on G and checks whether it correctly computes the 
        canonical form of G.
    """

    n_nodes = len(G.nodes)
    perm = np.random.permutation(n_nodes)
    perm_mapping = {i: perm[i] for i in range(n_nodes)}
    perm_edges = permute_edges(perm_mapping, G.edges)
    perm_G = nx.Graph()
    perm_G.add_nodes_from(G.nodes)
    perm_G.add_edges_from(perm_edges)

    labeling_1, automorphisms_1 = graph_canon(G, Q, traces)
    print(f"len of automorphisms {len(automorphisms_1)}")
    labeling_2, automorphisms_2 = graph_canon(perm_G, Q, traces)
    print(f"len of automorphisms {len(automorphisms_2)}")

    G_canon_edge = permute_edges(labeling_1, G.edges)
    perm_G_canon_edge = permute_edges(labeling_2, perm_G.edges)
    
    G_canon_adj = create_adjacency(n_nodes, G_canon_edge)
    perm_G_canon_adj = create_adjacency(n_nodes, perm_G_canon_edge)

    return np.array_equal(G_canon_adj, perm_G_canon_adj)
    

def first_non_trivial(partition):
    for part in partition:
        if len(part) != 1:
            return part
    return []

def rightmost_first_non_trivial(partition):
    for part in reversed(partition):
        if len(part) != 1:
            return part
    return []


# if __name__ == "__main__":

    # graph = nx.Graph()
    # graph.add_nodes_from([i for i in range(9)])
    # graph.add_edges_from([(0, 1), (0, 3), (1, 2), (1, 4), (2, 5), (3, 4), (3, 6), (4, 5), (4, 7), (5, 8), (6, 7), (7, 8)])

    # res, autos = graph_canon(graph, first_non_trivial, True)
    # print(len(autos))

    # for i in range(100):
    #     g = nx.dense_gnm_random_graph(20, 170)

    #     res = test_canon(g, first_non_trivial, True)
    #     if not res:
    #         break
    #     else:
    #         print(res)


    # new_graph = nx.Graph()
    # new_graph.add_nodes_from([i for i in range(9)])
    # new_graph.add_edges_from([(0, 1), (0, 3), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5), (3, 6), (3, 7), (4, 7), (5, 7), (5, 8), (6, 7), (7, 8)])

    # # frucht_graph = nx.frucht_graph()

    # res, autos = graph_canon(frucht_graph, first_non_trivial, True)
    # print(len(autos))
    # print(autos)