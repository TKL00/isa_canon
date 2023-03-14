import networkx as nx
from TreeNode import TreeNode
import copy
import numpy as np

def graph_canon(G):
    """
        Produces the canonical labelling of the Graph G using the Graph Canonicalization 
        algorithm suggested by Brendan McKay.
    """
    
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

        adj_matrix = nx.to_numpy_array(G)

        tau = p
        
        B = shatter_set(p, adj_matrix)

        while len(B) != 0:
            ## Choose individualizing the lexicographically smallest
            (i, j) = sorted(B)[0]
            tau = shatter_partition(tau, adj_matrix, i, j)
            ## Update B to see if partitions remains to be shattered
            B = shatter_set(tau, adj_matrix)

        return tau
   
    def generate_tree(root):

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

        def update_automorphisms(leaves, automorphisms):
            
            print(f"\t\t\t Updating automorphism list")
            new_labels = {}

            ## permute G to the newly added leaf
            new_leaf = leaves[-1]
            for i in range(len(new_leaf)):
                part = new_leaf[i]
                new_labels[part[0]] = i

            leaf_graph = nx.relabel_nodes(G, new_labels)
            leaf_adj = nx.to_numpy_array(leaf_graph)
            print(f"Newly added partition\t: {new_leaf}")

            ## permute G to all other leaves except the last (new) leaf
            for leaf in range(len(leaves) - 1):
                comp_leaf = leaves[leaf]
                
                for i in range(len(comp_leaf)):
                    part = comp_leaf[i]
                    ## In the partition, the node in part 'i' is mapped to node i.
                    new_labels[part[0]] = i
                comp_graph = nx.relabel_nodes(G, new_labels)
                comp_adj = nx.to_numpy_array(comp_graph)

                if np.array_equal(leaf_adj, comp_adj):
                    print(f"Adjacency matrix is equal to\t: {comp_leaf}")
                    ## NOTE: BEGIN FROM HERE TO COMPUTE AUTOMORPHISM. USE new_labels[i] = parts[0]
                    ## TO COMPUTE THE INVERSE OF A LEAF.
                    ## THEN DO "new_leaf^-1 * leaf" TO GET AUTOMORPHISM


        def generate_subtree(parent_node, partition, current_seq, to_indiv, leaves, automorphisms):
            """
                Subroutine to generate the rest of the tree spanning from this node

                `Parameters`:

                    parent_node (TreeNode): A TreeNode representing the parent node in the tree of the node being generated.

                    partition (list (list(int))): The current partition of the nodes in the graph.
                    
                    current_seq (list(int)): The current sequence of individualization that led to (including) this TreeNode's partition.

                    to_indiv (int): The element in the partition being individualized.

            """
            indiv_partition = individualize(partition, to_indiv)
            refinement = equitable_refinement(indiv_partition)
            
            ## If parent is root node, the current_seq is empty and will only contain this individualized node
            if not current_seq:
                new_sequence = [to_indiv]
            else:
                new_sequence = copy.deepcopy(current_seq).append(to_indiv)

            this_node = TreeNode(refinement, parent_node, new_sequence)

            ## Check for automorphism/node invariants -> if to be pruned, just return.

            ## Find first non-trivial part of the refined partition
            children_list = []
            for i in range(len(refinement)):
                if len(refinement[i]) != 1:
                    children_list = copy.deepcopy(refinement[i])
                    break
            ## If no choices available, this is a leaf node
            if len(children_list) == 0:
                leaves.append(this_node.get_partition())
                update_automorphisms(leaves, automorphisms)
                ## Check for automorphism between other leaf partitions
            else:
                this_node.set_children(children_list)
                for child in children_list:
                    ## NOTE: Check for pruning
                    generate_subtree(this_node, refinement, new_sequence, child, leaves, automorphisms)

                                        ## ROOT
        ## list of automorphism 
        automorphisms = []
        ## list of leaf partitions
        leaves = []

        ## Generate the root node's canonical
        for child in root.get_children():
            ## NOTE: Check for pruning here
            generate_subtree(root, root.get_partition(), [], child, leaves, automorphisms)

    root_partition = [list(G.nodes)]
    init_refinement = equitable_refinement(root_partition)   
    root_node = TreeNode(init_refinement, None, [])

    ## Find first non-trivial part of the refined partition
    children_list = []
    for i in range(len(init_refinement)):
        if len(init_refinement[i]) != 1:
            children_list = copy.deepcopy(init_refinement[i])
            break

    root_node.set_children(children_list)

    generate_tree(root_node)

    return None


if __name__ == "__main__":

    graph = nx.Graph()
    graph.add_nodes_from([i for i in range(9)])
    graph.add_edges_from([(0, 1), (0, 3), (1, 2), (1, 4), (2, 5), (3, 4), (3, 6), (4, 5), (4, 7), (5, 8), (6, 7), (7, 8)])

    print("G's adjacency list")
    print(nx.to_numpy_array(graph))

    graph_canon(graph)