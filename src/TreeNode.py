import copy

class TreeNode:

    """
        Data representation of a node object in the canonicalization search three.
    """

    def __init__(self, partition, parent, travers_seq) -> None:
        """
            `Parameters`:

                partition (list(list(Nodes))): The partition of nodes in a graph G traversed in this structure.
            
                parent (TreeNode): A reference to the parent node of this node.

                travers_seq (list(TreeNode)): A list of nodes such that travers_seq[i] was the i'th individualized node

                children (list(TreeNode)): A list of possible branches to make depending on choices of individualization
        """
        self._partition = partition
        self._parent = parent
        self._travers_seq = travers_seq
        self._children = []
    
    def get_partition(self):
        return self._partition

    def get_parent(self):
        return self._parent
    
    def get_travers_seq(self):
        return self._travers_seq
    
    def get_children(self):
        return self._children
    
    def set_children(self, c):
        self._children = c
    

    
    