"""
Created on Oct 6, 2011

@author: Mike Schachter (mike.schachter@gmail.com)
"""

import networkx as nx

from sum_product import UPGraph, SumProduct


def create_2node_graph():
    """ Simple 2 node undirected graph """
    
    g = UPGraph()
    g.add_node(1, single_potential_1, [0, 1])
    g.add_node(2, single_potential_2, [0, 1])
    g.connect(1, 2, pairwise_potential)
    
    return g

def create_3node_graph():
    """ Simple 3 node undirected graph """
    
    g = UPGraph()
    g.add_node(1, single_potential_1, [0, 1])
    g.add_node(2, single_potential_2, [0, 1])
    g.add_node(3, single_potential_1, [0, 1])
    g.connect(1, 2, pairwise_potential)
    g.connect(1, 3, pairwise_potential)
    
    return g


def create_problem3b_graph():
    """" Creates the graph structure, and associates potentials with graph structure """
    g = UPGraph()
    
    #create 6 nodes in graph
    for k in range(6):
        id = k+1
        if id in (1, 3, 5):
            pfunc = single_potential_1
        else:
            pfunc = single_potential_2
        g.add_node(id, pfunc, [0, 1])
    
    #add edges between relevant nodes for figure 1
    g.connect(1, 2, pairwise_potential)
    g.connect(1, 3, pairwise_potential)
    g.connect(2, 4, pairwise_potential)
    g.connect(2, 5, pairwise_potential)
    g.connect(3, 6, pairwise_potential)
    
    return g


def pairwise_potential(node1_val, node2_val):
    """ Implementation of the pairwise potential function """
    mat = [[1, 0.5], [0.5, 1]]
    return mat[node1_val][node2_val]

        
def single_potential_1(node_val):
    """ Implementation of the single node potentials for nodes 1, 3, 5 """
    mat = [0.7, 0.3]
    return mat[node_val]


def single_potential_2(node_val):
    """ Implementation of the single node potentials for nodes 2, 4, 6 """
    mat = [0.1, 0.9]
    return mat[node_val]



def simple_example():
    
    simple_upg = create_2node_graph()
    print simple_upg  
    sp = SumProduct(simple_upg)    
    
    while not sp.blocked:
        sp.iterate()
    
    sp.print_marginals()
    
def medium_example():
    
    upg = create_3node_graph()
    print upg  
    sp = SumProduct(upg)    
    
    while not sp.blocked:
        sp.iterate()
    
    sp.print_marginals()
    
    
    
def problem_3b():
    
    upg = create_problem3b_graph()
    print upg
        
    sp = SumProduct(upg)    
    
    while not sp.blocked:
        sp.iterate()
    
    sp.print_marginals()
    
    
if __name__ == '__main__':
    problem_3b()
