
import numpy as np

from ipf import UndirectedBinaryGraph, FunctionTable, IPF, get_all_binary_combos


def create_simple_graph():
    
    ubg = UndirectedBinaryGraph()
    
    #construct graph
    ubg.add_node(1)
    
    #construct and set potential functions
    config = {}
    config[(0,)] = 0.8
    config[(1,)] = 0.2
    ftable = FunctionTable(config)    
    ubg.add_clique([1], ftable)
    
    #create data to go with it
    N = 100
    data = []    
    rnums = np.random.rand(N)
    for r in rnums:
        if r <= 0.30:
            data.append( [0] )
        else:
            data.append( [1] )
        
    return (ubg, data)

def read_data_file(fpath):
    data = []
    f = open(fpath, 'r')
    for ln in f.readlines():
        s = ln.split(' ')
        r = [float(x.strip()) for x in s if len(x.strip()) > 0]
        data.append(r)
    f.close()
    return data

def create_prob2a_tree(single_node_cliques=True):
    
    ubg = UndirectedBinaryGraph()
    
    #construct graph
    ubg.add_node(1)
    ubg.add_node(2)
    ubg.add_node(3)
    ubg.add_node(4)    
    
    if single_node_cliques:
        #construct single node potentials
        for n in ubg.g.nodes():
            config = {}
            rnum = np.random.rand()
            config[(0,)] = rnum
            config[(1,)] = 1 - rnum
            ftable = FunctionTable(config)    
            ubg.add_clique([n], ftable)
    
    #construct and set potential functions
    cliques = [ (1, 2), (2, 3), (3, 4)]
    for c in cliques:
        bcombs = get_all_binary_combos(len(c))
        rnums = np.random.rand(len(bcombs))
        probs = rnums / rnums.sum()
        config = {}            
        for k,bc in enumerate(bcombs):            
            config[tuple(bc)] = probs[k]
        ftable = FunctionTable(config)
        ubg.add_clique(c, ftable)
    
    data = read_data_file('../IPF.dat')
        
    return (ubg, data)


def create_prob2a_full():
    
    ubg = UndirectedBinaryGraph()
    
    #construct graph
    ubg.add_node(1)
    ubg.add_node(2)
    ubg.add_node(3)
    ubg.add_node(4)    
    
    #construct single node potentials    
    for n in ubg.g.nodes():
        config = {}
        rnum = np.random.rand()
        config[(0,)] = rnum
        config[(1,)] = 1 - rnum
        ftable = FunctionTable(config)    
        ubg.add_clique([n], ftable)    
    
    #construct and set potential functions
    #cliques = [ (1, 2, 3, 4)]
    cliques = [ (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    for c in cliques:
        bcombs = get_all_binary_combos(len(c))
        rnums = np.random.rand(len(bcombs))
        probs = rnums / rnums.sum()
        config = {}            
        for k,bc in enumerate(bcombs):            
            config[tuple(bc)] = probs[k]
        ftable = FunctionTable(config)
        ubg.add_clique(c, ftable)
    
    data = read_data_file('../IPF.dat')
        
    return (ubg, data)


def run_simple_graph():
    
    (ubg, data) = create_simple_graph()    
    ipf = IPF(ubg, data)
    ipf.iterate()
    ipf.iterate()
    
    return ipf
    

def run_prob2a_tree(single_node_cliques=True):
    (ubg, data) = create_prob2a_tree(single_node_cliques)    
    ipf = IPF(ubg, data)
    
    for k in range(5):
        ipf.iterate()
    
    return ipf

def run_prob2a_full():
    (ubg, data) = create_prob2a_full()    
    ipf = IPF(ubg, data)
    
    return ipf
    
        