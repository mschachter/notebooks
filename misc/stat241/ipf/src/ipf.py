
import copy

import numpy as np
import networkx as nx

class FunctionTable():
    
    def __init__(self, config_table):
        self.config_table = config_table
    
    def eval(self, vals):
        vkey = tuple(vals)
        if vkey not in self.config_table:
            print 'ERROR: configuration not in config_table: %s' % str(vals)
        return self.config_table[vkey]
    

class UndirectedBinaryGraph():
    """ This class represents an undirected probablistic graph whose nodes are binary-valued """

    def __init__(self):
        #self.g is the networkx graph that embodies the nodes and connections        
        self.g = nx.Graph()
        
        #self.values is a map where the key is a node id and the value is a
        #list of possible values that n can take on
        self.values = {}
        
        #self.potentials is a map, where the key is a tuple representing the
        #node #s that are arguments of the potential, and the value is a
        #FunctionTable object
        self.potentials = {}
        
        #the normalizing constant
        self.Z = None
        
        #relevant cliques with associated potential functions
        self.cliques = []
        
        #associates a node ID with an index in self.nodes
        self.node_index = {}
    
    def add_node(self, id):
        """ Add node, given it's id, potential function """
        self.g.add_node(id)
        self.node_index[id] = len(self.node_index)       
        
    def add_clique(self, clique, potential):
        """ Add a list of nodes as a clique """
        self.cliques.append(clique)
        self.potentials[tuple(clique)] = potential
        for n1 in clique:
            for n2 in clique:
                e = [n1, n2]
                e.sort()
                self.g.add_edge(e[0], e[1])
            
    def probability(self, x):
        """ Compute the probability of a given set of values.
            
            Argument x is a map, key is the node id, value is the value for that node.
        """
        prod = 1
        for c in self.cliques:
            cval = [x[self.node_index[n]] for n in c]
            pfunc = self.potentials[tuple(c)] 
            prod *= pfunc.eval(cval)
        Z = self.normalizing_constant()
        return prod / Z
    
    def summary(self):
        """ Print a summary of all clique potentials """
        for c in self.cliques:
            print 'Clique: %s' % str(c)
            bcombs = get_all_binary_combos(len(c))
            for bc in bcombs:
                print '\t%s: %0.3f' % (str(bc), self.potentials[tuple(c)].eval(bc))
            
        
    def normalizing_constant(self):
        """ Compute the partition function for this graph """
        Z = 0
        #assumes all nodes are binary valued...
        bcombs = get_all_binary_combos(len(self.g.nodes()))
        for bc in bcombs:            
            cprod = 1
            for clique in self.cliques:
                bval = [bc[ self.node_index[n] ] for n in clique]                    
                pfunc = self.potentials[tuple(clique)]
                cprod *= pfunc.eval(bval)

            Z += cprod
            
        self.Z = Z
        
        return self.Z


    def __str__(self):
        rstr = ''
        rstr += 'Undirected Probability Tree:\n'
        rstr += '\tNodes: %s\n' % [x for x in self.g.nodes()]
        rstr += '\tEdges: %s\n' % [(x[0], x[1]) for x in self.g.edges()]
        rstr += '\tPotentials: %s\n' % self.potentials.keys()
        rstr += '\tValues: %s\n' % self.values
        
        return rstr


class IPF():
    
    def __init__(self, g, data):
        self.ubg = g
        self.data = data
        self.empirical_marginals = self.compute_empirical_marginals(data)
        self.N = len(data)
        
    def compute_empirical_marginals(self, data):
        """ Compute the count of each clique configuration in the data """
        marginals = {}
        for row in data:
            for clique in self.ubg.cliques:
                vals = []
                for n in clique:
                    vals.append(row[ self.ubg.node_index[n] ])
                ckey = tuple(clique)
                if ckey not in marginals:
                    marginals[ckey] = {}
                vkey = tuple(vals)
                if vkey not in marginals[ckey]:
                    marginals[ckey][vkey] = 0
                marginals[ckey][vkey] += 1
                
        #fill in values for any configurations not in the data
        for ckey,cmap in marginals.iteritems():
            bcombs = get_all_binary_combos(len(ckey))
            for bc in bcombs:
                vkey = tuple(bc)
                if vkey not in cmap:
                    cmap[vkey] = 0
            
        return marginals
        
    
    def compute_clique_marginal(self, clique):
        """ Marginalize a clique's potential with respect to all nodes not in the clique """
        Z = self.ubg.normalizing_constant()
        clique_bcombs = get_all_binary_combos(len(clique))        
        ctable = {}
        for clique_vals in clique_bcombs:
            num_other_nodes = len(self.ubg.g.nodes()) - len(clique)
            if num_other_nodes == 0:
                x = clique_vals
                psum = self.ubg.probability(x)                
            else:
                other_node_bcombs = get_all_binary_combos(num_other_nodes)
                psum = 0
                for onbc in other_node_bcombs:                
                    x = [0]*len(self.ubg.g.nodes())
                    last_index = 0
                    for k,n in enumerate(self.ubg.g.nodes()):
                        if n in clique:
                            cindx = clique.index(n)
                            x[k] = clique_vals[cindx]
                        else:
                            x[k] = onbc[last_index]
                            last_index += 1 
                    psum += self.ubg.probability(x)
            ctable[tuple(clique_vals)] = psum
        
        cmarg = FunctionTable(ctable) 
        return cmarg         
    
    
    def iterate(self):
        """ Run a single iteration of IPF across each clique """
        
        Zold = self.ubg.normalizing_constant()
        for clique in self.ubg.cliques:
            #get current potential function
            ckey = tuple(clique)
            old_pfunc = self.ubg.potentials[ckey]
            
            #compute the clique marginal (a table) for this clique            
            pc = self.compute_clique_marginal(clique)
            
            #update the value of the potential function for each
            #possible configuration of the clique            
            ctable = {}
            bcombs = get_all_binary_combos(len(clique))
            for bc in bcombs:
                vkey = tuple(bc)
                old_pfunc_val = old_pfunc.eval(bc)                 
                pemp = float(self.empirical_marginals[ckey][vkey]) / float(self.N)
                pc_val = pc.eval(vkey)
                print '%s, P_hat[%s]=%0.3f, P[%s]=%0.3f' % (str(clique), str(bc), pemp, str(bc), pc_val)
                update_val = 0.0
                if pemp > 0.0 and pc_val > 0.0:
                    update_val = old_pfunc_val*(pemp / pc_val)
                ctable[vkey] = update_val
            new_pfunc = FunctionTable(ctable)
            self.ubg.potentials[ckey] = new_pfunc
        Znew = self.ubg.normalizing_constant()
        print 'Zold=%0.4f, Znew=%0.4f' % (Zold, Znew)            
            
    def log_likelihood(self):
        L = 0
        for x in self.data:
            p = self.ubg.probability(x)
            if p > 0.0:
                L += np.log(p)
        return L
    
    def check(self):
        """ Check the empirical marginals against the current potentials """
        for clique in self.ubg.cliques:
            print 'Clique: %s' % str(clique)
            ckey = tuple(clique)
            em_table = self.empirical_marginals[ckey]
            for x,mempirical in em_table.iteritems():
                pempirical = float(mempirical) / float(self.N)
                pmodel = self.ubg.potentials[ckey].eval(x)                
                print '\t%s: empirical=%0.3f, model=%0.3f' % (str(x), pempirical, pmodel)
        
        
def get_all_binary_combos(nvals):
    """ Gets all possible combinations for a binary vector of length nvals """
    
    clist = []
    indx2change = 0
    cvals = [0]*nvals
    _get_all_binary_combos(clist, cvals, indx2change)
    
    return clist
    
def _get_all_binary_combos(clist, cvals, indx2change):
    """ Recursive helper function for get_all_binary_combos """
    
    for b in (0, 1):
        cv = copy.deepcopy(cvals)
        cv[indx2change] = b
        if cv not in clist: #cheap fix...
            clist.append(cv)
        #print 'cvals=%s, cv=%s, indx2change=%d, b=%d' % (cvals, cv, indx2change, b)
        if indx2change < len(cvals)-1:
            _get_all_binary_combos(clist, cv, indx2change+1)
