""" Ising Model code for Stat 241 """

import copy as pycopy
import time

import numpy as np
import networkx as nx


class IsingModel():
    
    def __init__(self, mean_field=False):
        self.g = nx.Graph()
        self.potentials = {}
        self.state = {}
        self.mean_field = mean_field
    
    def add_node(self, id, potential):
        self.g.add_node(id)
        self.potentials[id] = potential
        
        if not self.mean_field:
            #random initialization with binary state
            r = np.random.rand()
            if r < 0.50:
                self.state[id] = -1
            else:
                self.state[id] = 1
        else:
            #random initialization with value in [-1, 1]
            r = np.random.rand()*2 - 1
            self.state[id] = r
    
    
    def connect(self, id1, id2, edge_potential):
        e = [id1, id2]
        e.sort()
        ekey = (e[0], e[1])
        if ekey not in self.potentials:
            self.g.add_edge(e[0], e[1])        
            self.potentials[ekey] = edge_potential
            
    def compute_node_potential(self, n, val=None):
        if n in self.potentials:
            if val is None:
                val = self.state[n]
            return self.potentials[n](val)
        return 0.0
    
    
    def compute_edge_potential(self, n1, n2, val1=None, val2=None):
        e = [n1, n2]
        e.sort()
        ekey = (e[0], e[1])
        if ekey in self.potentials:
            if val1 is None:
                val1 = self.state[n1]
            if val2 is None:
                val2 = self.state[n2]
            return self.potentials[ekey]([val1, val2])
        return 0.0
    
    
    def sample(self, nsamps, nburn=1000):
        
        unsampled = self.state.keys()
        #burn in time...
        for k in range(nburn):
            if len(unsampled) == 0:
                unsampled = self.state.keys()
            rindx = np.random.randint(len(unsampled))
            nid = unsampled[rindx]
            self.update_sample_gibbs(nid)
        
        samples = []        
        for k in range(nsamps):
            if len(unsampled) == 0:
                unsampled = self.state.keys()
            rindx = np.random.randint(len(unsampled))
            nid = unsampled[rindx]
            self.update_sample_gibbs(nid)
            samples.append(pycopy.copy(self.state))
        
        return samples        
        
    
    def update_sample_gibbs(self, nid):
        #compute energy given current state and x_nid = 1
        e = self.compute_node_energy(nid, 1)
        
        if not self.mean_field:        
            #compute P[x_nid = 1 | x]
            p1 = 1.0 / (1.0 + np.exp(-e))
            
            #sample from this binary distribution
            rnum = np.random.rand()
            if rnum <= p1:
                s = 1
            else:
                s = -1
            #print 'nid=%d, e=%0.3f, p1=%0.3f, s=%d' % (nid, e, p1, s)    
            
        else:
            s = e
        
        self.state[nid] = s
    
    def update_sample_mean_field(self, nid):
        #compute energy given current state and x_nid = 1
        e1 = self.compute_energy(nid)
        
        #compute energy given current state and x_nid = 0
        e0 =  self.compute_energy(nid, 0)        
        
        #compute P[x_nid = 1 | x]
        p1 = e1 / (e0 + e1)
        
        #sample from this binary distribution
        rnum = np.random.rand()
        if rnum < p1:
            s = 1
        else:
            s = -1
            
        #print 'nid=%d, e0=%0.3f, e1=%0.3f, p1=%0.3f, s=%d' % (nid, e0, e1, p1, s)
        self.state[nid] = s
        
    
    def compute_node_energy(self, nid, nval, clamped_states={}):
        """ Compute energy of node nid given it's neighbors. """
        
        e = self.compute_node_potential(nid, nval)
        for nbid in self.g[nid].keys():
            nbval = None
            if nbid in clamped_states:
                nbval = clamped_states[nbid]            
            e += self.compute_edge_potential(nid, nbid, nval, nbval)
            
        return e
        


def create_prob54_torous(num_rows, num_cols, node_potential=0.2, edge_potential=0.2, mean_field=False):
    
    imodel = IsingModel(mean_field=mean_field)
    
    #create nodes
    for j in range(num_cols):
        for i in range(num_rows):
            s = grid2id(i, j, num_rows, num_cols)
            if (s+1) % 2 == 0:
                pfunc = lambda (x): node_potential + 1 
            else:
                pfunc = lambda (x): node_potential - 1
            imodel.add_node(s, pfunc)
    
    #connect nodes    
    for j in range(num_cols):
        for i in range(num_rows):
            s = grid2id(i, j, num_rows, num_cols)
            nbrs = get_torous_neighbors(i, j, num_rows, num_cols)
            efunc = lambda(x): edge_potential
            for nid in nbrs:                
                imodel.connect(s, nid, efunc)
    
    return imodel


def grid2id(i, j, num_rows, num_cols):
    return num_cols*i + j
    

def get_torous_neighbors(i, j, num_rows, num_cols):
    nbrs = []
    
    #get neighbor above
    ni = (i-1) % num_rows 
    nj = j
    ns = grid2id(ni, nj, num_rows, num_cols)
    nbrs.append(ns)
    
    #get neighbor below
    ni = (i+1) % num_rows 
    nj = j
    ns = grid2id(ni, nj, num_rows, num_cols)
    nbrs.append(ns)
    
    #get neighbor to left
    ni = i 
    nj = (j-1) % num_cols
    ns = grid2id(ni, nj, num_rows, num_cols)
    nbrs.append(ns)
    
    #get neighbor to right
    ni = i 
    nj = (j+1) % num_cols
    ns = grid2id(ni, nj, num_rows, num_cols)
    nbrs.append(ns)
    
    return nbrs    


def show_samples(nsamps=500, mean_field=False):
    
    import matplotlib.cm as cm
    from pylab import *
    
    num_rows = 7
    num_cols = 7
    imodel = create_prob54_torous(num_rows, num_cols, node_potential=0.2, edge_potential=0.2, mean_field=mean_field)
    
    samps = imodel.sample(nsamps, nburn=0)
    
    num_nodes = num_rows*num_cols
    all_samps = np.zeros([nsamps, num_nodes])    
    for i,samp in enumerate(samps):
        for j,val in samp.iteritems():            
            all_samps[i, j] = val
    
    ion()
    for k in range(nsamps):
        mat = np.reshape(all_samps[k, :], [num_rows, num_cols])
        imshow(mat, cmap=cm.gray, interpolation='nearest')
        title('Sample %d' % k)
        time.sleep(0.500)
        draw()
        
def compute_gibbs_means(ntries=10):

    num_rows = 7
    num_cols = 7
    imodel = create_prob54_torous(num_rows, num_cols, node_potential=0.2, edge_potential=0.2, mean_field=False)

    gibbs_means = np.zeros([ntries, num_rows, num_cols])
    
    for k in range(ntries):
        nsamps = 1000
        samps = imodel.sample(nsamps, nburn=1000)
        num_nodes = num_rows*num_cols
        all_samps = np.zeros([nsamps, num_nodes])    
        for i,samp in enumerate(samps):
            for j,val in samp.iteritems():            
                all_samps[i, j] = val
                
        gibbs_mean = np.reshape(all_samps.mean(axis=0), [num_rows, num_cols])

        gibbs_means[k, :, :] = gibbs_mean

    return gibbs_means

        
def compute_all_means():
    
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt 
    
    num_rows = 7
    num_cols = 7
    
    gibbs_means = compute_gibbs_means(10)
    gibbs_mean = gibbs_means.mean(axis=0)
    gibbs_std = gibbs_means.std(axis=0)
    
    imodel_mf = create_prob54_torous(num_rows, num_cols, node_potential=0.2, edge_potential=0.2, mean_field=True)
    samps_mf = imodel_mf.sample(1, nburn=1000)
    
    samp_mf = np.zeros([num_rows*num_cols, 1])
    for j,val in samps_mf[0].iteritems():
        samp_mf[j] = val
    mf_mean = np.reshape(samp_mf, [num_rows, num_cols])

    abs_diffs = np.abs(gibbs_mean - mf_mean) 
        
    ax = plt.subplot(221)
    plt.imshow(gibbs_mean, cmap=cm.gray, interpolation='nearest', vmin=-1.0, vmax=1.0)
    plt.colorbar()
    ax.set_title('Gibbs Sampled Mean')

    ax = plt.subplot(222)
    plt.imshow(gibbs_std, cmap=cm.jet, interpolation='nearest')
    plt.colorbar()
    ax.set_title('Gibbs Standard Deviation')
    
    ax = plt.subplot(223)
    plt.imshow(mf_mean, cmap=cm.gray, interpolation='nearest', vmin=-1.0, vmax=1.0)
    plt.colorbar()
    ax.set_title('Mean-field Approximation')

    ax = plt.subplot(224)
    plt.imshow(mf_mean, cmap=cm.jet, interpolation='nearest')
    plt.colorbar()
    ax.set_title('Absolute Difference: mean=%0.2f' % abs_diffs.mean())
    
    plt.show()
    
