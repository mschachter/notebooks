"""
Created on Oct 6, 2011

@author: Mike Schachter (mike.schachter@gmail.com)
"""

import copy
import operator

import numpy as np
import networkx as nx

class UPGraph():
    """ This class represents an undirected probablistic tree graph with potential functions """

    def __init__(self):
        #self.g is the networkx graph that embodies the nodes and connections        
        self.g = nx.Graph()
        
        #self.values is a map where the key is a node id and the value is a
        #list of possible values that n can take on
        self.values = {}
        
        #self.potentials is a map, where the key is a tuple representing the
        #arguments of the potential, and the value is the actual potential
        self.potentials = {}
        
        #the normalizing constant
        self.Z = None
    
    def add_node(self, id, potential, values):
        """ Add node, given it's id, potential function, and values it can take on """
        self.g.add_node(id)
        self.potentials[id] = potential
        self.values[id] = values
    
    def connect(self, i, j, potential):
        """ Add a connection with it's pairwise potential function """
        self.g.add_edge(i, j)
        self.potentials[(i, j)] = potential
        self.potentials[(j, i)] = potential
        
    def normalizing_constant(self):
        if self.Z is not None:
            return self.Z
        
        Z = 0
        #assumes all nodes are binary valued...
        bcombs = get_all_binary_combos(len(self.g.nodes()))
        for bc in bcombs:
            #set all the node values
            node_vals = {}
            for k,n in enumerate(self.g.nodes()):
                node_vals[n] = bc[k]
            prod = 1
            #evaluate all single node potentials
            for n in self.g.nodes():
                prod *= self.potentials[n](node_vals[n])
            #evaluate all pairwise node potentials
            for (n1, n2) in self.g.edges():
                prod *= self.potentials[(n1, n2)](node_vals[n1], node_vals[n2])
            Z += prod
            
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


class Message():
    """ This class represents a message functional from node i to j. It's
        recursive and contains nested child messages that are evaluated
        when the message is evaluated. """ 
    
    def __init__(self, i, j, ivals, single_potential, pairwise_potential):
        self.single_potential = single_potential
        self.pairwise_potential = pairwise_potential
        self.vals = ivals
        self.to_node = j
        self.from_node = i        
        self.child_messages = []
    
    def add_child(self, cmsg):
        self.child_messages.append(cmsg)
    
    def evaluate(self, jval):
        #marginalize across i
        sum = 0
        for ival in self.vals:            
            prod = self.single_potential(ival)
            #BUG: because pairwise potentials are symmetric, we don't care if
            #the order of ival and jval are reversed            
            prod *= self.pairwise_potential(ival, jval)
            
            #multiply out the child messages
            for cmsg in self.child_messages:
                prod *= cmsg.evaluate(ival)
            
            sum += prod            
        return sum
    
    def __str__(self):
        rstr = '[Message: from=%d, to=%d | children=%s]' % \
               (self.from_node, self.to_node, ','.join([str(x) for x in self.child_messages]))
        return rstr

class SPNode():
    """ This class is a datastructure for a tree node, it contains the messages sent to
        and by the node, and a list of neighbors. It's meant to be a local structure that
        can exist on it's on processor in a cluster of computers.
    """
               
    def __init__(self, id, neighbors):
        self.id = id        
        self.neighbors = neighbors
        self.received_from = []
        self.sent_to = []
        self.is_leaf_node = len(neighbors) == 1
        self.messages = []
    
    def decide_recipients(self):
        """ Given the neighbors, messages sent, and messages received, this
            function determines what neighbors should be sent messages to, if any.
        """
        
        rlst = []
        if self.is_leaf_node:
            if len(self.sent_to) == 0 and len(self.received_from) == 0:
                rlst.append(self.neighbors[0])
        else:
            if len(self.sent_to) == 0:
                if len(self.received_from) == len(self.neighbors):                    
                    #send messages to all neighbors that we haven't sent messages to
                    rlst = [n for n in self.neighbors if n not in self.sent_to]
                elif len(self.received_from) == len(self.neighbors)-1:
                    #send messages to neighbor who we haven't received messages from
                    rlst = [n for n in self.neighbors if n not in self.received_from]
                        
            elif len(self.received_from) == len(self.neighbors) and \
                 len(self.sent_to) < len(self.neighbors):
                #send messages to neighbors who haven't been sent messages yet
                rlst = [n for n in self.neighbors if n not in self.sent_to]
        
        return rlst
                
        

    
class SumProduct():
    """ This class represents an instance of the SUM-PRODUCT algorithm using
        synchronous pseudo-parallel updates. """
        
    def __init__(self, upgraph):
        self.upg = upgraph        
        self.message_queue = []
        self.iteration = 0
        self.blocked = False
        self.sp_nodes = {}
        
        for n in self.upg.g.nodes():
            nbrs = self.upg.g[n].keys() #neighbors of n
            sp = SPNode(n, nbrs)
            self.sp_nodes[n] = sp
    
    
    def iterate(self):
        """ Go through each node, let each node determine which other nodes to
            send messages to, queue those messages, and then send them after
            all nodes are polled. """
        
        if self.blocked:
            print 'SUM-PRODUCT is blocked!'
            return
        
        self.iteration += 1
        print 'Iteration %d:' % self.iteration
        
        #go through each SPNode, and ask it whether
        #it has any messages to send
        for spn in self.sp_nodes.values():            
            print '\tUpdating node %d' % spn.id
            
            #get a list of messages to be sent and queue them
            rlst = spn.decide_recipients()
            if len(rlst) > 0:
                for r in rlst:
                    self.queue_message(spn.id, r)
            
            print '\t  Neighbors: %s' % spn.neighbors
            print '\t  Sent To: %s' % spn.sent_to
            print '\t  Received From: %s' % spn.received_from
            print '\t  # of Messages: %d' % len(spn.messages)
            print '\t  Targets: %s' % rlst
        
        #send all the messages in the queue.
        msgs_sent = 0
        for msg in self.message_queue:
            print 'Sending message from %d to %d' % (msg.from_node, msg.to_node)
            print msg
            self.send_message(msg)
            msgs_sent += 1
                
        self.message_queue = []        
        
        if msgs_sent == 0:
            self.blocked = True
    
    
    def queue_message(self, i, j):
        """ Queues a message function from node i to node j """
        
        #get all the messages sent to i
        spn_i = self.sp_nodes[i]
        child_msgs = spn_i.messages
        
        #append all messages sent to i to message that will be sent to j
        msg = Message(i, j, self.upg.values[i], self.upg.potentials[i], self.upg.potentials[(i, j)])
        for child_msg in child_msgs:
            if child_msg.from_node != j:
                msg.add_child(child_msg)
        
        #add the message to j's queue
        self.message_queue.append(msg)        
    
    def send_message(self, msg):
        """ Actually sends a message """
        
        spn_to = self.sp_nodes[msg.to_node]
        spn_from = self.sp_nodes[msg.from_node]
        
        spn_to.messages.append(msg)
        spn_to.received_from.append(spn_from.id)        
        spn_from.sent_to.append(spn_to.id)
    
    
    def compute_marginal(self, n):
        """ Computes the marginal of node n """
        
        spn = self.sp_nodes[n]
        vals = self.upg.values[n]        
        #take product of single node potential and messages sent to node
        prod = np.zeros([len(vals), 1])
        for k,v in enumerate(vals):
            prod[k] = self.upg.potentials[n](v)
            for msg in spn.messages:
                prod[k] *= msg.evaluate(v)
            
        return prod / self.upg.normalizing_constant()
    
    def print_marginals(self):
        """ Compute and print all marginals for the graph. """ 
        
        for n in self.upg.g.nodes():
            m = self.compute_marginal(n)
            print 'P[x%d] = (%0.2f, %0.2f)' % (n, m[0], m[1])
        
    
        
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
    