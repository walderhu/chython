# A dependency-free version of networkx's implementation of Johnson's cycle finding algorithm
# Original implementation: https://github.com/networkx/networkx/blob/master/networkx/algorithms/cycles.py#L109
# Original paper: Donald B Johnson. "Finding all the elementary circuits of a directed graph." SIAM Journal on Computing. 1975.

from collections import OrderedDict, defaultdict
import copy

def simple_cycles(G):
    # Yield every elementary cycle in python graph G exactly once
    # Expects a dictionary mapping from vertices to iterables of vertices
    def _unblock(thisnode, blocked, B):
        stack = {thisnode}
        while stack:
            node = stack.pop()
            if node in blocked:
                blocked.remove(node)
                stack.update(B[node])
                B[node].clear()
    G = {v: set(nbrs) for (v, nbrs) in G.items()} # make a copy of the graph
    sccs = strongly_connected_components(G)
    while sccs:
        scc = sccs.pop()
        startnode = scc.pop()
        path = [startnode]
        blocked = set()
        closed = set()
        blocked.add(startnode)
        B = defaultdict(set)
        stack = [ (startnode, list(G[startnode])) ]
        while stack:
            thisnode, nbrs = stack[-1]
            if nbrs:
                nextnode = nbrs.pop()
                if nextnode == startnode:
                    yield path[:]
                    closed.update(path)
                elif nextnode not in blocked:
                    path.append(nextnode)
                    stack.append( (nextnode,list(G[nextnode])) )
                    closed.discard(nextnode)
                    blocked.add(nextnode)
                    continue
            if not nbrs:
                if thisnode in closed:
                    _unblock(thisnode,blocked,B)
                else:
                    for nbr in G[thisnode]:
                        if thisnode not in B[nbr]:
                            B[nbr].add(thisnode)
                stack.pop()
                path.pop()
        remove_node(G, startnode)
        H = subgraph(G, set(scc))
        sccs.extend(strongly_connected_components(H))


def strongly_connected_components(graph):
    # Tarjan's algorithm for finding SCC's
    # Robert Tarjan. "Depth-first search and linear graph algorithms." SIAM journal on computing. 1972.
    # Code by Dries Verdegem, November 2012
    # Downloaded from http://www.logarithmic.net/pfh/blog/01208083168

    index_counter = [0]
    stack = []
    lowlink = {}
    index = {}
    result = []
    
    def _strong_connect(node):
        index[node] = index_counter[0]
        lowlink[node] = index_counter[0]
        index_counter[0] += 1
        stack.append(node)
    
        successors = graph[node]
        for successor in successors:
            if successor not in index:
                _strong_connect(successor)
                lowlink[node] = min(lowlink[node],lowlink[successor])
            elif successor in stack:
                lowlink[node] = min(lowlink[node],index[successor])

        if lowlink[node] == index[node]:
            connected_component = []

            while True:
                successor = stack.pop()
                connected_component.append(successor)
                if successor == node: break
            result.append(connected_component[:])
    
    for node in graph:
        if node not in index:
            _strong_connect(node)
    
    return result


def remove_node(G, target):
    # Completely remove a node from the graph
    # Expects values of G to be sets
    del G[target]
    for nbrs in G.values():
        nbrs.discard(target)


def subgraph(G, vertices):
    # Get the subgraph of G induced by set vertices
    # Expects values of G to be sets
    return {v: G[v] & vertices for v in vertices}


class Cycles:
    def __init__(self, structure):
        self.find_unique_cycles(structure)
        self.make_microcycle_graph()

    def find_sssr(self):
        sorted_cycles = sorted(self.all_cycles, key=lambda x: len(x))
        sssr = []
        atoms = set()
        for cycle in sorted_cycles:
            add_cycle = False
            for atom in cycle:

                if atom not in atoms:
                    add_cycle = True
                    break
            if add_cycle:
                sssr.append(cycle)
                
            for atom in cycle:
                atoms.add(atom)

        return sssr

    def find_unique_cycles(self, structure):
        all_cycles = simple_cycles(structure.graph)
        self.all_cycles = []

        unique_cycles = set()
        for cycle in all_cycles:

            if len(cycle) > 2:
                self.all_cycles.append(cycle)
                cycle_components = sorted(cycle, key=lambda x: x.nr)
                cycle_components = tuple(cycle_components)
                if len(cycle_components) < 10:
                    unique_cycles.add(cycle_components)
        
        self.unique_cycles = unique_cycles

    def make_microcycle_graph(self):
        self.graph = {}

        for cycle_1 in self.unique_cycles:
            if len(cycle_1) < 10:
                self.graph[cycle_1] = []
            for cycle_2 in self.unique_cycles:
                if cycle_1 != cycle_2:
                    if len(cycle_2) < 10:
                        if len(set(cycle_1).intersection(set(cycle_2))) > 1:
                            self.graph[cycle_1].append(cycle_2)