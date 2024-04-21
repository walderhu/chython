import logging


class Node:
    index = 0
    def __init__(self):
        self.atom = None
        self.neighbors = []
        self.is_visited = False
        self.parent = None
        self.mate = None
        self.index = Node.index
        Node.index += 1

    def __repr__(self):
        return str(self.atom)

    def set_atom(self, atom):
        self.atom = atom

class Path:

    def __init__(self):
        self.nodes = []

class Match:

    def __init__(self, nodes):
        self.nodes = nodes
        self.freenodes = []
        for node in nodes:
            self.freenodes.append(node)
        self.supernodes = []

    @staticmethod
    def from_structure(structure):
        pi_subgraph = structure.find_pi_subgraph()
        node_nr = len(pi_subgraph.keys())

        atom_to_node = {}
        nodes = [Node() for i in range(node_nr)]

        for i, atom in enumerate(pi_subgraph):
            atom_to_node[atom] = i

        for atom, neighbours in pi_subgraph.items():
            index_1 = atom_to_node[atom]
            nodes[index_1].set_atom(atom)

            for neighbour in neighbours:
                index_2 = atom_to_node[neighbour]
                nodes[index_1].neighbors.append(nodes[index_2])

        return Match(nodes)

    def clear_nodes(self):
        for node in self.nodes:
            node.is_visited = False
            node.parent = None

    def find_augmenting_path(self, root):
        self.clear_nodes()
        queue = [root]
        while len(queue) > 0:
            cur_node = queue.pop(0)
            cur_node.is_visited = True
            for node in cur_node.neighbors:
                if node == cur_node.parent:
                    continue

                elif node.is_visited:
                    cycle = self.find_cycles(node, cur_node)
                    if len(cycle) % 2 == 1:
                        logging.debug('blossom: {}'.format(cycle))
                        snode = self.shrink_blossom(cycle)
                        self.supernodes.append(snode)
                        for v in cycle:
                            if v in queue:
                                queue.remove(v)
                            if v.is_visited:
                                snode.is_visited = True
                                snode.parent = v.parent
                        queue.append(snode)
                        break

                elif node.mate is None:
                    node.parent = cur_node
                    return self.construct_augmenting_path(node)

                else:
                    node.is_visited = True
                    node.mate.is_visited = True
                    node.parent = cur_node
                    node.mate.parent = node
                    queue.append(node.mate)
        raise Exception('cannot find an augmenting path')

    def unmatched_nodes(self):
        self.maximum_matching()

        count = 0
        for node in self.nodes:
            if node.mate is not None:
                count += 1

        return len(self.nodes) - count

    def maximum_matching(self):
        while len(self.freenodes) > 0:
            logging.debug('freenodes: {}'.format(self.freenodes))

            for node in self.freenodes:
                try:
                    path = self.find_augmenting_path(node)
                    logging.debug('augmenting path: {}'.format(path.nodes))
                    self.invert_path(path)
                    self.freenodes.remove(path.nodes[0])
                    self.freenodes.remove(path.nodes[-1])
                    break
                except Exception as e:
                    logging.info(e)
            else:
                logging.info('Tried all free nodes, no more augmenting path.')

                break

    def invert_path(self, path):
        assert len(path.nodes) % 2 == 0
        for i in range(0, len(path.nodes), 2):
            path.nodes[i].mate = path.nodes[i + 1]
            path.nodes[i + 1].mate = path.nodes[i]

    def construct_augmenting_path(self, node):
        path = Path()
        path.nodes.append(node)
        node = node.parent
        path.nodes.append(node)
        while node.mate is not None:
            node = node.parent
            path.nodes.append(node)

        while len(self.supernodes) > 0:
            snode = self.supernodes.pop()
            self.expand_supernode(snode)
            if snode == path.head():
                path.replace_head()
            elif snode == path.tail():
                path.replace_tail()

        while path.nodes[0].mate is not None:
            path.nodes.insert(path.nodes[0].parent, 0)

        while path.nodes[-1].mate is not None:
            path.nodes.append(path.nodes[-1].parent)

        return path