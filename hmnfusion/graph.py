import copy
import itertools

import networkx as nx

from networkx.readwrite import json_graph

class Graph():
    """Class Fusion. Provide attributes and methods to manipulate Fusion items"""
    node_number = 0

    def __init__(self, consensus_interval=500):
        """Construct Fusion object with a software name (optional)"""
        self._g = nx.Graph()
        self.update_graph_metadata({'consensus_interval':consensus_interval})

    # Getters Setters
    @property
    def graph(self):
        return self._g
    @graph.setter
    def graph(self, g):
        self._g = g

    # Others.
    def update_graph_metadata(self, data):
        for k, v in data.items():
            self.graph.graph[k] = v

    def add_node(self, fusion, level=0, is_consensus=False, is_interest=False):
        fusion = copy.deepcopy(fusion)
        self.node_number += 1
        fusion.number = self.node_number
        self.graph.add_node(self.node_number, fusion=fusion, level=level, is_consensus=is_consensus, is_interest=is_interest)
        return self.node_number

    def add_nodes(self, fusions):
        for fusion in fusions:
            self.add_node(fusion)

    def _add_node_consensus(self, nl, nr, level=1):
        fusion = None
        if level == 1:
            if self.graph.nodes[nl]['fusion'] > self.graph.nodes[nr]['fusion']:
                fusion = copy.deepcopy(self.graph.nodes[nl]['fusion'])
            else:
                fusion = copy.deepcopy(self.graph.nodes[nr]['fusion'])
        else:
            fusion = copy.deepcopy(self.graph.nodes[nl]['fusion'])
        fusion.is_consensus = True
        node_fusion = self.add_node(fusion, level, True)
        self.graph.add_edge(nl, node_fusion)
        self.graph.add_edge(nr, node_fusion)

    def _grapp_subgraph(self, software=''):
        """ """
        g = nx.subgraph(self.graph, [x for x in self.graph.nodes if self.graph.nodes[x]['fusion'].software == software])
        cons = [x for x in g.nodes if g.nodes[x]['is_consensus']]
        alone = [x for x in g.nodes if self.graph.degree[x] == 0]
        return (g, cons, alone)

    def consensus_single(self):
        """Create fusion consensus for each software"""
        softwares = set([self.graph.nodes[x]['fusion'].software for x in self.graph.nodes])
        softwares = sorted(list(softwares))
        for software in softwares:
            g, cons, alone = self._grapp_subgraph(software)
            nodes = g.nodes
            for fl, fr in itertools.combinations(nodes, 2):
                if self.graph.nodes[fl]['fusion'].is_near(self.graph.nodes[fr]['fusion'], self.graph.graph['consensus_interval']):
                    if self.graph.degree[fl] > 0 or self.graph.degree[fr] > 0:
                        for nc in list(nx.neighbors(self.graph, fl)) + list(nx.neighbors(self.graph, fr)):
                            self.graph.add_edge(nc, fl)
                            self.graph.add_edge(nc, fr)
                            # Update.
                            if self.graph.nodes[fl]['fusion'] > self.graph.nodes[nc]['fusion']:
                                self.graph.nodes[nc]['fusion'].update(self.graph.nodes[fl]['fusion'])
                            if self.graph.nodes[fr]['fusion'] > self.graph.nodes[nc]['fusion']:
                                self.graph.nodes[nc]['fusion'].update(self.graph.nodes[fr]['fusion'])
                    else:
                        self._add_node_consensus(fl, fr)

    def consensus_genefuse_lumpy(self):
        """ """
        g_genefuse, n_genefuse_consensus, n_genefuse_alone = self._grapp_subgraph('genefuse')
        # Link lumpy to genefuse consensus.
        for n_genefuse in n_genefuse_consensus:
            g_lumpy, n_lumpy_consensus, n_lumpy_alone = self._grapp_subgraph('lumpy')
            for n_lumpy in n_lumpy_consensus + n_lumpy_alone:
                if self.graph.nodes[n_genefuse]['fusion'].is_near(self.graph.nodes[n_lumpy]['fusion'], self.graph.graph['consensus_interval']):
                    self._add_node_consensus(n_genefuse, n_lumpy, 2)
        # Link lumpy to genefuse alone.
        for n_genefuse in n_genefuse_alone:
            g_lumpy, n_lumpy_consensus, n_lumpy_alone = self._grapp_subgraph('lumpy')
            for n_lumpy in n_lumpy_consensus + n_lumpy_alone:
                if self.graph.nodes[n_genefuse]['fusion'].is_near(self.graph.nodes[n_lumpy]['fusion'], self.graph.graph['consensus_interval']):
                    self._add_node_consensus(n_genefuse, n_lumpy, 2)

    def select_node_interest(self):
        return [x for x in self.graph.nodes if self.graph[x]['is_interest']]

    def define_node_interest(self):
        """ """
        g_lumpy, n_lumpy_consensus, n_lumpy_alone = self._grapp_subgraph('lumpy')
        n_lumpy_consensus = sorted(n_lumpy_consensus, reverse=True, key=lambda x: self.graph.nodes[x]['fusion'])
        n_lumpy_alone = sorted(n_lumpy_alone, reverse=True, key=lambda x: self.graph.nodes[x]['fusion'])
        nodes = list(self.graph.nodes)

        interest = 0
        for n in nodes:
            if self.graph.nodes[n]['level'] == 2:
                self.graph.nodes[n]['is_interest'] = True
                interest += 1
        if interest == 0:
            for n in nodes:
                if self.graph.nodes[n]['level'] == 1 and self.graph.nodes[n]['is_consensus'] and self.graph.nodes[n]['fusion'].software == 'genefuse':
                    self.graph.nodes[n]['is_interest'] = True
                    interest += 1
        if interest == 0:
            for n in nodes:
                if self.graph.nodes[n]['level'] == 0 and self.graph.nodes[n]['fusion'].software == 'genefuse':
                    self.graph.nodes[n]['is_interest'] = True
                    interest += 1
        if interest == 0:
            if len(n_lumpy_consensus) > 0:
                self.graph.nodes[n_lumpy_consensus[0]]['is_interest'] = True
                interest += 1
        if interest == 0:
            if len(n_lumpy_alone) > 0:
                self.graph.nodes[n_lumpy_alone[0]]['is_interest'] = True
                interest += 1
        return interest

    def trim_node(self):
        """Delete fusion in the graph which have on their breakpoints the same chromosome."""
        nodes = list(self.graph.nodes)
        for n in nodes:
            if self.graph.nodes[n]['fusion'].is_same_chrom():
                self.graph.remove_node(n)
        nodes = list(self.graph.nodes)
        for n in nodes:
            if self.graph.nodes[n]['is_consensus'] and self.graph.degree[n] == 0:
                self.graph.remove_node(n)

    def label_build_from(self, n):
        labels = []
        if self.graph.nodes[n]['level'] > 0:
            candidates = nx.neighbors(self.graph, n)
            candidates = [x for x in candidates if self.graph.nodes[x]['level'] < self.graph.nodes[n]['level']]
            labels = [self.graph.nodes[x]['fusion'].get_name() for x in candidates]
        labels.sort()
        return labels

    # Import Export
    def to_dict(self):
        """Export Fusion as dict"""
        return nx.readwrite.json_graph.cytoscape_data(self.graph)

    @classmethod
    def from_dict(cls, data):
        """Construct a Fusion object from a dict"""
        g = Graph()
        g.graph = nx.readwrite.json_graph.cytoscape_graph(data)
        return g

    # Meta functions
    def __key(self):
        return (self._graph, self.graph.graph['consensus_interval'])

    def __repr__(self):
        return 'Graph %s %s'%(self._graph, self.graph.graph['consensus_interval'])

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        if isinstance(other, A):
            return self.__key() == other.__key()
        return NotImplemented
