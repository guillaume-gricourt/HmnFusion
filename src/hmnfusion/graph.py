import copy
import itertools
from typing import Dict, List

import networkx as nx
from hmnfusion import fusion as hmn_fusion


class Graph(object):
    """A graph store all fusion as node and edges are used to built consensus from them.

    Attributes
    ----------
    g : nx.Graph
        A graph object to deal data.

    Methods
    -------
    __init__(consensus_interval=500)
        Build a Graph object and store consensus_interval value in metadata.
    """

    node_number = 0

    def __init__(self, consensus_interval=500):
        self._g = nx.Graph()
        self.update_graph_metadata({"consensus_interval": consensus_interval})

    # Getters Setters
    @property
    def graph(self) -> nx.Graph:
        return self._g

    @graph.setter
    def graph(self, g: nx.Graph) -> None:
        self._g = g

    # Others.
    def update_graph_metadata(self, data: Dict) -> None:
        """Add metadata into the graph.

        Parameters
        ----------
        data : Union[str, int]
            Data to add.

        Return
        ------
        None
        """
        for k, v in data.items():
            self.graph.graph[k] = v

    def add_node(
        self,
        fusion: hmn_fusion.Fusion,
        level: int = 0,
        is_consensus: bool = False,
        is_interest: bool = False,
    ) -> int:
        """Add a node (eg. a fusion) to the graph.

        Parameters
        ----------
        fusion : fusion.Fusion
            A fusion event
        level: int, default: 0
            The level of the fusion
              - 0: built directly from a software : Genefuse, Lumpy
              - 1: a consensus between fusions coming from a unique software
              - 2: a consensus between fusions coming from consensus previously defined
        is_consensus: bool, default: False
            If the fusion is a consensus of, at least, two fusions
        is_interest: bool, default: False
            If the fusion is enough interesting to quantify it.

        Return
        ------
        int
            Tht nth fusion in the graph
        """
        fusion = copy.deepcopy(fusion)
        self.node_number += 1
        fusion.number = self.node_number
        self.graph.add_node(
            self.node_number,
            fusion=fusion,
            level=level,
            is_consensus=is_consensus,
            is_interest=is_interest,
        )
        return self.node_number

    def add_nodes(self, fusions: List[hmn_fusion.Fusion]) -> None:
        """Add several nodes (eg. a fusion) to the graph.

        Parameters
        ----------
        fusions : List[fusion.Fusion]
            Some fusions to add

        Return
        ------
        None

        See Also
        --------
        add_node()
        """
        for fusion in fusions:
            self.add_node(fusion)

    def _add_node_consensus(self, nl, nr, level=1):
        fusion = None
        if level == 1:
            if self.graph.nodes[nl]["fusion"] > self.graph.nodes[nr]["fusion"]:
                fusion = copy.deepcopy(self.graph.nodes[nl]["fusion"])
            else:
                fusion = copy.deepcopy(self.graph.nodes[nr]["fusion"])
        else:
            fusion = copy.deepcopy(self.graph.nodes[nl]["fusion"])
        fusion.is_consensus = True
        node_fusion = self.add_node(fusion, level, True)
        self.graph.add_edge(nl, node_fusion)
        self.graph.add_edge(nr, node_fusion)
        return node_fusion

    def _grapp_subgraph(self, software=""):
        g = nx.subgraph(
            self.graph,
            [
                x
                for x in self.graph.nodes
                if self.graph.nodes[x]["fusion"].software == software
            ],
        )
        cons = [x for x in g.nodes if g.nodes[x]["is_consensus"]]
        alone = [x for x in g.nodes if self.graph.degree[x] == 0]
        return (g, cons, alone)

    def consensus_single(self) -> None:
        """Create consensus nodes from nodes already here in the graph.

        Return
        ------
        None
        """
        softwares = set(
            [self.graph.nodes[x]["fusion"].software for x in self.graph.nodes]
        )
        for software in sorted(list(softwares)):
            g, cons, alone = self._grapp_subgraph(software)
            nodes = g.nodes
            nodes_added = []
            for fl, fr in itertools.combinations(nodes, 2):
                if self.graph.nodes[fl]["fusion"].is_near(
                    self.graph.nodes[fr]["fusion"],
                    self.graph.graph["consensus_interval"],
                ):
                    if self.graph.degree[fl] > 0 or self.graph.degree[fr] > 0:
                        neighbors = list(nx.neighbors(self.graph, fl))
                        neighbors += list(nx.neighbors(self.graph, fr))
                        for nc in neighbors:
                            self.graph.add_edge(nc, fl)
                            self.graph.add_edge(nc, fr)
                            # Update.
                            if (
                                self.graph.nodes[fl]["fusion"]
                                > self.graph.nodes[nc]["fusion"]
                            ):
                                self.graph.nodes[nc]["fusion"].update(
                                    self.graph.nodes[fl]["fusion"]
                                )
                            if (
                                self.graph.nodes[fr]["fusion"]
                                > self.graph.nodes[nc]["fusion"]
                            ):
                                self.graph.nodes[nc]["fusion"].update(
                                    self.graph.nodes[fr]["fusion"]
                                )
                    else:
                        nodes_added.append(self._add_node_consensus(fl, fr))
            # Refined for circular nodes.
            while len(nodes_added) > 1:
                combs = list(itertools.combinations(nodes_added, 2))
                loop = -1
                for ix, (fl, fr) in enumerate(combs):
                    if self.graph.nodes[fl]["fusion"].is_near(
                        self.graph.nodes[fr]["fusion"],
                        self.graph.graph["consensus_interval"],
                    ):
                        node_kept = -1
                        node_rm = -1
                        neighbors = []
                        # Select node to keep.
                        if (
                            self.graph.nodes[fl]["fusion"]
                            > self.graph.nodes[fr]["fusion"]
                        ):
                            node_kept = fl
                            node_rm = fr
                        else:
                            node_kept = fr
                            node_rm = fl
                        neighbors = set(nx.neighbors(self.graph, node_rm)).difference(set(nx.neighbors(self.graph, node_kept)))
                        for nc in neighbors:
                            self.graph.add_edge(nc, node_kept)
                        nodes_added.remove(node_rm)
                        self.graph.remove_node(node_rm)
                        break
                    loop = ix
                # If no update was done
                if loop + 1 == len(combs):
                    break

    def consensus_genefuse_lumpy(self) -> None:
        """Create consensus nodes from nodes already built as consensus.

        Return
        ------
        None
        """
        g_genefuse, n_genefuse_consensus, n_genefuse_alone = self._grapp_subgraph(
            "genefuse"
        )
        # Link lumpy to genefuse consensus.
        for n_genefuse in n_genefuse_consensus:
            g_lumpy, n_lumpy_consensus, n_lumpy_alone = self._grapp_subgraph("lumpy")
            for n_lumpy in n_lumpy_consensus + n_lumpy_alone:
                if self.graph.nodes[n_genefuse]["fusion"].is_near(
                    self.graph.nodes[n_lumpy]["fusion"],
                    self.graph.graph["consensus_interval"],
                ):
                    self._add_node_consensus(n_genefuse, n_lumpy, 2)
        # Link lumpy to genefuse alone.
        for n_genefuse in n_genefuse_alone:
            g_lumpy, n_lumpy_consensus, n_lumpy_alone = self._grapp_subgraph("lumpy")
            for n_lumpy in n_lumpy_consensus + n_lumpy_alone:
                if self.graph.nodes[n_genefuse]["fusion"].is_near(
                    self.graph.nodes[n_lumpy]["fusion"],
                    self.graph.graph["consensus_interval"],
                ):
                    self._add_node_consensus(n_genefuse, n_lumpy, 2)

    def select_node_interest(self):
        return [x for x in self.graph.nodes if self.graph[x]["is_interest"]]

    def define_node_interest(self) -> int:
        """Set nodes in the graph if they are of interest

        Return
        ------
        int:
            Number of node flag as interest
        """
        g_lumpy, n_lumpy_consensus, n_lumpy_alone = self._grapp_subgraph("lumpy")
        n_lumpy_consensus = sorted(
            n_lumpy_consensus, reverse=True, key=lambda x: self.graph.nodes[x]["fusion"]
        )
        n_lumpy_alone = sorted(
            n_lumpy_alone, reverse=True, key=lambda x: self.graph.nodes[x]["fusion"]
        )
        nodes = list(self.graph.nodes)

        interest = 0
        for n in nodes:
            if self.graph.nodes[n]["level"] == 2:
                self.graph.nodes[n]["is_interest"] = True
                interest += 1
        if interest == 0:
            for n in nodes:
                if (
                    self.graph.nodes[n]["level"] == 1
                    and self.graph.nodes[n]["is_consensus"]
                    and self.graph.nodes[n]["fusion"].software == "genefuse"
                ):
                    self.graph.nodes[n]["is_interest"] = True
                    interest += 1
        if interest == 0:
            for n in nodes:
                if (
                    self.graph.nodes[n]["level"] == 0
                    and self.graph.nodes[n]["fusion"].software == "genefuse"
                ):
                    self.graph.nodes[n]["is_interest"] = True
                    interest += 1
        if interest == 0:
            if len(n_lumpy_consensus) > 0:
                self.graph.nodes[n_lumpy_consensus[0]]["is_interest"] = True
                interest += 1
        if interest == 0:
            if len(n_lumpy_alone) > 0:
                self.graph.nodes[n_lumpy_alone[0]]["is_interest"] = True
                interest += 1
        return interest

    def trim_node(self) -> None:
        """Delete fusion in the graph which have on their breakpoints
        the same chromosome.

        Return
        ------
        None
        """
        nodes = list(self.graph.nodes)
        for n in nodes:
            if self.graph.nodes[n]["fusion"].is_same_chrom():
                self.graph.remove_node(n)
        nodes = list(self.graph.nodes)
        for n in nodes:
            if self.graph.nodes[n]["is_consensus"] and self.graph.degree[n] == 0:
                self.graph.remove_node(n)

    def label_build_from(self, n) -> List[str]:
        """Get name of all nodes used to build the selected node.

        Parameters
        ----------
        n:
            id of a node selected

        Return
        ------
        List[str]
            Labels of nodes
        """
        labels = []
        if self.graph.nodes[n]["level"] > 0:
            candidates = nx.neighbors(self.graph, n)
            candidates = [
                x
                for x in candidates
                if self.graph.nodes[x]["level"] < self.graph.nodes[n]["level"]
            ]
            labels = [self.graph.nodes[x]["fusion"].get_name() for x in candidates]
        labels.sort()
        return labels

    # Import Export
    def to_dict(self) -> Dict:
        """Export the g (eg. graph) attribute as a Json
        Return
        ------
        Dict
            The graph attribute represented as a dictionary.
        """
        return nx.readwrite.json_graph.cytoscape_data(self.graph)

    @classmethod
    def from_dict(cls, data: Dict) -> "Graph":
        """Build object from a dictionary

        Parameters
        ----------
        data: Dict
            A dictionary with all data needed to build a Graph

        Return
        ------
        Graph
            A novel object.
        """
        g = Graph()
        g.graph = nx.readwrite.json_graph.cytoscape_graph(data)
        nodes = g.graph.nodes
        for n in nodes:
            g.graph.nodes[n]["fusion"] = hmn_fusion.Fusion.from_dict(
                g.graph.nodes[n]["fusion"]
            )
        return g

    # Meta functions
    def __key(self):
        return (self.graph, self.graph.graph["consensus_interval"])

    def __repr__(self) -> str:
        return "Graph %s %s" % (self.graph, self.graph.graph["consensus_interval"])

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Graph):
            return self.__key() == other.__key()
        return NotImplemented
