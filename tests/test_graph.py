import unittest

import matplotlib.pyplot as plt
import networkx as nx
from hmnfusion import fusion, graph


class TestGraph(unittest.TestCase):
    """Class to test Graph class"""

    def setUp(self):
        """Initialize objects"""
        # Evidence.
        # 0
        self.e1 = dict(raw=0, split=0, mate=0, clipped=0, depth=0)
        # 30
        self.e2 = dict(raw=3, split=10, mate=10, clipped=10, depth=100)
        # 60
        self.e3 = dict(raw=6, split=50, mate=0, clipped=10, depth=200)
        # 300
        self.e4 = dict(raw=30, split=100, mate=100, clipped=100, depth=300)
        # 500
        self.e5 = dict(raw=50, split=200, mate=200, clipped=100, depth=1000)
        # 600
        self.e6 = dict(raw=60, split=350, mate=150, clipped=100, depth=1000)
        # 700
        self.e7 = dict(raw=70, split=100, mate=300, clipped=300, depth=1000)
        # 800
        self.e8 = dict(raw=80, split=200, mate=250, clipped=350, depth=1000)
        # 900
        self.e9 = dict(raw=90, split=500, mate=300, clipped=100, depth=1000)

        # Region - Cons.
        self.r1 = dict(chrom="chr9", orientation="left", position=500)
        self.r2 = dict(chrom="chr22", orientation="rigth", position=10000)
        self.r3 = dict(chrom="chr9", orientation="left", position=1200)
        self.r4 = dict(chrom="chr22", orientation="right", position=10200)
        self.r5 = dict(chrom="chr9", orientation="rigth", position=600)
        self.r6 = dict(chrom="chr22", orientation="left", position=10400)
        self.r7 = dict(chrom="chr9", orientation="right", position=800)
        self.r8 = dict(chrom="chr22", orientation="left", position=9500)

        self.r9 = dict(chrom="chr10", orientation="left", position=400)
        self.r10 = dict(chrom="chr11", orientation="right", position=600)
        self.r11 = dict(chrom="chr10", orientation="left", position=500)
        self.r12 = dict(chrom="chr11", orientation="right", position=1000)

        # Region - Alone.
        self.r13 = dict(chrom="chr10", orientation="right", position=4000)
        self.r14 = dict(chrom="chr11", orientation="left", position=6000)
        self.r15 = dict(chrom="chr12", orientation="left", position=400)
        self.r16 = dict(chrom="chr13", orientation="right", position=600)
        self.r17 = dict(chrom="chr12", orientation="right", position=4000)
        self.r18 = dict(chrom="chr13", orientation="left", position=6000)

        # Region - Same chr.
        self.r19 = dict(chrom="chr10", orientation="left", position=400)
        self.r20 = dict(chrom="chr10", orientation="right", position=6000)

        # Region - Tricky
        self.r21 = dict(chrom="chr9", orientation="left", position=1760)
        self.r22 = dict(chrom="chr22", orientation="right", position=2004)
        self.r23 = dict(chrom="chr9", orientation="left", position=1100)
        self.r24 = dict(chrom="chr22", orientation="right", position=2736)
        self.r25 = dict(chrom="chr9", orientation="left", position=1242)
        self.r26 = dict(chrom="chr22", orientation="right", position=2594)
        self.r27 = dict(chrom="chr9", orientation="left", position=1616)
        self.r28 = dict(chrom="chr22", orientation="right", position=2146)

        # Fusion.
        self.f1 = fusion.Fusion()
        # Same chr
        self.f10g = fusion.Fusion.from_dict(
            dict(first=self.r19, second=self.r20, evidence=self.e2, software="genefuse")
        )
        self.f10l = fusion.Fusion.from_dict(
            dict(first=self.r20, second=self.r19, evidence=self.e4, software="lumpy")
        )
        # Cons
        self.f1g = fusion.Fusion.from_dict(
            dict(first=self.r1, second=self.r2, evidence=self.e2, software="genefuse")
        )
        self.f2g = fusion.Fusion.from_dict(
            dict(first=self.r3, second=self.r4, evidence=self.e3, software="genefuse")
        )
        self.f3g = fusion.Fusion.from_dict(
            dict(first=self.r6, second=self.r5, evidence=self.e6, software="genefuse")
        )
        self.f4g = fusion.Fusion.from_dict(
            dict(first=self.r7, second=self.r8, evidence=self.e5, software="genefuse")
        )
        self.f5g = fusion.Fusion.from_dict(
            dict(first=self.r9, second=self.r10, evidence=self.e6, software="genefuse")
        )
        # Alone
        self.f6g = fusion.Fusion.from_dict(
            dict(first=self.r12, second=self.r11, evidence=self.e8, software="genefuse")
        )
        self.f7g = fusion.Fusion.from_dict(
            dict(first=self.r13, second=self.r14, evidence=self.e7, software="genefuse")
        )
        self.f8g = fusion.Fusion.from_dict(
            dict(first=self.r15, second=self.r16, evidence=self.e9, software="genefuse")
        )
        self.f9g = fusion.Fusion.from_dict(
            dict(first=self.r18, second=self.r17, evidence=self.e5, software="genefuse")
        )
        # Cons
        self.f1l = fusion.Fusion.from_dict(
            dict(first=self.r1, second=self.r2, evidence=self.e2, software="lumpy")
        )
        self.f2l = fusion.Fusion.from_dict(
            dict(first=self.r3, second=self.r4, evidence=self.e3, software="lumpy")
        )
        self.f3l = fusion.Fusion.from_dict(
            dict(first=self.r6, second=self.r5, evidence=self.e6, software="lumpy")
        )
        self.f4l = fusion.Fusion.from_dict(
            dict(first=self.r7, second=self.r8, evidence=self.e5, software="lumpy")
        )
        self.f5l = fusion.Fusion.from_dict(
            dict(first=self.r9, second=self.r10, evidence=self.e6, software="lumpy")
        )
        # Alone
        self.f6l = fusion.Fusion.from_dict(
            dict(first=self.r12, second=self.r11, evidence=self.e8, software="lumpy")
        )
        self.f7l = fusion.Fusion.from_dict(
            dict(first=self.r13, second=self.r14, evidence=self.e7, software="lumpy")
        )
        self.f8l = fusion.Fusion.from_dict(
            dict(first=self.r15, second=self.r16, evidence=self.e9, software="lumpy")
        )
        self.f9l = fusion.Fusion.from_dict(
            dict(first=self.r18, second=self.r17, evidence=self.e9, software="lumpy")
        )
        # Tricky
        self.f11l = fusion.Fusion.from_dict(
            dict(first=self.r21, second=self.r22, evidence=self.e2, software="lumpy")
        )
        self.f12l = fusion.Fusion.from_dict(
            dict(first=self.r23, second=self.r24, evidence=self.e3, software="lumpy")
        )
        self.f13l = fusion.Fusion.from_dict(
            dict(first=self.r25, second=self.r26, evidence=self.e4, software="lumpy")
        )
        self.f14l = fusion.Fusion.from_dict(
            dict(first=self.r27, second=self.r28, evidence=self.e5, software="lumpy")
        )

        """
        --------
        - Plan -
        --------

        genefuse, lumpy
          - 1 node
            1a : 1 alone (f1g)
          - 2 nodes
            2a : 2 alone (f1g, f5g)
            2b : 2 -> 1 cons (f1g, f2g)
          - 3 nodes
            3a : 3 alone (f1g, f7g, f8g)
            3b : 1 alone, 2 -> 1 cons (f7g, f3g, f4g)
            3c : 3 -> 1 cons (f2g, f3g, f4g)
          - 4 nodes
            4a : 2 alone, 2 -> 1 cons (f3g, f4g, f7g, f8g)
            4b : 2 -> 1 cons, 2 -> 1 cons (f1g, f6g, f4g, f5g)
            4c : 1 alone, 3 -> 1 cons (f6g, f2g, f3g, f4g)
            4d : 4 -> 1 cons (f1g, f2g, f3g, f4g)
            4e : 4 alone (f6g, f7g, f8g, f9g)
        """

        # Graph.
        self.g1 = graph.Graph()
        self.g2 = graph.Graph(0)
        self.nx1 = nx.Graph([(1, 2)])
        self.nx2 = nx.Graph([(4, "a")])

        # Build graph genefuse.
        self.g1ag = graph.Graph()
        self.g2ag = graph.Graph()
        self.g2bg = graph.Graph()
        self.g3ag = graph.Graph()
        self.g3bg = graph.Graph()
        self.g3cg = graph.Graph()
        self.g4ag = graph.Graph()
        self.g4bg = graph.Graph()
        self.g4cg = graph.Graph()
        self.g4dg = graph.Graph()
        self.g4eg = graph.Graph()

        self.g1ag.add_node(self.f1g)
        self.g2ag.add_nodes([self.f1g, self.f5g])
        self.g2bg.add_nodes([self.f1g, self.f2g])
        self.g3ag.add_nodes([self.f1g, self.f7g, self.f8g])
        self.g3bg.add_nodes([self.f7g, self.f3g, self.f4g])
        self.g3cg.add_nodes([self.f2g, self.f3g, self.f4g])
        self.g4ag.add_nodes([self.f3g, self.f4g, self.f7g, self.f8g])
        self.g4bg.add_nodes([self.f1g, self.f6g, self.f4g, self.f5g])
        self.g4cg.add_nodes([self.f6g, self.f2g, self.f3g, self.f4g])
        self.g4dg.add_nodes([self.f1g, self.f2g, self.f3g, self.f4g])
        self.g4eg.add_nodes([self.f6g, self.f7g, self.f8g, self.f9g])
        self.gg = [
            self.g1ag,
            self.g2ag,
            self.g2bg,
            self.g3ag,
            self.g3bg,
            self.g3cg,
            self.g4ag,
            self.g4bg,
            self.g4cg,
            self.g4dg,
            self.g4eg,
        ]

        # Build graph lumpy.
        self.g1al = graph.Graph()
        self.g2al = graph.Graph()
        self.g2bl = graph.Graph()
        self.g3al = graph.Graph()
        self.g3bl = graph.Graph()
        self.g3cl = graph.Graph()
        self.g4al = graph.Graph()
        self.g4bl = graph.Graph()
        self.g4cl = graph.Graph()
        self.g4dl = graph.Graph()
        self.g4el = graph.Graph()
        self.g4fl = graph.Graph()

        self.g1al.add_node(self.f1l)
        self.g2al.add_nodes([self.f1l, self.f5l])
        self.g2bl.add_nodes([self.f1l, self.f2l])
        self.g3al.add_nodes([self.f1l, self.f7l, self.f8l])
        self.g3bl.add_nodes([self.f7l, self.f3l, self.f4l])
        self.g3cl.add_nodes([self.f2l, self.f3l, self.f4l])
        self.g4al.add_nodes([self.f3l, self.f4l, self.f7l, self.f8l])
        self.g4bl.add_nodes([self.f1l, self.f6l, self.f4l, self.f5l])
        self.g4cl.add_nodes([self.f6l, self.f2l, self.f3l, self.f4l])
        self.g4dl.add_nodes([self.f1l, self.f2l, self.f3l, self.f4l])
        self.g4el.add_nodes([self.f6l, self.f7l, self.f8l, self.f9l])
        self.g4fl.add_nodes([self.f11l, self.f12l, self.f13l, self.f14l])
        self.gl = [
            self.g1al,
            self.g2al,
            self.g2bl,
            self.g3al,
            self.g3bl,
            self.g3cl,
            self.g4al,
            self.g4bl,
            self.g4cl,
            self.g4dl,
            self.g4el,
            self.g4fl,
        ]

        # Build graph genefuse-lumpy.
        self.g1agl = graph.Graph()
        self.g2agl = graph.Graph()
        self.g2bgl = graph.Graph()
        self.g3agl = graph.Graph()
        self.g3bgl = graph.Graph()
        self.g3cgl = graph.Graph()
        self.g4agl = graph.Graph()
        self.g4bgl = graph.Graph()
        self.g4cgl = graph.Graph()
        self.g4dgl = graph.Graph()
        self.g4egl = graph.Graph()

        self.g1agl.add_nodes([self.f1g, self.f1l])
        self.g2agl.add_nodes([self.f1g, self.f5g, self.f1l, self.f5l])
        self.g2bgl.add_nodes([self.f1g, self.f2g, self.f1l, self.f2l])
        self.g3agl.add_nodes(
            [self.f1g, self.f7g, self.f8g, self.f1l, self.f7l, self.f8l]
        )
        self.g3bgl.add_nodes(
            [self.f7g, self.f3g, self.f4g, self.f7l, self.f3l, self.f4l]
        )
        self.g3cgl.add_nodes(
            [self.f2g, self.f3g, self.f4g, self.f2l, self.f3l, self.f4l]
        )
        self.g4agl.add_nodes(
            [
                self.f3g,
                self.f4g,
                self.f7g,
                self.f8g,
                self.f3l,
                self.f4l,
                self.f7l,
                self.f8l,
            ]
        )
        self.g4bgl.add_nodes(
            [
                self.f1g,
                self.f6g,
                self.f4g,
                self.f5g,
                self.f1l,
                self.f6l,
                self.f4l,
                self.f5l,
            ]
        )
        self.g4cgl.add_nodes(
            [
                self.f6g,
                self.f2g,
                self.f3g,
                self.f4g,
                self.f6l,
                self.f2l,
                self.f3l,
                self.f4l,
            ]
        )
        self.g4dgl.add_nodes(
            [
                self.f1g,
                self.f2g,
                self.f3g,
                self.f4g,
                self.f1l,
                self.f2l,
                self.f3l,
                self.f4l,
            ]
        )
        self.g4egl.add_nodes(
            [
                self.f6g,
                self.f7g,
                self.f8g,
                self.f9g,
                self.f6l,
                self.f7l,
                self.f8l,
                self.f9l,
            ]
        )
        self.ggl = [
            self.g1agl,
            self.g2agl,
            self.g2bgl,
            self.g3agl,
            self.g3bgl,
            self.g3cgl,
            self.g4agl,
            self.g4bgl,
            self.g4cgl,
            self.g4dgl,
            self.g4egl,
        ]

        # Build graph genefuse-lumpy.
        self.g1at = graph.Graph()
        self.g1at.add_nodes(
            [
                self.f6g,
                self.f7g,
                self.f8g,
                self.f9g,
                self.f6l,
                self.f7l,
                self.f8l,
                self.f9l,
                self.f10g,
                self.f10l,
            ]
        )
        self.gt = [self.g1at]

        # All graphs.
        self.graph_all = self.gg + self.gl + self.ggl + self.gt

    def helper_build_consensus_single(self):
        for g in self.graph_all:
            g.consensus_single()

    def helper_build_consensus_genefuse_lumpy(self):
        for g in self.graph_all:
            g.consensus_genefuse_lumpy()

    def helper_build_define_node_interest(self):
        self.helper_build_consensus_single()
        self.helper_build_consensus_genefuse_lumpy()
        for g in self.graph_all:
            g.define_node_interest()

    def helper_check_quanti(self, g, n_nodes, n_edges, n_lev0, n_lev1, n_lev2):
        """Helper fuction for test_consensus"""
        self.assertEqual(g.graph.number_of_nodes(), n_nodes)
        self.assertEqual(g.graph.number_of_edges(), n_edges)
        self.assertEqual(
            len([x for x in g.graph.nodes if g.graph.nodes[x]["level"] == 0]), n_lev0
        )

        self.assertEqual(
            len([x for x in g.graph.nodes if g.graph.nodes[x]["level"] == 1]), n_lev1
        )
        self.assertEqual(
            len([x for x in g.graph.nodes if g.graph.nodes[x]["level"] == 2]), n_lev2
        )

    def helper_check_quali(self, g, cons_lev1=[], cons_lev2=[]):
        cons_lev1 = set(cons_lev1)
        th_lev1 = set(
            g.graph.nodes[x]["fusion"]
            for x in g.graph.nodes
            if g.graph.nodes[x]["level"] == 1
        )
        cons_lev2 = set(cons_lev2)
        th_lev2 = set(
            g.graph.nodes[x]["fusion"]
            for x in g.graph.nodes
            if g.graph.nodes[x]["level"] == 2
        )
        self.assertEqual(cons_lev1, th_lev1)
        self.assertEqual(cons_lev2, th_lev2)

    def helper_check_interest(self, g, n_lev0, n_lev1, n_lev2, fusions):
        """Helper fuction for test_consensus"""
        self.assertEqual(
            len(
                [
                    x
                    for x in g.graph.nodes
                    if g.graph.nodes[x]["level"] == 0
                    and g.graph.nodes[x]["is_interest"]
                ]
            ),
            n_lev0,
        )
        self.assertEqual(
            len(
                [
                    x
                    for x in g.graph.nodes
                    if g.graph.nodes[x]["level"] == 1
                    and g.graph.nodes[x]["is_interest"]
                ]
            ),
            n_lev1,
        )
        self.assertEqual(
            len(
                [
                    x
                    for x in g.graph.nodes
                    if g.graph.nodes[x]["level"] == 2
                    and g.graph.nodes[x]["is_interest"]
                ]
            ),
            n_lev2,
        )
        fusions = set(fusions)
        th_fusions = set(
            g.graph.nodes[x]["fusion"]
            for x in g.graph.nodes
            if g.graph.nodes[x]["is_interest"]
        )
        self.assertEqual(fusions, th_fusions)

    def helper_plot_graph(self, g, path=None):
        plt.subplot()
        nx.draw(g.graph, with_labels=True, font_weight="bold")
        if path:
            plt.savefig(path)
        else:
            plt.show()

    def test_init(self):
        """Check initialize attributes values"""
        self.assertTrue(nx.is_empty(self.g1.graph))

    def test_getters(self):
        """Test getters attributes"""
        self.assertTrue(nx.is_empty(self.g1.graph))
        self.assertTrue(nx.is_empty(self.g2.graph))

    def test_setters(self):
        """Test setters attributes"""
        # Simple
        self.g1.graph = self.nx1
        self.assertEqual(self.g1.graph, self.nx1)
        self.g2.graph = self.nx2
        self.assertEqual(self.g2.graph, self.nx2)

    def test_update_graph_metadata(self):
        """Test update_graph_metadata(self, data)"""
        self.g1.update_graph_metadata({"consensus_interval": 0})
        self.assertEqual(self.g1.graph.graph["consensus_interval"], 0)
        self.g1.update_graph_metadata({"consensus_interval": 100, "test": True})
        self.assertEqual(self.g1.graph.graph["consensus_interval"], 100)
        self.assertTrue(self.g1.graph.graph["test"])

    def test_add_node(self):
        """Test add_node(self, fusion)"""
        self.assertEqual(self.g2ag.node_number, 2)
        self.assertTrue(
            self.f1g
            in [self.g2ag.graph.nodes[x]["fusion"] for x in self.g2ag.graph.nodes]
        )
        self.assertFalse(
            self.f9g
            in [self.g2ag.graph.nodes[x]["fusion"] for x in self.g2ag.graph.nodes]
        )
        self.assertEqual(self.g4ag.node_number, 4)
        self.assertTrue(
            self.f3g
            in [self.g4ag.graph.nodes[x]["fusion"] for x in self.g4ag.graph.nodes]
        )
        self.assertFalse(
            self.f9g
            in [self.g4ag.graph.nodes[x]["fusion"] for x in self.g4ag.graph.nodes]
        )

    def test_consensus_single(self):
        """Test consensus_single(self)"""
        # single
        self.helper_build_consensus_single()

        self.helper_plot_graph(self.g4fl, "/tmp/graph.png")
        self.helper_check_quanti(self.g1ag, 1, 0, 1, 0, 0)
        self.helper_check_quanti(self.g2ag, 2, 0, 2, 0, 0)
        self.helper_check_quanti(self.g2bg, 3, 2, 2, 1, 0)
        # --
        self.helper_check_quanti(self.g3ag, 3, 0, 3, 0, 0)
        self.helper_check_quanti(self.g3bg, 4, 2, 3, 1, 0)
        self.helper_check_quanti(self.g3cg, 4, 3, 3, 1, 0)
        # --
        self.helper_check_quanti(self.g4ag, 5, 2, 4, 1, 0)
        self.helper_check_quanti(self.g4bg, 6, 4, 4, 2, 0)
        self.helper_check_quanti(self.g4cg, 5, 3, 4, 1, 0)
        # --
        self.helper_check_quanti(self.g4dg, 5, 4, 4, 1, 0)
        self.helper_check_quanti(self.g4eg, 4, 0, 4, 0, 0)

        self.helper_check_quanti(self.g1al, 1, 0, 1, 0, 0)
        self.helper_check_quanti(self.g2al, 2, 0, 2, 0, 0)
        self.helper_check_quanti(self.g2bl, 3, 2, 2, 1, 0)
        # --
        self.helper_check_quanti(self.g3al, 3, 0, 3, 0, 0)
        self.helper_check_quanti(self.g3bl, 4, 2, 3, 1, 0)
        self.helper_check_quanti(self.g3cl, 4, 3, 3, 1, 0)
        # --
        self.helper_check_quanti(self.g4al, 5, 2, 4, 1, 0)
        self.helper_check_quanti(self.g4bl, 6, 4, 4, 2, 0)
        self.helper_check_quanti(self.g4cl, 5, 3, 4, 1, 0)
        # --
        self.helper_check_quanti(self.g4dl, 5, 4, 4, 1, 0)
        self.helper_check_quanti(self.g4el, 4, 0, 4, 0, 0)
        self.helper_check_quanti(self.g4fl, 5, 4, 4, 1, 0)

        self.helper_check_quanti(self.g1agl, 2, 0, 2, 0, 0)
        self.helper_check_quanti(self.g2agl, 4, 0, 4, 0, 0)
        self.helper_check_quanti(self.g2bgl, 6, 4, 4, 2, 0)
        # --
        self.helper_check_quanti(self.g3agl, 6, 0, 6, 0, 0)
        self.helper_check_quanti(self.g3bgl, 8, 4, 6, 2, 0)
        self.helper_check_quanti(self.g3cgl, 8, 6, 6, 2, 0)
        # --
        self.helper_check_quanti(self.g4agl, 10, 4, 8, 2, 0)
        self.helper_check_quanti(self.g4bgl, 12, 8, 8, 4, 0)
        self.helper_check_quanti(self.g4cgl, 10, 6, 8, 2, 0)
        # --
        self.helper_check_quanti(self.g4dgl, 10, 8, 8, 2, 0)
        self.helper_check_quanti(self.g4egl, 8, 0, 8, 0, 0)

        self.helper_check_quali(self.g2bg, [self.f2g])
        self.helper_check_quali(self.g3bg, [self.f3g])
        self.helper_check_quali(self.g3cg, [self.f3g])
        self.helper_check_quali(self.g4ag, [self.f3g])
        self.helper_check_quali(self.g4bg, [self.f4g, self.f6g])
        self.helper_check_quali(self.g4cg, [self.f3g])
        self.helper_check_quali(self.g4dg, [self.f3g])

        self.helper_check_quali(self.g2bl, [self.f2l])
        self.helper_check_quali(self.g3bl, [self.f3l])
        self.helper_check_quali(self.g3cl, [self.f3l])
        self.helper_check_quali(self.g4al, [self.f3l])
        self.helper_check_quali(self.g4bl, [self.f4l, self.f6l])
        self.helper_check_quali(self.g4cl, [self.f3l])
        self.helper_check_quali(self.g4dl, [self.f3l])
        self.helper_check_quali(self.g4fl, [self.f14l])

        self.helper_check_quali(self.g2bgl, [self.f2g, self.f2l])
        self.helper_check_quali(self.g3bgl, [self.f3g, self.f3l])
        self.helper_check_quali(self.g3cgl, [self.f3g, self.f3l])
        self.helper_check_quali(self.g4agl, [self.f3g, self.f3l])
        self.helper_check_quali(self.g4bgl, [self.f4g, self.f6g, self.f4l, self.f6l])
        self.helper_check_quali(self.g4cgl, [self.f3g, self.f3l])
        self.helper_check_quali(self.g4dgl, [self.f3g, self.f3l])

    def test_consensus_genefuse_lumpy(self):
        """Test consensus_genefuse_lumpy(self)"""
        # single
        self.helper_build_consensus_single()
        # genefuse_lumpy
        # self.helper_plot_graph(self.g2bgl)
        self.helper_build_consensus_genefuse_lumpy()

        self.helper_check_quanti(self.g1ag, 1, 0, 1, 0, 0)
        self.helper_check_quanti(self.g2ag, 2, 0, 2, 0, 0)
        self.helper_check_quanti(self.g2bg, 3, 2, 2, 1, 0)
        # --
        self.helper_check_quanti(self.g3ag, 3, 0, 3, 0, 0)
        self.helper_check_quanti(self.g3bg, 4, 2, 3, 1, 0)
        self.helper_check_quanti(self.g3cg, 4, 3, 3, 1, 0)
        # --
        self.helper_check_quanti(self.g4ag, 5, 2, 4, 1, 0)
        self.helper_check_quanti(self.g4bg, 6, 4, 4, 2, 0)
        self.helper_check_quanti(self.g4cg, 5, 3, 4, 1, 0)
        # --
        self.helper_check_quanti(self.g4dg, 5, 4, 4, 1, 0)
        self.helper_check_quanti(self.g4eg, 4, 0, 4, 0, 0)

        self.helper_check_quanti(self.g1al, 1, 0, 1, 0, 0)
        self.helper_check_quanti(self.g2al, 2, 0, 2, 0, 0)
        self.helper_check_quanti(self.g2bl, 3, 2, 2, 1, 0)
        # --
        self.helper_check_quanti(self.g3al, 3, 0, 3, 0, 0)
        self.helper_check_quanti(self.g3bl, 4, 2, 3, 1, 0)
        self.helper_check_quanti(self.g3cl, 4, 3, 3, 1, 0)
        # --
        self.helper_check_quanti(self.g4al, 5, 2, 4, 1, 0)
        self.helper_check_quanti(self.g4bl, 6, 4, 4, 2, 0)
        self.helper_check_quanti(self.g4cl, 5, 3, 4, 1, 0)
        # --
        self.helper_check_quanti(self.g4dl, 5, 4, 4, 1, 0)
        self.helper_check_quanti(self.g4el, 4, 0, 4, 0, 0)

        self.helper_check_quanti(self.g1agl, 3, 2, 2, 0, 1)
        # --
        self.helper_check_quanti(self.g2agl, 6, 4, 4, 0, 2)
        self.helper_check_quanti(self.g2bgl, 7, 6, 4, 2, 1)
        # --
        self.helper_check_quanti(self.g3agl, 9, 6, 6, 0, 3)
        self.helper_check_quanti(self.g3bgl, 10, 8, 6, 2, 2)
        self.helper_check_quanti(self.g3cgl, 9, 8, 6, 2, 1)
        # --
        self.helper_check_quanti(self.g4agl, 13, 10, 8, 2, 3)
        self.helper_check_quanti(self.g4bgl, 14, 12, 8, 4, 2)
        self.helper_check_quanti(self.g4cgl, 12, 10, 8, 2, 2)
        self.helper_check_quanti(self.g4dgl, 11, 10, 8, 2, 1)
        self.helper_check_quanti(self.g4egl, 12, 8, 8, 0, 4)

        self.helper_check_quali(self.g2bg, [self.f2g])
        self.helper_check_quali(self.g3bg, [self.f3g])
        self.helper_check_quali(self.g3cg, [self.f3g])
        self.helper_check_quali(self.g4ag, [self.f3g])
        self.helper_check_quali(self.g4bg, [self.f4g, self.f6g])
        self.helper_check_quali(self.g4cg, [self.f3g])
        self.helper_check_quali(self.g4dg, [self.f3g])

        self.helper_check_quali(self.g2bl, [self.f2l])
        self.helper_check_quali(self.g3bl, [self.f3l])
        self.helper_check_quali(self.g3cl, [self.f3l])
        self.helper_check_quali(self.g4al, [self.f3l])
        self.helper_check_quali(self.g4bl, [self.f4l, self.f6l])
        self.helper_check_quali(self.g4cl, [self.f3l])
        self.helper_check_quali(self.g4dl, [self.f3l])

        self.helper_check_quali(self.g1agl, [], [self.f1g])
        self.helper_check_quali(self.g2agl, [], [self.f1g, self.f5g])
        self.helper_check_quali(self.g2bgl, [self.f2l, self.f2g], [self.f2g])
        self.helper_check_quali(self.g3agl, [], [self.f1g, self.f7g, self.f8g])
        self.helper_check_quali(self.g3bgl, [self.f3g, self.f3l], [self.f3g, self.f7g])
        self.helper_check_quali(self.g3cgl, [self.f3g, self.f3l], [self.f3g])
        self.helper_check_quali(
            self.g4agl, [self.f3g, self.f3l], [self.f7g, self.f8g, self.f3g]
        )
        self.helper_check_quali(
            self.g4bgl, [self.f4g, self.f6g, self.f4l, self.f6l], [self.f4g, self.f6g]
        )
        self.helper_check_quali(self.g4cgl, [self.f3g, self.f3l], [self.f3g, self.f6g])
        self.helper_check_quali(self.g4dgl, [self.f3g, self.f3l], [self.f3g])
        self.helper_check_quali(
            self.g4egl, [], [self.f6g, self.f7g, self.f8g, self.f9g]
        )

    def test_define_node_interest(self):
        """Test define_node_interest(self)"""
        self.helper_build_define_node_interest()

        self.helper_check_interest(self.g1ag, 1, 0, 0, [self.f1g])
        self.helper_check_interest(self.g2ag, 2, 0, 0, [self.f1g, self.f5g])
        self.helper_check_interest(self.g2bg, 0, 1, 0, [self.f2g])
        self.helper_check_interest(self.g3ag, 3, 0, 0, [self.f1g, self.f7g, self.f8g])
        self.helper_check_interest(self.g3bg, 0, 1, 0, [self.f3g])
        self.helper_check_interest(self.g3cg, 0, 1, 0, [self.f3g])
        self.helper_check_interest(self.g4ag, 0, 1, 0, [self.f3g])
        self.helper_check_interest(self.g4bg, 0, 2, 0, [self.f4g, self.f6g])
        self.helper_check_interest(self.g4cg, 0, 1, 0, [self.f3g])
        self.helper_check_interest(self.g4dg, 0, 1, 0, [self.f3g])
        self.helper_check_interest(
            self.g4eg, 4, 0, 0, [self.f6g, self.f7g, self.f8g, self.f9g]
        )

        self.helper_check_interest(self.g1al, 1, 0, 0, [self.f1l])
        self.helper_check_interest(self.g2al, 1, 0, 0, [self.f5l])
        self.helper_check_interest(self.g2bl, 0, 1, 0, [self.f2l])
        self.helper_check_interest(self.g3al, 1, 0, 0, [self.f8l])
        self.helper_check_interest(self.g3bl, 0, 1, 0, [self.f3l])
        self.helper_check_interest(self.g3cl, 0, 1, 0, [self.f3l])
        self.helper_check_interest(self.g4al, 0, 1, 0, [self.f3l])
        self.helper_check_interest(self.g4bl, 0, 1, 0, [self.f6l])
        self.helper_check_interest(self.g4cl, 0, 1, 0, [self.f3l])
        self.helper_check_interest(self.g4dl, 0, 1, 0, [self.f3l])
        self.helper_check_interest(self.g4el, 1, 0, 0, [self.f8l])

        self.helper_check_interest(self.g1agl, 0, 0, 1, [self.f1g])
        self.helper_check_interest(self.g2agl, 0, 0, 2, [self.f1g, self.f5g])
        self.helper_check_interest(self.g2bgl, 0, 0, 1, [self.f2g])
        self.helper_check_interest(self.g3agl, 0, 0, 3, [self.f1g, self.f7g, self.f8g])
        self.helper_check_interest(self.g3bgl, 0, 0, 2, [self.f3g, self.f7g])
        self.helper_check_interest(self.g3cgl, 0, 0, 1, [self.f3g])
        self.helper_check_interest(self.g4agl, 0, 0, 3, [self.f7g, self.f8g, self.f3g])
        self.helper_check_interest(self.g4bgl, 0, 0, 2, [self.f4g, self.f6g])
        self.helper_check_interest(self.g4cgl, 0, 0, 2, [self.f3g, self.f6g])
        self.helper_check_interest(self.g4dgl, 0, 0, 1, [self.f3g])
        self.helper_check_interest(
            self.g4egl, 0, 0, 4, [self.f6g, self.f7g, self.f8g, self.f9g]
        )

    def test_trim_node(self):
        """Test trim_node(self)"""
        self.helper_build_define_node_interest()
        self.helper_check_quanti(self.g1at, 14, 8, 10, 0, 4)
        self.helper_check_interest(
            self.g1at, 0, 0, 4, [self.f6g, self.f7g, self.f8g, self.f9g]
        )

        self.g1at.trim_node()

        self.helper_check_quanti(self.g1at, 12, 8, 8, 0, 4)
        self.helper_check_interest(
            self.g1at, 0, 0, 4, [self.f6g, self.f7g, self.f8g, self.f9g]
        )

    def test_label_build_from(self):
        """Test label_build_from(self)"""
        self.helper_build_consensus_single()
        self.helper_build_consensus_genefuse_lumpy()
        # self.f6g (1) self.f2g (2) self.f3g (3) self.f4g (4) self.f6l (5)
        # self.f2l (6) self.f3l (7) self.f4l (8)
        # lev 1 : cons_genefuse (9) cons_lumpy (10)
        # lev 2 : cons (11) cons alone (12)

        self.assertEqual(len(self.g4cgl.label_build_from(1)), 0)
        self.assertEqual(len(self.g4cgl.label_build_from(2)), 0)
        self.assertEqual(len(self.g4cgl.label_build_from(3)), 0)
        self.assertEqual(len(self.g4cgl.label_build_from(4)), 0)
        self.assertEqual(len(self.g4cgl.label_build_from(5)), 0)
        self.assertEqual(len(self.g4cgl.label_build_from(6)), 0)
        self.assertEqual(len(self.g4cgl.label_build_from(7)), 0)
        self.assertEqual(len(self.g4cgl.label_build_from(8)), 0)
        self.assertEqual(self.g4cgl.label_build_from(9), ["GEN_2", "GEN_3", "GEN_4"])
        self.assertEqual(self.g4cgl.label_build_from(10), ["LUM_6", "LUM_7", "LUM_8"])
        self.assertEqual(self.g4cgl.label_build_from(11), ["HMN_10", "HMN_9"])
        self.assertEqual(self.g4cgl.label_build_from(12), ["GEN_1", "LUM_5"])

    def test_get_name_fusion(self):
        """Test get_name(self) from fusion"""
        self.helper_build_consensus_single()
        self.helper_build_consensus_genefuse_lumpy()
        # self.f6g (1) self.f2g (2) self.f3g (3) self.f4g (4) self.f6l (5)
        # self.f2l (6) self.f3l (7) self.f4l (8)
        # lev 1 : cons_genefuse (9) cons_lumpy (10)
        # lev 2 : cons (11) cons alone (12)

        self.assertEqual(self.g4cgl.graph.nodes[1]["fusion"].get_name(), "GEN_1")
        self.assertEqual(self.g4cgl.graph.nodes[2]["fusion"].get_name(), "GEN_2")
        self.assertEqual(self.g4cgl.graph.nodes[3]["fusion"].get_name(), "GEN_3")
        self.assertEqual(self.g4cgl.graph.nodes[4]["fusion"].get_name(), "GEN_4")
        self.assertEqual(self.g4cgl.graph.nodes[5]["fusion"].get_name(), "LUM_5")
        self.assertEqual(self.g4cgl.graph.nodes[6]["fusion"].get_name(), "LUM_6")
        self.assertEqual(self.g4cgl.graph.nodes[7]["fusion"].get_name(), "LUM_7")
        self.assertEqual(self.g4cgl.graph.nodes[8]["fusion"].get_name(), "LUM_8")
        self.assertEqual(self.g4cgl.graph.nodes[9]["fusion"].get_name(), "HMN_9")
        self.assertEqual(self.g4cgl.graph.nodes[10]["fusion"].get_name(), "HMN_10")
        self.assertEqual(self.g4cgl.graph.nodes[11]["fusion"].get_name(), "HMN_11")
        self.assertEqual(self.g4cgl.graph.nodes[12]["fusion"].get_name(), "HMN_12")


if __name__ == "__main__":
    unittest.main()
