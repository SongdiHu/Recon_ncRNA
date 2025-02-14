# Purpose of this package is to prepare the input needed by the tree-diet algorithm.
# Expected input format: dot-bracket.
# Detailed functions:
# a. create edges; b. create bags depending on the edges.
import sys
import networkx as nx
import matplotlib.pyplot as plt

from graph_classes import Graph, Bag
from tree_diet import tree_diet
from MinimumHTree import MHGraph


class StructureMerge:

    def __init__(self, sequences=None):
        self.sequences = sequences
        self.mh_tree_decomps = []
        self.nx_g_bb_list = None
        self.nx_g = None
        self.tree_decomp_width = None
        self.tree_decomp = None

    def add_sequence(self, sequence):
        self.sequences.append(sequence)
        return self.sequences

    def compose_nx_graphs_with_weights(self, nx_g_bb_list=None):
        save = False
        composed_nx_g = None
        if not nx_g_bb_list:
            if not self.nx_g_bb_list:
                self.nx_g_bb_list = self.create_graphs_with_colors()
            nx_g_bb_list = self.nx_g_bb_list
            save = True

        try:
            composed_nx_g = nx_g_bb_list[0][0]
            for nx_g, bb in nx_g_bb_list[1:]:
                composed_nx_g = nx.compose(composed_nx_g, nx_g)
                if save:
                    self.nx_g = composed_nx_g
        except IndexError:
            print("Input list need to be a list of tuples (nx_graph, backbones). "
                  "It needs to contain at least two tuples.")
        return composed_nx_g

    def create_graphs_with_colors(self, list_of_dot_bracket_str=None):
        if not list_of_dot_bracket_str:
            if not self.sequences:
                sys.exit("Not enough sequences. (at least two sequences are required)")
            list_of_dot_bracket_str = self.sequences

        nx_g_bb_list = []
        for seq in list_of_dot_bracket_str:
            nx_g, bb = self.create_graph(seq)
            nx_g_bb_list.append((nx_g, bb))

        # Color the structures
        color_i = 0
        for nx_g, bb in nx_g_bb_list:
            for u, v, d in nx_g.edges(data=True):
                if (u, v) not in bb:
                    d['weight'] = 2 + color_i
                else:
                    d['weight'] = 1
            color_i += 1

        return nx_g_bb_list

    def create_graph(self, dot_bracket_str):
        # add bases & backbones
        nx_graph = nx.Graph()
        nx_graph.add_node(0)
        stackings = []
        for i in range(1, len(dot_bracket_str)):
            nx_graph.add_node(i)
            nx_graph.add_edge(i - 1, i)
            stackings.append((i - 1, i))

        # add pairings
        # is it similar to bracket validation?
        left_stack = []
        ptr_c = len(dot_bracket_str)  # moves between close brackets (backwards)
        for i in range(len(dot_bracket_str)):
            if dot_bracket_str[i] == "(":
                left_stack.append(i)
            elif dot_bracket_str[i] == ")":
                nx_graph.add_edge(left_stack.pop(), i)
            else:
                pass

        return nx_graph, stackings

    def root_undirected_tree(self, tree_decomp=None):
        if tree_decomp is None:
            if not self.tree_decomp:
                self.tree_decomp = self.get_tree_decomp()[1]
            tree_decomp = self.tree_decomp

        list_of_rooted_trees = []
        # print("depth:", nx.shortest_path_length(self.tree_decomp, list(self.tree_decomp.nodes)[0]))
        mhg = MHGraph(tree_decomp.nodes)

        for edge in tree_decomp.edges:
            mhg.addEdge(edge[0], edge[1])

        root_cands = mhg.rootForMinimumHeight()
        for cand in root_cands:
            list_of_rooted_trees.append(nx.bfs_tree(tree_decomp, mhg.get_bag(cand)))
        if tree_decomp is None:
            self.mh_tree_decomps = list_of_rooted_trees

        # print("MHT nodes:", mh_tree_decomps[0].nodes())
        # print("MHT edges:", mh_tree_decomps[0].edges())
        # print("depth:", nx.shortest_path_length(mh_tree_decomps[0], list(self.mh_tree_decomps[0].edges())[0][0]))
        
        return list_of_rooted_trees

    def get_tree_diet(self, tw, balanced_or_not=None):
        adj = self.to_td_graph().adj
        if balanced_or_not:
            if not self.mh_tree_decomps:
                self.mh_tree_decomps = self.root_undirected_tree()[0]
            root_bag, b_tags = self.to_td_bags(self.mh_tree_decomps)
        else:
            if not self.tree_decomp:
                self.tree_decomp = self.get_tree_decomp()[1]
            root_bag, b_tags = self.to_td_bags(self.tree_decomp)
        # print("backbone?", self.nx_g_bb_list[0][1])
        return tree_diet(root_bag, adj, tw, [], b_tags)


    # ====================================== CONVERSIONS ======================================
    def to_td_graph(self, nx_g=None):
        if not nx_g:
            if not self.nx_g:
                self.nx_g = self.compose_nx_graphs_with_weights()
            nx_g = self.nx_g
        self.td_g = Graph()

        # add vertices
        for i in nx_g.nodes:
            self.td_g.add_vertex(int(i))

        # add edges
        for edge in nx_g.edges:
            self.td_g.add_edge(edge[0], edge[1])

        # print("td graph:", self.td_g.adj)
        return self.td_g

    def to_td_bags(self, tree):
        nodes_list_original = list(tree.nodes)
        # print("bag order", nodes_list_original)
        edges_list = list(tree.edges)
        bags_list = [Bag(list(bag)) for bag in nodes_list_original]

        R = Bag([])
        R.add_child(bags_list[0])
        for edge in edges_list:
            bags_list[nodes_list_original.index(edge[0])].add_child(bags_list[nodes_list_original.index(edge[1])])

        b_tags = dict()
        b_tags.update({R: -1})
        for i in range(len(bags_list)):
            b_tags.update({bags_list[i]: i})

        return R, b_tags

# ====================================== GETTERS ======================================
    def get_adj(self, g):
        adj = {}

        for bag in g.adj:
            adj[bag] = set()
            if len(g.adj[bag]) > 0:
                for item in g.adj[bag]:
                    adj[bag].add(item)

        return adj

    # test run on "min-degree"
    def get_tree_decomp(self, nx_g=None, method="min-degree"):
        tree_width, tree = None, None
        td_self = False

        if nx_g is None:
            td_self = True
            if not self.nx_g:
                self.nx_g = self.compose_nx_graphs_with_weights()
            nx_g = self.nx_g

        try:
            match method:
                case "min-degree":
                    tree_width, tree = nx.approximation.treewidth_min_degree(nx_g)
                case "min-fill-in":
                    tree_width, tree = nx.approximation.treewidth_min_fill_in(nx_g)
        except TypeError:
            print("Please input a valid nx graph.")

        if td_self:
            self.tree_decomp = tree
            self.tree_decomp_width = tree_width

        return tree_width, tree

    def get_mh_tree(self):
        return self.mh_tree_decomps[0], self.tree_decomp_width

    def get_bags_in_order(self, tree):
        bags_in_order = list(tree.nodes)
        bags_in_order = [tuple(positions) for positions in bags_in_order]
        return bags_in_order

# ====================================== RESULTS ======================================
    def print_and_plot_outputs(self):
        if self.nx_g:
            print("\nNetworkX graph vertices:", self.nx_g.nodes)
            print("NetworkX graph edges:", self.nx_g.edges)

            seed_pos = nx.spring_layout(self.nx_g, seed=100)
            plt.title("Network (original edges)")
            edges, weights = zip(*nx.get_edge_attributes(self.nx_g, 'weight').items())
            nx.draw(self.nx_g, pos=seed_pos, edge_color=weights, width=2, with_labels=True)
            # plt.margins(x=0.4)
            plt.show()
        # plt.savefig("/home/songdi/Desktop/coloring.png")

        if self.tree_decomp:
            print("\ntree width:", self.tree_decomp_width)
            print("# of bags:", len(self.tree_decomp.nodes))
            print("tree-decomp bags:", self.tree_decomp.nodes)
            print("tree-decomp adjs:", self.tree_decomp.adj)

            seed_pos = nx.spring_layout(self.tree_decomp, seed=100)
            plt.title("Tree decomposition (original bags)")
            nx.draw(self.tree_decomp, pos=seed_pos, width=2, with_labels=True)
            # plt.margins(x=0.4)
            plt.show()

        if self.mh_tree_decomps:
            print("\nMinimum-height tree width:", self.mh_tree_decomps[0])
            print("# of bags:", len(self.mh_tree_decomps[0].nodes))
            print("Minimum-height tree bags:", self.mh_tree_decomps[0].nodes)
            print("Minimum-height tree adjs:", self.mh_tree_decomps[0].adj)

            seed_pos = nx.spring_layout(self.mh_tree_decomps[0], seed=100)
            plt.title("Balanced & Rooted Tree decomposition (original bags)")
            nx.draw(self.mh_tree_decomps[0], pos=seed_pos, width=2, with_labels=True)
            # plt.margins(x=0.4)
            plt.show()
