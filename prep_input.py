# Purpose of this program is to prepare the input needed by the tree-diet algorithm.
# Expected input format: dot-bracket.
# Detailed functions:
# a. create edges; b. create bags depending on the edges.

from graph_classes import Graph, Bag
from tree_diet import tree_diet
import random
import networkx as nx
import matplotlib.pyplot as plt
import pydot
from networkx.drawing.nx_pydot import graphviz_layout


# def create_graph(G, dot_bracket_str):
#     # add bases & backbones
#     G.add_node(0)
#     for i in range(1, len(dot_bracket_str)):
#         G.add_node(i)
#         G.add_edge(i - 1, i)
#
#     # add pairings
#     # is it similar to bracket validation?
#     ptr_c = len(dot_bracket_str)  # moves between close brackets (backwards)
#     for i in range(len(dot_bracket_str)):
#         if dot_bracket_str[i] == "(":
#             for j in reversed(range(ptr_c)):
#                 if dot_bracket_str[j] == ")":
#                     G.add_edge(i, j)
#                     ptr_c = j
#                     break
#     return G

def create_graph(G, dot_bracket_str):
    # add bases & backbones
    G.add_node(0)
    stackings = []
    for i in range(1, len(dot_bracket_str)):
        G.add_node(i)
        G.add_edge(i - 1, i)
        stackings.append((i-1, i))


    # add pairings
    # is it similar to bracket validation?
    left_stack = []
    ptr_c = len(dot_bracket_str)  # moves between close brackets (backwards)
    for i in range(len(dot_bracket_str)):
        if dot_bracket_str[i] == "(":
            left_stack.append(i)
        elif dot_bracket_str[i] == ")":
            G.add_edge(left_stack.pop(), i)
        else:
            pass

    return G, stackings


# # Assume two input strings have equal lengths
# def create_graph(G, dot_bracket_str1, dot_bracket_str2):
#     # add bases & backbones
#     G.add_node(0)
#     stackings = []
#
#     for i in range(len(dot_bracket_str1+dot_bracket_str2)):
#
#
#     return G, stackings


def convert_graph(nx_g):
    td_g = Graph()

    # add vertices
    for node in nx_g.nodes:
        td_g.add_vertex(node)

    # add edges
    for edge in nx_g.edges:
        td_g.add_edge(edge[0], edge[1])

    return td_g


# TODO: treewidth_min_degree, treewidth_min_fill_in, OR junction tree?
# Need to dig deeper, hard to tell which is better in our case.

def convert_tree_decomp(tree_decomp):
    rt = Bag([])
    bags_dict = {}

    for node in tree_decomp.nodes:
        # print(node, tree_decomp1.adj[node])
        bags_dict[node] = Bag(list(node))

    solved_item = []
    for bag_frozen in bags_dict:
        if len(tree_decomp.adj[bag_frozen]) > 0:
            for item in tree_decomp.adj[bag_frozen]:
                if item not in solved_item:
                    bags_dict[bag_frozen].add_child(bags_dict[item])
        solved_item.append(bag_frozen)

    # print(bags_dict)

    rt.add_child(bags_dict[list(bags_dict)[0]])

    return rt


def get_adj(g):
    adj = {}

    for bag in g.adj:
        adj[bag] = set()
        if len(g.adj[bag]) > 0:
            for item in g.adj[bag]:
                adj[bag].add(item)

    return adj


if __name__ == '__main__':
    nx_g_1 = nx.Graph()
    nx_g_2 = nx.Graph()

    # nx_g = nx.Graph()

    # G = Graph()

    # ================================================== Input =========================================================

    # Original seqs
    seq1 = "(((..)))"
    seq2 = ".(((.)))"

    nx_g_1, bb1 = create_graph(nx_g_1, seq1)
    nx_g_2, bb2 = create_graph(nx_g_2, seq2)
    print("seq1:", nx_g_1.adj, bb1)
    print("seq2:", nx_g_2.adj, bb2)

    # Coloring
    for u, v, d in nx_g_1.edges(data=True):
        if (u, v) not in bb1:
            d['weight'] = 2
        else:
            d['weight'] = 1

    for u, v, d in nx_g_2.edges(data=True):
        if (u, v) not in bb2:
            d['weight'] = 3
        else:
            d['weight'] = 1

    nx_g = nx.compose(nx_g_1, nx_g_2)
    print("seqC:", nx_g.adj)

    # nx_g, backbone = create_graph(nx_g, "(.(.(.).).(.).)")
    # print("networkx adj:", nx_g.adj)
    # print("backbone:", backbone)

    td_g = convert_graph(nx_g)
    print("tree-diet adj:", td_g.adj, "\n")
    #
    # # print(G, G.nodes, G.edges)
    # # print(G.get_adj())
    #
    # # tree_decomp0 = nx.junction_tree(nx_g)
    # # print("tree0 nodes:", tree_decomp0.nodes)
    # # print("tree-decomp0 adj:", tree_decomp0.adj, "\n")
    #
    # # ==================================================================================================================
    #
    tree_decomp1_width, tree_decomp1 = nx.approximation.treewidth_min_degree(nx_g)
    print("tree1 width:", tree_decomp1_width)
    print("tree-decomp1 nodes:", tree_decomp1.nodes)
    print("tree-decomp1 adjs:", tree_decomp1.adj)
    # tree_decomp1_width, tree_decomp1 = nx.approximation.treewidth_min_fill_in(nx_g)
    # print("tree2 width:", tree_decomp1_width)
    # print("tree-decomp2 nodes:", tree_decomp1.nodes)
    #
    R = convert_tree_decomp(tree_decomp1)

    #
    # nx.draw(nx_g, with_labels=True)
    # nx.draw_networkx(nx_g)
    seed_pos = nx.spring_layout(nx_g, seed=100)
    edges, weights = zip(*nx.get_edge_attributes(nx_g, 'weight').items())
    nx.draw(nx_g, pos=seed_pos, edge_color=weights, width=2, with_labels=True)
    plt.margins(x=0.4)
    plt.show()
    # plt.savefig("/home/songdi/Desktop/coloring.png")
    seed_pos = nx.spring_layout(tree_decomp1, seed=100)
    nx.draw(tree_decomp1, pos=seed_pos, width=2, with_labels=True)
    plt.margins(x=0.4)
    plt.show()
    #
    OPT, real_edges, color_dictionary = tree_diet(R, td_g.adj, 1, bb1)
    print(OPT, real_edges, color_dictionary)
    print(bb1)
