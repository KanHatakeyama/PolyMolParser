import networkx as nx
import matplotlib.pyplot as plt


def draw_chem_graph(g):
    """
    draw chemicals from networkX object

    Parameters
    ----------
    g : networkX object
    """
    plt.figure()
    pos = nx.spring_layout(g)
    #node_labels = nx.get_node_attributes(g,'atomic_num')
    #nx.draw_networkx_labels(g, pos, labels = node_labels,font_size=10, font_color="b")
    node_labels = nx.get_node_attributes(g, 'unit_info')
    nx.draw_networkx_labels(g, pos, labels=node_labels,
                            font_size=16, font_color="r")
    node_labels = nx.get_node_attributes(g, 'block_info')
    nx.draw_networkx_labels(g, pos, labels=node_labels,
                            font_size=16, font_color="g")
    nx.draw_networkx_labels(g, pos, font_size=9, font_color="r")
    nx.draw_networkx_nodes(g, pos, node_size=100, node_color="w")
    nx.draw_networkx_edges(g, pos, width=1)
