from . utils import *
import networkx as nx
import copy
from rdkit import Chem
from .mol_parse import mol_to_nx, apply_unit_info


def prepare_dict(text):

    mol = Chem.MolFromMolBlock(text)

    # prepare netowork x object
    g = mol_to_nx(mol)

    # set unit info
    apply_unit_info(g, text)

    # share connection types of neiboghring [r] nodes
    share_connect_type(g, unit_key="block_info")
    share_connect_type(g)

    # share Mw, Mn, ... as global variables
    # share_global_keywords(g,unit_key="block_info")
    # share_global_keywords(g)

    # fragmentate by block
    block_dict = fragmentate_graph(g, dict_key="block_info")

    # fragmentate by polymer unit
    fragment_dict = {}
    for i, block_id in enumerate(block_dict):
        fragment_dict[i] = fragmentate_graph(
            block_dict[block_id]["graph"], dict_key="unit_info")

    # save unit info to dict
    update_unit_info(fragment_dict)

    # separaete end groups
    separate_end_groups(fragment_dict)

    replace_unit_nodes(fragment_dict)

    fragment_dict["general"] = {}
    return fragment_dict


def has_specific_neigbor(g, node_id, unit_key):
    for neighbor_id in nx.neighbors(g, node_id):
        if unit_key in g.nodes[neighbor_id]:
            return neighbor_id
    return -1


# share connect information with neighboring [r] nodes
def share_connect_type(g, unit_key="unit_info"):
    for node_id in g.nodes:
        if unit_key in g.nodes[node_id]:
            if "c" in g.nodes[node_id][unit_key]:
                connect_type = g.nodes[node_id][unit_key]["c"]

                neighbor_id = has_specific_neigbor(g, node_id, unit_key)
                if neighbor_id != -1:
                    g.nodes[neighbor_id][unit_key]["c"] = connect_type


# share global keywords
def share_global_keywords(g, unit_key="unit_info", header="mutual_", keywords=GLOBAL_KEYWORDS):
    for keyword in keywords:
        for node_id in g.nodes:
            if unit_key in g.nodes[node_id]:
                if keyword in list(g.nodes[node_id][unit_key]):
                    for apply_node_id in g.nodes:
                        if unit_key in g.nodes[apply_node_id]:
                            g.nodes[apply_node_id][unit_key][f"{header}{keyword}"] = g.nodes[node_id][unit_key][keyword]

                    g.nodes[node_id][unit_key].pop(keyword)


def fragmentate_graph(g, dict_key="block_info"):
    g = copy.deepcopy(g)

    # search for neighboring blocks
    for node_id in g.nodes:
        if dict_key in g.nodes[node_id]:
            neighbor_id = has_specific_neigbor(g, node_id, dict_key)
            if neighbor_id != -1:
                g.remove_edge(node_id, neighbor_id)

    block_dict = {}
    for i, c in enumerate(nx.connected_components(g)):
        mini_dict = {}
        mini_dict["graph"] = nx.Graph(g.subgraph(c))
        block_dict[i] = mini_dict

    return block_dict


# collect local keywords
def get_local_polymer_dict(fragment_g):
    local_polymer_dict = {}
    for node_id in fragment_g.nodes:
        for info_key in ["unit_info", "block_info"]:
            if info_key in fragment_g.nodes[node_id]:
                for k, v in fragment_g.nodes[node_id][info_key].items():
                    append_key(local_polymer_dict, f"{info_key}_{k}", v)

    return local_polymer_dict


def update_unit_info(fragment_dict):
    for block_id in fragment_dict:
        for fragment_id in fragment_dict[block_id]:
            component_dict = fragment_dict[block_id][fragment_id]
            fragment_g = component_dict["graph"]
            component_dict["unit_data"] = get_local_polymer_dict(fragment_g)


# search for endgroups of polymer units


def get_repeating_unit_numbers(fragment_g, current_node, id_mode=False):
    global explored_nodes
    explored_nodes = set()
    global unit_count
    unit_count = 0

    def inner_loop(fragment_g, current_node, id_mode):
        global unit_count
        global explored_nodes
        explored_nodes.add(current_node)

        next_cand_nodes = set(fragment_g.neighbors(current_node))
        next_cand_nodes = next_cand_nodes-explored_nodes

        for next_node in next_cand_nodes:
            # stop searching in the case of unit node
            unit_flag = False
            for info_key in ["unit_info", "block_info"]:
                if info_key in fragment_g.nodes[next_node]:
                    unit_flag = True
                    break

            if id_mode:
                if unit_flag:
                    return current_node, next_node, explored_nodes
                else:
                    return inner_loop(fragment_g, next_node, id_mode)

            if unit_flag:
                unit_count += 1
            else:
                inner_loop(fragment_g, next_node, id_mode)

    if id_mode:
        return inner_loop(fragment_g, current_node, id_mode)
    else:
        inner_loop(fragment_g, current_node, id_mode)
        return unit_count


def obtain_end_node_list(fragment_g):
    # select end candidates
    end_cand_nodes = []
    for node_id in fragment_g.nodes:
        if len(list(fragment_g.neighbors(node_id))) == 1:
            neighbor_flag = True
            for info_key in ["unit_info", "block_info"]:
                if info_key in fragment_g.nodes[node_id]:
                    neighbor_flag = False

            if neighbor_flag:
                end_cand_nodes.append(node_id)

    end_nodes = []
    for cand in end_cand_nodes:

        # end group must reach only one unit node
        num_units = get_repeating_unit_numbers(fragment_g, cand)
        if num_units == 1:
            end_nodes.append(cand)

    return end_nodes


# separate end groups
def separate_end_groups(fragment_dict):
    for i in fragment_dict:
        for j in fragment_dict[i]:
            fragment_g = fragment_dict[i][j]["graph"]

            end_node_list = obtain_end_node_list(fragment_g)
            end_graph_dict = {}

            for num, end_node in enumerate(end_node_list):
                temp_dict = {}
                # search for cutting part and endgroup nodes
                edge1, edge2, end_ids = get_repeating_unit_numbers(
                    fragment_g, end_node, id_mode=True)

                fragment_g.remove_edge(edge1, edge2)
                end_graph = nx.Graph(fragment_g.subgraph(end_ids))
                temp_dict["graph"] = end_graph
                fragment_dict[i][j]["graph"] = nx.Graph(
                    fragment_g.subgraph(set(fragment_g.nodes)-end_ids))
                fragment_g = nx.Graph(fragment_g.subgraph(
                    set(fragment_g.nodes)-end_ids))

                end_graph_dict[num] = temp_dict

            fragment_dict[i][j]["end_groups"] = end_graph_dict

# change unit node info


def replace_unit_nodes(fragment_dict):
    for i in fragment_dict:
        for j in fragment_dict[i]:
            cut_g = fragment_dict[i][j]["graph"]

            for node_id in list(cut_g.nodes):
                for info_key in ["unit_info", "block_info"]:
                    if info_key in cut_g.nodes[node_id]:
                        # cut_g.remove_node(node_id)
                        cut_g.nodes[node_id]["atomic_num"] = 0
                        # break
