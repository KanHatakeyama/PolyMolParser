
from . utils import *
import re
import networkx as nx
from rdkit import Chem

# https://github.com/maxhodak/keras-molecules/pull/32/files


def mol_to_nx(mol):
    G = nx.Graph()

    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(),
                   atomic_num=atom.GetAtomicNum(),
                   formal_charge=atom.GetFormalCharge(),
                   chiral_tag=atom.GetChiralTag(),
                   hybridization=atom.GetHybridization(),
                   num_explicit_hs=atom.GetNumExplicitHs(),
                   is_aromatic=atom.GetIsAromatic())

    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType())
    return G


# https://github.com/maxhodak/keras-molecules/pull/32/files
def nx_to_mol(G):
    mol = Chem.RWMol()
    atomic_nums = nx.get_node_attributes(G, 'atomic_num')
    chiral_tags = nx.get_node_attributes(G, 'chiral_tag')
    formal_charges = nx.get_node_attributes(G, 'formal_charge')
    node_is_aromatics = nx.get_node_attributes(G, 'is_aromatic')
    node_hybridizations = nx.get_node_attributes(G, 'hybridization')
    num_explicit_hss = nx.get_node_attributes(G, 'num_explicit_hs')
    node_to_idx = {}
    for node in G.nodes():
        a = Chem.Atom(atomic_nums[node])
        a.SetChiralTag(chiral_tags[node])
        a.SetFormalCharge(formal_charges[node])
        a.SetIsAromatic(node_is_aromatics[node])
        a.SetHybridization(node_hybridizations[node])
        a.SetNumExplicitHs(num_explicit_hss[node])
        idx = mol.AddAtom(a)
        node_to_idx[node] = idx

    bond_types = nx.get_edge_attributes(G, 'bond_type')
    for edge in G.edges():
        first, second = edge
        ifirst = node_to_idx[first]
        isecond = node_to_idx[second]
        bond_type = bond_types[first, second]
        mol.AddBond(ifirst, isecond, bond_type)

    Chem.SanitizeMol(mol)
    return mol


def get_special_atom_dict(mol_text):

    lines = mol_text.split("\n")

    summary_line = lines[3]
    mol_info = summary_line.split(" ")
    mol_info = [i for i in mol_info if i != ""]
    n_atoms = mol_info[0]
    n_bonds = mol_info[1]
    n_bonds = int(re.findall(r'\d+', n_bonds)[0])
    n_atoms = int(re.findall(r'\d+', n_atoms)[0])

    # get property lines
    prop_lines = lines[4+n_atoms+n_bonds:-1]

    key_lines = prop_lines[::2]
    val_lines = prop_lines[1::2]

    special_atom_dict = {}
    for k, v in zip(key_lines, val_lines):
        key = int(re.findall(r'\d+', k)[0])
        special_atom_dict[key] = v

    return special_atom_dict


def parse_unit_keys(original_str_keywords, unit_symbol=UNIT_SYMBOL):
    keyword_dict = {}
    if original_str_keywords.find(unit_symbol) == -1:
        return None

    str_keys = original_str_keywords.replace(unit_symbol, "").replace(" ", "")
    key_list = str_keys.split(",")

    for query in key_list:
        if query.find("=") == -1:
            continue
        k, v = query.split("=")
        keyword_dict[k] = v

    return keyword_dict


def apply_unit_info(g, text):
    special_atom_dict = get_special_atom_dict(text)

    for k in special_atom_dict:
        parsed_dict = parse_unit_keys(
            special_atom_dict[k], unit_symbol=BLOCK_SYMBOL)
        if parsed_dict is not None:
            g.nodes[k-1]["block_info"] = parsed_dict

    for k in special_atom_dict:
        parsed_dict = parse_unit_keys(special_atom_dict[k])
        if parsed_dict is not None:
            g.nodes[k-1]["unit_info"] = parsed_dict
