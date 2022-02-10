from rdkit.Chem.Descriptors import MolWt
from .mol_parse import nx_to_mol
from rdkit import Chem
from .utils import *
import copy
from .graph_parse import prepare_dict


def parse_mol_text(text):
    fragment_dict = {}
    try:
        fragment_dict = prepare_dict(text)
    except:
        fragment_dict["status"] = "Error initial parse!"
        return fragment_dict

    try:
        attribute_mw_and_smiles(fragment_dict)
    except:
        fragment_dict["status"] = "Error parsing smiles."
        return fragment_dict

    try:
        calc_MW_from_n(fragment_dict)
    except:
        fragment_dict["status"] = "Error Calculating molecular weight"
        return fragment_dict

    try:
        get_global_keywords(fragment_dict)
    except:
        fragment_dict["status"] = "Error parsing Mw, Mn, pdi"
        return fragment_dict

    try:
        calc_MW_from_x(fragment_dict)
    except:
        fragment_dict["status"] = "Error parsing x"
        return fragment_dict

    try:
        calc_average_weight(fragment_dict)
    except:
        fragment_dict["status"] = "error calculating average MW"
        return fragment_dict

    fragment_dict["status"] = "successful parsing!"
    return fragment_dict


def calc_average_weight(fragment_dict):

    block_n_list = []
    unit_n_list = []
    for i in fragment_dict:
        if i == "general":
            continue

        block = fragment_dict[i]

        block_n = -1

        if "n" in block:
            block_n = float(block["n"])
            block_n_list.append(block_n)

        for j in block:
            unit = block[j]
            if type(unit) is dict:
                unit_n_list.append(float(unit["unit_data"]["unit_info_n"]))

                if block_n == -1:
                    block_n = float(unit["unit_data"]["block_info_n"])
                    block_n_list.append(block_n)

    fragment_dict["general"]["Average_MW_per_unit"] = fragment_dict["general"]["Mn"] / \
        sum(unit_n_list)
    fragment_dict["general"]["Average_MW_per_block"] = fragment_dict["general"]["Mn"] / \
        sum(block_n_list)


def attribute_mw_and_smiles(fragment_dict):
    for i in fragment_dict:
        for j in fragment_dict[i]:
            current_chem_data = fragment_dict[i][j]
            sm, mw = get_sm_and_mw(current_chem_data["graph"])
            current_chem_data["unit_MW"] = mw
            current_chem_data["SMILES"] = sm

            for end_chem_id in current_chem_data["end_groups"]:
                sm, mw = get_sm_and_mw(
                    current_chem_data["end_groups"][end_chem_id]["graph"])
                current_chem_data["end_groups"][end_chem_id]["unit_MW"] = mw
                current_chem_data["end_groups"][end_chem_id]["SMILES"] = sm


def get_sm_and_mw(g):

    mol = nx_to_mol(g)
    sm = Chem.MolToSmiles(mol)
    sm = sm.replace("**", "*")
    mw = MolWt(mol)

    return sm, mw


def calc_repeated_MW(chem_data, key="unit_info_n"):
    # calc MW of polymer repeated unit
    if not key in chem_data["unit_data"]:
        if key == "unit_info_n":
            chem_data["unit_data"][key] = 1

    if key in chem_data["unit_data"]:
        chem_data["repeated_MW"] = float(
            chem_data["unit_data"][key])*chem_data["unit_MW"]

        # end group
        for end_id in chem_data["end_groups"]:
            chem_data["repeated_MW"] += chem_data["end_groups"][end_id]["unit_MW"]


def calc_MW_from_n(fragment_dict, mode="n"):
    # calculate MW from n
    total_MW = 0
    for i in fragment_dict:
        block_MW = 0
        block_n = 0
        block_list = []
        for j in fragment_dict[i]:
            if type(j) is not int:
                continue

            chem_data = fragment_dict[i][j]
            calc_repeated_MW(chem_data, key=f"unit_info_{mode}")
            if "repeated_MW" in chem_data:
                block_MW += chem_data["repeated_MW"]

            if f"block_info_{mode}" not in chem_data["unit_data"]:
                chem_data["unit_data"][f"block_info_{mode}"] = DEFAULT_BLOCK_N

            block_n = chem_data["unit_data"][f"block_info_{mode}"]
            block_list.append(block_n)

        if block_MW != 0:
            block_data = fragment_dict[i]
            block_data["unit_MW"] = block_MW
            if block_n != 0:
                # block_data[f"{mode}"]=float(block_n)
                block_data[f"{mode}"] = max([float(i) for i in block_list])
                block_data["repeated_MW"] = block_data[f"{mode}"]*block_MW
                total_MW += block_data["repeated_MW"]

    if total_MW != 0:
        fragment_dict["general"]["Mn"] = total_MW


def get_global_keywords(fragment_dict):
    for i in fragment_dict:
        if i == "general":
            continue
        for j in fragment_dict[i]:
            if type(j) is not int:
                continue
            if "unit_data" in fragment_dict[i][j]:
                for search_key in GLOBAL_KEYWORDS:
                    for current_key in fragment_dict[i][j]["unit_data"]:
                        if current_key.find(search_key) != -1:
                            fragment_dict["general"][search_key] = float(
                                fragment_dict[i][j]["unit_data"][current_key])
                            break


# in the case of x
def calc_MW_from_x(fragment_dict):
    temp_graph_dict = copy.deepcopy(fragment_dict)

    # tempolarily estimate Mn from x
    calc_MW_from_n(temp_graph_dict, mode="x")
    estimated_Mn = temp_graph_dict["general"]["Mn"]

    if "Mn" in fragment_dict["general"]:
        target_mw = fragment_dict["general"]["Mn"]
    elif "Mw" in fragment_dict["general"]:
        target_mw = fragment_dict["general"]["Mw"]
    else:
        target_mw = DEFAULT_MW

    ratio = target_mw/estimated_Mn

    # update n
    for i in fragment_dict:
        if i == "general":
            continue
        for j in fragment_dict[i]:
            if type(j) is not int:
                continue
            if "unit_data" in list(fragment_dict[i][j]):
                for k in list(fragment_dict[i][j]["unit_data"]):
                    if k in ["unit_info_x", "block_info_x"]:
                        fragment_dict[i][j]["unit_data"][k.replace("_x", "_n")] = float(
                            fragment_dict[i][j]["unit_data"][k])*ratio
