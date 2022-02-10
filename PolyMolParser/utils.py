
UNIT_SYMBOL = "[r]"
BLOCK_SYMBOL = "{r}"
GLOBAL_KEYWORDS = ["Mw", "Mn", "pdi"]
DEFAULT_MW = 10000
DEFAULT_BLOCK_N = 1


def append_key(dict_data, key, val, header=""):
    target_key = key+header
    if target_key not in dict_data:
        dict_data[target_key] = val
    else:
        append_key(dict_data, key, val, header="_")
