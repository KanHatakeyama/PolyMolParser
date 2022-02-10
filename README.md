# PolyMolParser
- Parse custom polymer-compatible MOL file

# Install
- pip install git+https://github.com/KanHatakeyama/PolyMolParser.git
- RDKit and networkX are required

# Quick start 
- See Tutorial.ipynb for details

1. Prepare molecules
<img src="pics/1.png">

2. Copy as MOL text
<img src="pics/2.png">

3. Run codes
- from PolyMolParser.dict_parse import parse_mol_text
- parse_mol_text(mol_text)
- (`mol_text` is string object of mol data. See https://en.wikipedia.org/wiki/Chemical_table_file)

- Results

```
{0: {0: {'graph': <networkx.classes.graph.Graph at 0x1f5b9862288>,
   'unit_data': {'block_info_n': '4', 'unit_info_n': 1},
   'end_groups': {},
   'unit_MW': 81.13799999999999,
   'SMILES': '*C1CCC(*)C(*)C1',
   'repeated_MW': 81.13799999999999},
  1: {'graph': <networkx.classes.graph.Graph at 0x1f5b98dc088>,
   'unit_data': {'unit_info_n': '10', 'block_info_n': 1},
   'end_groups': {},
   'unit_MW': 46.094,
   'SMILES': '*CS*',
   'repeated_MW': 460.94},
  'unit_MW': 542.078,
  'n': 4.0,
  'repeated_MW': 2168.312},
 'general': {'Mn': 2168.312,
  'Average_MW_per_unit': 197.11927272727272,
  'Average_MW_per_block': 542.078},
 'status': 'successful parsing!'}
```

# Todo
- Refactoring & Bug check

# History
- Prototype version 2022.2.10

# Author
- Kan Hatakeyama-Sato
- Waseda University


