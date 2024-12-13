import json, glob

from pymatgen.core.surface import *


def load_structures():
    
    all_ads_dict = {'OOH': [], 'OH': [], 'O': []}
    for f in glob.glob('../pmg_structures/adsorbates/*'):
        if '.xyz' not in f:
            continue
        n = f.split('/')[-1].replace('.xyz', '')
        if 'OOH' in n:
            all_ads_dict['OOH'].append(Molecule.from_file(f))
        elif 'OH' in n:
            all_ads_dict['OH'].append(Molecule.from_file(f))
        else:
            all_ads_dict['O'].append(Molecule.from_file(f))
            
    slabdict = json.load(open('../pmg_structures/RuO2_slabs.json', 'r'))
    slabdict = {i: Slab.from_dict(slabdict[i]) for i in slabdict.keys()}
    bulk = Structure.from_dict(json.load(open('../pmg_structures/RuO2_bulk.json', 'r')))
    
    return all_ads_dict, slabdict, bulk


def doped_RuO2_OER_workflow(dopants_list, nlayers_list):

    all_ads_dict, slabdict, bulk = load_structures()
    return None
