import itertools, copy, json

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.surface import *

def get_Ru_type(slab):
    """
    Main used for discerning between Ru-atoms w/ 
    an ontop (Ru_a) or bridge (Ru_b) configuration
    """
    
    site_type = []
    for site in slab:
        if site.species_string == 'O':
            site_type.append('O')
        else:
            Onn = [nn for nn in slab.get_neighbors(site, 2.1) if nn.species_string == 'O']
            if len([nn for nn in Onn if nn.coords[2]-site.coords[2] > 0.1]) == 2:
                site_type.append('Ru_b')
            elif len([nn for nn in Onn if nn.coords[2]-site.coords[2] > 0.1]) == 1:
                site_type.append('Ru_a')
            elif len([nn for nn in Onn if nn.coords[2]-site.coords[2] > 0.1]) == 0:
                if len([nn for nn in Onn if site.coords[2]-nn.coords[2] > 0.1]) == 2:
                    site_type.append('Ru_b')
                elif len([nn for nn in Onn if site.coords[2]-nn.coords[2] > 0.1]) == 1:
                    site_type.append('Ru_a')
                    
    return site_type

def get_nlayers(slab):
    """
    Identify the nth layer of each layer of atoms as a 
    site property. Layer 0 is the topmost layer of atoms
    """
    
    rounded_fracc = []
    actual_fracc = []
    for site in slab:
        if site.species_string == 'Ru':
            coord = round(site.coords[2], 0)
            if coord not in rounded_fracc:
                rounded_fracc.append(coord)
                actual_fracc.append(site.coords[2])
                
    # reverse = True, the topmost layer is 0th
    actual_fracc = sorted(actual_fracc, reverse=True)
    layer_ranges = []
    for i, f in enumerate(actual_fracc):
        if i == 0:
            layer_ranges.append([f-0.1, f+3])
        else:
            layer_ranges.append([f-0.1, actual_fracc[i-1]-0.1])
    nlayers = []
    for i, site in enumerate(slab):
        nth = len(layer_ranges)
        for ii, r in enumerate(layer_ranges):
            if r[0] < site.coords[2] < r[1]:
                nth = ii
        nlayers.append(nth)
    
    return nlayers

def get_superslabs():
    RuO2 = Structure.from_dict(json.load(open('../pmg_structures/RuO2_bulk.json', 'r')))
    vsize = 15
    ssize = 15
    slabgen = SlabGenerator(RuO2, (1,1,0), ssize, vsize, center_slab=True, 
                            lll_reduce=True, max_normal_search=1)
    slabdict = {}
    for i, slab in enumerate(slabgen.get_slabs(symmetrize=True)):
        superslab = slab.copy()
        superslab.make_supercell([2,1,1])
        superslab.make_supercell([[1,1,0], [-1,1,0], [0,0,1]])
        Rutype = get_Ru_type(superslab)
        superslab.add_site_property('Rutype', Rutype)
        nlayers = get_nlayers(superslab)
        superslab.add_site_property('nlayer', nlayers)
        slabdict[i] = superslab
        
    return slabdict

def get_dope_slabs(slab, nlayers, dopants):
    
    dope_indices = []
    for i, site in enumerate(slab):
        for nlayer in nlayers:
            if site.species_string == 'Ru':
                if site.nlayer == nlayer:
                    dope_indices.append(i)
    doped_slabs = []
    for indices in itertools.permutations(dope_indices, len(dopants)):
        s = slab.copy()
        for n, i in enumerate(indices):
            s.replace(i, dopants[n], properties={'selective_dynamics': [True]*3, 
                                                 'Rutype': slab[i].Rutype,
                                                 'nlayer': slab[i].nlayer})
        doped_slabs.append(s)
        
    sm = StructureMatcher()
    return [g[0] for g in sm.group_structures(doped_slabs)]

def get_selective_dynamics(slab, nlayers=3):
    
    selective_dynamics = []
    for site in slab:
        if site.nlayer <= nlayers:
            selective_dynamics.append([True]*3)
        else:
            selective_dynamics.append([False]*3)
    
    return selective_dynamics

def get_adslabs(slab, adscoords, adsconfigs):
    
    all_adslabs = []
    for adscoord in adscoords:
        for ads in adsconfigs:
            adslab = slab.copy()
            adsites = []
            for a in ads:
                adslab.append(a.species_string, np.array(adscoord)+a.coords, 
                              coords_are_cartesian=True, 
                              properties={'selective_dynamics': [True]*3, 
                                          'sitetype': 'adsorbate'})
            all_adslabs.append(adslab)
        
    sm = StructureMatcher()
    if not all_adslabs:
        return []
    return [g[0] for g in sm.group_structures(all_adslabs)]

def get_vac_slabs(slab):
    
    rounded_fracc = []
    actual_fracc = []
    for site in slab:
        if site.species_string == 'O':
            coord = round(site.coords[2], 0)
            if coord not in rounded_fracc:
                rounded_fracc.append(coord)
                actual_fracc.append(site.coords[2])
                
    actual_fracc = sorted(actual_fracc)
    vac_slabs = []
    for i, site in enumerate(slab):
        if site.species_string == 'O':
            if actual_fracc[-1]-0.25 < site.coords[2] < actual_fracc[-1]+0.25:
                vac_slab = slab.copy()
                vacsite = slab[i]
                vac_slab.remove_sites([i])
                setattr(vac_slab, 'vacsite', vacsite) 
                vac_slabs.append(vac_slab)

    sm = StructureMatcher()
    if not vac_slabs:
        return []
    return [g[0] for g in sm.group_structures(vac_slabs)]
