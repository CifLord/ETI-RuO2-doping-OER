import itertools, copy

from pymatgen.analysis.structure_matcher import StructureMatcher


def get_dope_slabs(slab, nlayers, dopants):
    
    rounded_fracc = []
    actual_fracc = []
    for site in slab:
        if site.species_string == 'Ru':
            coord = round(site.coords[2], 0)
            if coord not in rounded_fracc:
                rounded_fracc.append(coord)
                actual_fracc.append(site.coords[2])
                
    actual_fracc = sorted(actual_fracc)
    dope_indices = []
    
    for i, site in enumerate(slab):
        for nlayer in nlayers:
            if site.species_string == 'Ru':
                if actual_fracc[-nlayer]-0.5 < site.coords[2] < actual_fracc[-nlayer]+0.5:
                    dope_indices.append(i)
            
    revdopants = copy.deepcopy(dopants)
    doped_slabs = []
    for indices in itertools.combinations(dope_indices, len(dopants)):
        s = slab.copy()
        for n, i in enumerate(indices):
            s.replace(i, dopants[n], properties={'selective_dynamics': [True]*3})
        doped_slabs.append(s)
        if len(dopants) > 1:
            for n, i in enumerate(indices):
                s.replace(i, revdopants[n],  properties={'selective_dynamics': [True]*3})
            doped_slabs.append(s)
        
    sm = StructureMatcher()
    return [g[0] for g in sm.group_structures(doped_slabs)]

def get_selective_dynamics(slab, nlayers=3):
    rounded_fracc = []
    actual_fracc = []
    for site in slab:
        if site.species_string != 'O':
            coord = round(site.coords[2], 0)
            if coord not in rounded_fracc:
                rounded_fracc.append(coord)
                actual_fracc.append(site.coords[2])
                
    actual_fracc = sorted(actual_fracc)
    selective_dynamics = []
    for site in slab:
        if site.coords[2] > actual_fracc[-nlayers]-0.5:
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

