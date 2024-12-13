import json, glob

from pymatgen.core.surface import *


def get_adsorption_sites():
    
    slabdict = json.load(open('../pmg_structures/RuO2_slabs.json', 'r'))
    
    adspos = {0: []}

    slab = slabdict[1].copy()
    Opos = []
    for site in slab:
        if site.species_string == 'O':
            pos = float('%.2f' %(site.coords[2]))
            if pos not in Opos:
                Opos.append(pos)
    cpos = sorted(Opos)[-2]+(sorted(Opos)[-6] - sorted(Opos)[-8])
    ocoords = []
    for site in slab:
        if site.species_string == 'O':
            if sorted(Opos)[-6] == float('%.2f' %(site.coords[2])):
                ocoords.append([site.coords[0], site.coords[1], cpos])
    adspos[1] = ocoords

    slab = slabdict[2].copy()
    Opos = []
    for site in slab:
        if site.species_string == 'O':
            pos = float('%.2f' %(site.coords[2]))
            if pos not in Opos:
                Opos.append(pos)
    cpos = sorted(Opos)[-2]+(sorted(Opos)[-6] - sorted(Opos)[-8])
    ocoords = []
    for site in slab:
        if site.species_string == 'O':
            if sorted(Opos)[-6] == float('%.2f' %(site.coords[2])):
                ocoords.append([site.coords[0], site.coords[1], cpos])
    adspos[2] = ocoords

    return adspos

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
    adspos = json.load(open('adsorption_sites.json', 'r'))    

    return all_ads_dict, slabdict, bulk, adspos


def doped_RuO2_OER_workflow(dopants_list, nlayers_list):

    all_ads_dict, slabdict, bulk, adspos = load_structures()
    return None


if __name__ == "__main__":

    maind = '/project/grabow/rtran25/ETI_RuO2_doping/'

    for term in tqdm(slabdict.keys()):
        superslab = slabdict[term].copy()

        superslab.add_site_property('selective_dynamics', get_selective_dynamics(superslab))
        a = AseAtomsAdaptor.get_atoms(superslab)
        a.set_tags([1 if all(sel) else 0 for sel in superslab.site_properties['selective_dynamics']])
        magmoms = [0.6]*len(a)
        n = 'RuO2_110_term%s' %(term)
        new_f = os.path.join('redo_all_again/%s' %(n))
        os.mkdir(new_f)
        superslab.to(os.path.join(new_f, '%s.cif' %(n)))
        shutil.copyfile('submit_script', os.path.join(new_f, 'submit_script'))
        metadata = {'structure': superslab.as_dict()}
        json.dump(metadata, open(os.path.join(new_f, 'metadata.json'), 'w'))
        os.chdir(new_f)
        write_inputs(a, mags=magmoms)
        os.chdir(maind)

        for adsname in all_ads_dict.keys():
            all_adslabs = get_adslabs(superslab, adspos[term], all_ads_dict[adsname])
            for adspos_i, adslab in enumerate(all_adslabs):

                a = AseAtomsAdaptor.get_atoms(adslab)
                a.set_tags([1 if all(sel) else 0 for sel in adslab.site_properties['selective_dynamics']])
                magmoms = [0.6]*len(a)
                n = 'RuO2_110_term%s_%s%s' %(term, adsname, adspos_i)
                new_f = os.path.join('redo_all_again/%s' %(n))
                os.mkdir(new_f)
                adslab.to(os.path.join(new_f, '%s.cif' %(n)))
                shutil.copyfile('submit_script', os.path.join(new_f, 'submit_script'))
                metadata = {'structure': adslab.as_dict(), 'adsites': [site.as_dict() for site in adslab \
                                                                       if 'sitetype' in site.properties]}
                json.dump(metadata, open(os.path.join(new_f, 'metadata.json'), 'w'))
                os.chdir(new_f)
                write_inputs(a, mags=magmoms)
                os.chdir(maind)

        vac_slabs = get_vac_slabs(superslab)
        for vac_i, vac_slab in enumerate(vac_slabs):
            a = AseAtomsAdaptor.get_atoms(vac_slab)
            a.set_tags([1 if all(sel) else 0 for sel in vac_slab.site_properties['selective_dynamics']])
            magmoms = [0.6]*len(a)
            n = 'RuO2_110_term%s_vac%s' %(term, vac_i)
            new_f = os.path.join('redo_all_again/%s' %(n))
            os.mkdir(new_f)
            vac_slab.to(os.path.join(new_f, '%s.cif' %(n)))
            shutil.copyfile('submit_script', os.path.join(new_f, 'submit_script'))
            metadata = {'structure': vac_slab.as_dict(), 'vac_site': vac_slab.vacsite.as_dict()}
            json.dump(metadata, open(os.path.join(new_f, 'metadata.json'), 'w'))
            os.chdir(new_f)
            write_inputs(a, mags=magmoms)
            os.chdir(maind)

            for adsname in all_ads_dict.keys():
                if adsname == 'O':
                    continue

                all_adslabs = get_adslabs(vac_slab, [vac_slab.vacsite.coords], all_ads_dict[adsname])
                for adspos_i, adslab in enumerate(all_adslabs):

                    a = AseAtomsAdaptor.get_atoms(adslab)
                    a.set_tags([1 if all(sel) else 0 for sel in adslab.site_properties['selective_dynamics']])
                    magmoms = [0.6]*len(a)
                    n = 'RuO2_110_term%s_vac%s_%s%s' %(term, vac_i, adsname, adspos_i)
                    new_f = os.path.join('redo_all_again/%s' %(n))
                    os.mkdir(new_f)
                    adslab.to(os.path.join(new_f, '%s.cif' %(n)))
                    shutil.copyfile('submit_script', os.path.join(new_f, 'submit_script'))
                    metadata = {'structure': adslab.as_dict(), 
                                'adsites': [site.as_dict() for site in adslab \
                                            if 'sitetype' in site.properties],
                                'vac_site': vac_slab.vacsite.as_dict()}
                    json.dump(metadata, open(os.path.join(new_f, 'metadata.json'), 'w'))
                    os.chdir(new_f)
                    write_inputs(a, mags=magmoms)
                    os.chdir(maind)

    #     for dopants in [['Co'], ['Mn'], ['Zn'], ['Ti'], ['Nb'], 
    #                     ['Co', 'Mn'], ['Co', 'Nb'], ['Mn', 'Nb']]:
        for dopants in [['Mn', 'Nb']]:

            doped_slabs = get_dope_slabs(superslab, [1,2], dopants)

            for dopepos, superslab_dope in enumerate(doped_slabs):

                Mcomp = 'Ru%s' %(int(superslab_dope.composition['Ru']))
                for d in dopants:
                    Mcomp+='%s%s' %(d, int(superslab_dope.composition[d]))

                domags = ['nm', 'fm'] if any([el in dopants for el in ['Co', 'Mn']]) else ['nm']

                for mag in domags:

                    a = AseAtomsAdaptor.get_atoms(superslab_dope)
                    a.set_tags([1 if all(sel) else 0 for sel \
                                in superslab_dope.site_properties['selective_dynamics']])
                    magmoms = [0.6]*len(a) if mag == 'nm' else \
                    [0.6 if at.symbol not in ['Co', 'Mn'] else 1 for at in a]
                    n = 'RuO2_%s_dope%s_110_term%s_%s' %(Mcomp, dopepos, term, mag)
                    new_f = os.path.join('redo_all_again/%s' %(n))
                    os.mkdir(new_f)
                    superslab_dope.to(os.path.join(new_f, '%s.cif' %(n)))
                    shutil.copyfile('submit_script', os.path.join(new_f, 'submit_script'))
                    metadata = {'structure': superslab_dope.as_dict(),
                                'dopant_sites': [site.as_dict() for site in superslab_dope \
                                                 if site.species_string not in ['O', 'Ru', 'H']]}
                    json.dump(metadata, open(os.path.join(new_f, 'metadata.json'), 'w'))
                    os.chdir(new_f)
                    write_inputs(a, mags=magmoms)
                    os.chdir(maind)

                    vac_slabs = get_vac_slabs(superslab_dope)
                    for vac_i, vac_slab in enumerate(vac_slabs):
                        a = AseAtomsAdaptor.get_atoms(vac_slab)
                        a.set_tags([1 if all(sel) else 0 for sel \
                                    in vac_slab.site_properties['selective_dynamics']])
                        magmoms = [0.6]*len(a) if mag == 'nm' else \
                        [0.6 if at.symbol not in ['Co', 'Mn'] else 1 for at in a]
                        n = 'RuO2_%s_dope%s_110_term%s_%s_vac%s' %(Mcomp, dopepos, term, mag, vac_i)
                        new_f = os.path.join('redo_all_again/%s' %(n))
                        os.mkdir(new_f)
                        vac_slab.to(os.path.join(new_f, '%s.cif' %(n)))
                        shutil.copyfile('submit_script', os.path.join(new_f, 'submit_script'))
                        metadata = {'structure': vac_slab.as_dict(), 'vac_site': vac_slab.vacsite.as_dict(),
                                    'dopant_sites': [site.as_dict() for site in superslab_dope \
                                                     if site.species_string not in ['O', 'Ru', 'H']]}
                        json.dump(metadata, open(os.path.join(new_f, 'metadata.json'), 'w'))
                        os.chdir(new_f)
                        write_inputs(a, mags=magmoms)
                        os.chdir(maind)

                        for adsname in all_ads_dict.keys():
                            if adsname == 'O':
                                continue

                            all_adslabs = get_adslabs(vac_slab, [vac_slab.vacsite.coords], all_ads_dict[adsname])
                            for adspos_i, adslab in enumerate(all_adslabs):

                                a = AseAtomsAdaptor.get_atoms(adslab)
                                a.set_tags([1 if all(sel) else 0 \
                                            for sel in adslab.site_properties['selective_dynamics']])
                                magmoms = [0.6]*len(a) if mag == 'nm' else \
                                [0.6 if at.symbol not in ['Co', 'Mn'] else 1 for at in a]
                                n = 'RuO2_%s_dope%s_110_term%s_%s_vac%s_%s%s' \
                                %(Mcomp, dopepos, term, mag, vac_i, adsname, adspos_i)
                                new_f = os.path.join('redo_all_again/%s' %(n))
                                os.mkdir(new_f)
                                adslab.to(os.path.join(new_f, '%s.cif' %(n)))
                                shutil.copyfile('submit_script', os.path.join(new_f, 'submit_script'))
                                metadata = {'structure': adslab.as_dict(), 
                                            'adsites': [site.as_dict() for site in adslab \
                                                        if 'sitetype' in site.properties],
                                            'vac_site': vac_slab.vacsite.as_dict(),
                                            'dopant_sites': [site.as_dict() for site in superslab_dope \
                                                             if site.species_string not in ['O', 'Ru', 'H']]}
                                json.dump(metadata, open(os.path.join(new_f, 'metadata.json'), 'w'))
                                os.chdir(new_f)
                                write_inputs(a, mags=magmoms)
                                os.chdir(maind)

                    for adsname in all_ads_dict.keys():
                        all_adslabs = get_adslabs(superslab, adspos[term], all_ads_dict[adsname])
                        for adspos_i, adslab in enumerate(all_adslabs):

                            a = AseAtomsAdaptor.get_atoms(adslab)
                            a.set_tags([1 if all(sel) else 0 for sel \
                                        in adslab.site_properties['selective_dynamics']])
                            magmoms = [0.6]*len(a) if mag == 'nm' else \
                            [0.6 if at.symbol not in ['Co', 'Mn'] else 1 for at in a]
                            n = 'RuO2_%s_dope%s_110_term%s_%s_%s%s' %(Mcomp, dopepos, term,
                                                                      mag, adsname, adspos_i)
                            new_f = os.path.join('redo_all_again/%s' %(n))
                            os.mkdir(new_f)
                            adslab.to(os.path.join(new_f, '%s.cif' %(n)))
                            shutil.copyfile('submit_script', os.path.join(new_f, 'submit_script'))
                            metadata = {'structure': adslab.as_dict(), 
                                        'adsites': [site.as_dict() for site in adslab \
                                                    if 'sitetype' in site.properties], 
                                        'dopant_sites': [site.as_dict() for site in superslab_dope \
                                                         if site.species_string not in ['O', 'Ru', 'H']]}
                            json.dump(metadata, open(os.path.join(new_f, 'metadata.json'), 'w'))
                            os.chdir(new_f)
                            write_inputs(a, mags=magmoms)
                            os.chdir(maind)
