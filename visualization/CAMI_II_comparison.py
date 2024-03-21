"""
This script is used to reproduce the plot of the CAMI II results(compared to other binners).
"""
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
def get_hq_taxi(dataset = 'skin'):
    if dataset == 'skin':
        data_index = [1,13,14,15,16,17,18,19,20,28]
    else:
        data_index = [6,7,8,13,14,15,16,17,18,19]
    genome_path = os.path.join('Results/Simulated/CAMI_II/{0}/{0}_{1}'.format(dataset,data_index[0]), 'genome')
    method_list = {}
    species_list = {}
    genus_list = {}

    for root, dirs, files in os.walk(genome_path, topdown=False):
        for name in dirs:
            method_list[name] = []
            species_list[name] = []
            genus_list[name] = []

    for temp in data_index:
        taxi = pd.read_csv('Results/Simulated/CAMI_II/{0}/taxonomic_profile.txt'.format(dataset,temp),
                           sep='\t', skiprows=3, dtype={'@@TAXID': str, 'TAXPATH': str})
        taxi_genus = taxi[taxi['RANK'] == 'genus']['@@TAXID'].values.tolist()
        taxi_species = taxi[taxi['RANK'] == 'species']['@@TAXID'].values.tolist()
        method_path_list = []
        genome_path = os.path.join('Results/Simulated/CAMI_II/{0}/{0}_{1}'.format(dataset,temp), 'genome')

        for root, dirs, files in os.walk(genome_path, topdown=False):
            for name in dirs:
                method_path_list.append(os.path.join(root, name))

        for method_path in method_path_list:
            metric = pd.read_csv(os.path.join(method_path, 'metrics_per_bin.tsv'), sep='\t')
            com_90_pur_95 = metric[(metric['Completeness (bp)'].astype(float) > float(0.9)) & (
                    metric['Purity (bp)'].astype(float) >= float(0.95))]
            strain_list = com_90_pur_95['Most abundant genome'].values.tolist()
            method_list[method_path.split('/')[-1]].extend(strain_list)

            for temp_strain in strain_list:
                if temp_strain in taxi['_CAMI_GENOMEID'].values.tolist():
                    taxi_split = taxi[taxi['_CAMI_GENOMEID'] == temp_strain]['TAXPATH'].values[0].split('|')
                    if taxi_split[-2] in taxi_species:
                        species_list[method_path.split('/')[-1]].append(taxi_split[-2])
                    if taxi_split[-3] in taxi_genus:
                        genus_list[method_path.split('/')[-1]].append(taxi_split[-3])

    result = {'Metabat2':{'strain':list(set(method_list['Metabat2_single'])),'species':list(set(species_list['Metabat2_single'])),'genus':list(set(genus_list['Metabat2_single']))},
              'VAMB': {'strain':list(set(method_list['VAMB'])),'species':list(set(species_list['VAMB'])),'genus':list(set(genus_list['VAMB']))},
              'SemiBin': {'strain':list(set(method_list['SemiBin'])),'species':list(set(species_list['SemiBin'])),'genus':list(set(genus_list['SemiBin']))},
              'Metabat2_multi':{'strain':list(set(method_list['Metabat2_multi'])),'species':list(set(species_list['Metabat2_multi'])),'genus':list(set(genus_list['Metabat2_multi']))},
              'SemiBin_c':{'strain':list(set(method_list['SemiBin_c'])),'species':list(set(species_list['SemiBin_c'])),'genus':list(set(genus_list['SemiBin_c']))},
              'NoSemi':{'strain':list(set(method_list['NoSemi'])),'species':list(set(species_list['NoSemi'])),'genus':list(set(genus_list['NoSemi']))},
              'SemiBin_m':{'strain':list(set(method_list['SemiBin_m'])),'species':list(set(species_list['SemiBin_m'])),'genus':list(set(genus_list['SemiBin_m']))},
              'SemiBin_mc':{'strain':list(set(method_list['SemiBin_mc'])),'species':list(set(species_list['SemiBin_mc'])),'genus':list(set(genus_list['SemiBin_mc']))}}

    return result

def plot_hq_num():
    skin_result = get_hq_taxi()
    skin_strain_metabat2 = len(skin_result['Metabat2']['strain'])
    skin_species_metabat2 = len(skin_result['Metabat2']['species'])
    skin_genus_metabat2 = len(skin_result['Metabat2']['genus'])

    skin_strain_vamb= len(skin_result['VAMB']['strain'])
    skin_species_vamb = len(skin_result['VAMB']['species'])
    skin_genus_vamb = len(skin_result['VAMB']['genus'])

    skin_strain_SemiBin = len(skin_result['SemiBin']['strain'])
    skin_species_SemiBin = len(skin_result['SemiBin']['species'])
    skin_genus_SemiBin = len(skin_result['SemiBin']['genus'])

    print(skin_strain_metabat2,skin_species_metabat2,skin_genus_metabat2)
    print(skin_strain_vamb,skin_species_vamb,skin_genus_vamb)
    print(skin_strain_SemiBin,skin_species_SemiBin,skin_genus_SemiBin)

    print('Skin_strain_improvement(Metabat2, VAMB): ', (skin_strain_SemiBin-skin_strain_metabat2)/skin_strain_metabat2, (skin_strain_SemiBin-skin_strain_vamb)/skin_strain_vamb)
    print('Skin_species_improvement(Metabat2, VAMB):', (skin_species_SemiBin-skin_species_metabat2)/skin_species_metabat2, (skin_species_SemiBin-skin_species_vamb)/skin_species_vamb)
    print('Skin_genus_improvement(Metabat2, VAMB):', (skin_genus_SemiBin-skin_genus_metabat2)/skin_genus_metabat2, (skin_genus_SemiBin-skin_genus_vamb)/skin_genus_vamb)

    oral_result = get_hq_taxi(dataset='oral')
    oral_strain_metabat2 = len(oral_result['Metabat2']['strain'])
    oral_species_metabat2 = len(oral_result['Metabat2']['species'])
    oral_genus_metabat2 = len(oral_result['Metabat2']['genus'])

    oral_strain_vamb = len(oral_result['VAMB']['strain'])
    oral_species_vamb = len(oral_result['VAMB']['species'])
    oral_genus_vamb = len(oral_result['VAMB']['genus'])

    oral_strain_SemiBin = len(oral_result['SemiBin']['strain'])
    oral_species_SemiBin = len(oral_result['SemiBin']['species'])
    oral_genus_SemiBin = len(oral_result['SemiBin']['genus'])
    print(oral_strain_metabat2,oral_species_metabat2,oral_genus_metabat2)
    print(oral_strain_vamb,oral_species_vamb,oral_genus_vamb)
    print(oral_strain_SemiBin,oral_species_SemiBin,oral_genus_SemiBin)
    print('oral_strain_improvement(Metabat2, VAMB):', (oral_strain_SemiBin-oral_strain_metabat2)/oral_strain_metabat2, (oral_strain_SemiBin-oral_strain_vamb)/oral_strain_vamb)
    print('oral_species_improvement(Metabat2, VAMB):', (oral_species_SemiBin-oral_species_metabat2)/oral_species_metabat2, (oral_species_SemiBin-oral_species_vamb)/oral_species_vamb)
    print('oral_genus_improvement(Metabat2, VAMB):', (oral_genus_SemiBin-oral_genus_metabat2)/oral_genus_metabat2,  (oral_genus_SemiBin-oral_genus_vamb)/oral_genus_vamb)

    line_width = 1

    plt.figure(figsize=(4, 2))
    plt.plot(['genus', 'species', 'strain'],[skin_genus_metabat2, skin_species_metabat2, skin_strain_metabat2], label='Metabat2',color='#ec7014',linewidth = line_width, marker='o',)


    plt.plot(['genus', 'species', 'strain'],[skin_genus_vamb, skin_species_vamb, skin_strain_vamb], label='VAMB',color='#7570b3',linewidth = line_width, marker='o',)
    plt.plot(['genus', 'species', 'strain'], [skin_genus_SemiBin, skin_species_SemiBin, skin_strain_SemiBin], label='SemiBin',
             color='#1b9e77',linewidth = line_width, marker='o',)
    plt.legend()
    plt.xticks([])
    plt.savefig('CAMI_II_skin.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(4, 2))
    plt.plot(['genus', 'species', 'strain'],[oral_genus_metabat2, oral_species_metabat2, oral_strain_metabat2], label='Metabat2',color='#ec7014',linewidth = line_width, marker='o',)
    plt.plot(['genus', 'species', 'strain'],[oral_genus_vamb, oral_species_vamb, oral_strain_vamb], label='VAMB',color='#7570b3',linewidth = line_width, marker='o',)
    plt.plot(['genus', 'species', 'strain'], [oral_genus_SemiBin, oral_species_SemiBin, oral_strain_SemiBin], label='SemiBin',
             color='#1b9e77',linewidth = line_width, marker='o',)
    plt.legend()
    plt.savefig('CAMI_II_oral.pdf', dpi=300, bbox_inches='tight')
    plt.close()

def plot_venn_plot():
    skin_result = get_hq_taxi()
    skin_strain_metabat2 = skin_result['Metabat2']['strain']
    skin_species_metabat2 = skin_result['Metabat2']['species']
    skin_genus_metabat2 = skin_result['Metabat2']['genus']

    skin_strain_vamb= skin_result['VAMB']['strain']
    skin_species_vamb = skin_result['VAMB']['species']
    skin_genus_vamb = skin_result['VAMB']['genus']

    skin_strain_SemiBin = skin_result['SemiBin']['strain']
    skin_species_SemiBin = skin_result['SemiBin']['species']
    skin_genus_SemiBin = skin_result['SemiBin']['genus']

    oral_result = get_hq_taxi(dataset='oral')
    oral_strain_metabat2 = oral_result['Metabat2']['strain']
    oral_species_metabat2 = oral_result['Metabat2']['species']
    oral_genus_metabat2 = oral_result['Metabat2']['genus']

    oral_strain_vamb = oral_result['VAMB']['strain']
    oral_species_vamb = oral_result['VAMB']['species']
    oral_genus_vamb = oral_result['VAMB']['genus']

    oral_strain_SemiBin = oral_result['SemiBin']['strain']
    oral_species_SemiBin = oral_result['SemiBin']['species']
    oral_genus_SemiBin = oral_result['SemiBin']['genus']

    from matplotlib_venn import venn3_unweighted
    out = venn3_unweighted([set(skin_strain_metabat2),set(skin_strain_vamb),set(skin_strain_SemiBin)],set_labels=('Metabat2', 'VAMB','SemiBin'),normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("strain(Skin)", fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_Strain_Skin.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(skin_species_metabat2),set(skin_species_vamb),set(skin_species_SemiBin)],set_labels=('Metabat2', 'VAMB','SemiBin'),normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Skin)", fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_Species_Skin.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(skin_genus_metabat2),set(skin_genus_vamb),set(skin_genus_SemiBin)],set_labels=('Metabat2', 'VAMB','SemiBin'),normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Skin)", fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_Genus_Skin.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(oral_strain_metabat2),set(oral_strain_vamb),set(oral_strain_SemiBin)],set_labels=('Metabat2', 'VAMB','SemiBin'),normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("strain(Oral)", fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_Strain_Oral.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(oral_species_metabat2),set(oral_species_vamb),set(oral_species_SemiBin)],set_labels=('Metabat2', 'VAMB','SemiBin'),normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Oral)", fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_Species_Oral.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(oral_genus_metabat2),set(oral_genus_vamb),set(oral_genus_SemiBin)],set_labels=('Metabat2', 'VAMB','SemiBin'),normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Oral)", fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_Genus_Oral.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

def get_similarity_distribution(dataset = 'skin'):
    if dataset == 'skin':
        data_index = [1,13,14,15,16,17,18,19,20,28]
    else:
        data_index = [6,7,8,13,14,15,16,17,18,19]

    genome_to_id =  pd.read_csv('Results/Simulated/CAMI_II/{0}/genome_to_id.tsv'.format(dataset),sep='\t',header=None).values.tolist()
    genome_id = {}
    for temp in genome_to_id:
        genome_id[temp[1].split('/')[-1]] = temp[0]

    distribution = {}

    for sample in data_index:
        distribution[sample] = {90:[],95:[],97:[],98:[],99:[],99.5:[],99.9:[],100:[]}
        genome_ANI = pd.read_csv('Results/Simulated/CAMI_II/{0}/genome_{1}_ANI'.format(dataset,sample), header=None,
                                 sep='\t')
        genome_ANI.columns = ['genome_id_1', 'genome_id_2', 'ANI', 'A', 'B']
        genome_group = genome_ANI.groupby('genome_id_1')
        for genome, group in genome_group:
            group = group[group.genome_id_2 != genome]
            ANI = group['ANI'].values.tolist()
            if ANI != []:
                ANI_max = np.max(ANI)
            else:
                ANI_max = 0
            if ANI_max <= 90:
                distribution[sample][90].append(genome_id[genome.split('/')[-1]])
            elif ANI_max <= 95:
                distribution[sample][95].append(genome_id[genome.split('/')[-1]])
            elif ANI_max <= 97:
                distribution[sample][97].append(genome_id[genome.split('/')[-1]])
            elif ANI_max <= 98:
                distribution[sample][98].append(genome_id[genome.split('/')[-1]])
            elif ANI_max <= 99:
                distribution[sample][99].append(genome_id[genome.split('/')[-1]])
            elif ANI_max <= 99.5:
                distribution[sample][99.5].append(genome_id[genome.split('/')[-1]])
            elif ANI_max <= 99.9:
                distribution[sample][99.9].append(genome_id[genome.split('/')[-1]])
            else:
                distribution[sample][100].append(genome_id[genome.split('/')[-1]])

    return distribution

def plot_similarity_distribution():
    skin_distribution = get_similarity_distribution()
    x = ['0~90','90~95','95~97','97~98','98~99','99~99.5','99.5~99.9','99.9~100']
    y = {}
    for sample in [1,13,14,15,16,17,18,19,20,28]:
        y[sample] = []
        for value in skin_distribution[sample]:
            y[sample].append(len(skin_distribution[sample][value]))

    plt.plot(x,y[1], label='2017.12.04_18.56.22_sample_1',color='#e6ab02')
    plt.plot(x, y[13], label='2017.12.04_18.56.22_sample_13',color='#e7298a')
    plt.plot(x, y[14], label='2017.12.04_18.56.22_sample_14',color='#7570b3')
    plt.plot(x, y[15], label='2017.12.04_18.56.22_sample_15',color='#d95f02')
    plt.plot(x, y[16], label='2017.12.04_18.56.22_sample_16',color='#1b9e77')
    plt.xticks(rotation=50,color = 'black')
    plt.yticks(color='black')
    plt.legend()
    plt.title('Genome distribution of Skin',fontsize = 15)

    plt.xlabel('ANI',fontsize = 15,color = 'black')
    plt.ylabel('Num',fontsize = 15,color = 'black')
    plt.savefig('CAMI_II_distribution_Skin_1.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()


    plt.plot(x,y[17], label='2017.12.04_18.56.22_sample_17',color='#e6ab02')
    plt.plot(x, y[18], label='2017.12.04_18.56.22_sample_18',color='#e7298a')
    plt.plot(x, y[19], label='2017.12.04_18.56.22_sample_19',color='#7570b3')
    plt.plot(x, y[20], label='2017.12.04_18.56.22_sample_20',color='#d95f02')
    plt.plot(x, y[28], label='2017.12.04_18.56.22_sample_28',color='#1b9e77')
    plt.xticks(rotation=50,color = 'black')
    plt.yticks(color='black')
    plt.legend()
    plt.title('Genome distribution of Skin',fontsize = 15)

    plt.xlabel('ANI',fontsize = 15,color = 'black')
    plt.ylabel('Num',fontsize = 15,color = 'black')
    plt.savefig('CAMI_II_distribution_Skin_2.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    oral_distribution = get_similarity_distribution(dataset='oral')

    x = ['0~90','90~95','95~97','97~98','98~99','99~99.5','99.5~99.9','99.9~100']
    y = {}
    for sample in  [6,7,8,13,14,15,16,17,18,19]:
        y[sample] = []
        for value in oral_distribution[sample]:
            y[sample].append(len(oral_distribution[sample][value]))

    plt.plot(x,y[6], label='2017.12.04_18.56.22_sample_6',color='#e6ab02')
    plt.plot(x, y[7], label='2017.12.04_18.56.22_sample_7',color='#e7298a')
    plt.plot(x, y[8], label='2017.12.04_18.56.22_sample_8',color='#7570b3')
    plt.plot(x, y[13], label='2017.12.04_18.56.22_sample_13',color='#d95f02')
    plt.plot(x, y[14], label='2017.12.04_18.56.22_sample_14',color='#1b9e77')
    plt.xticks(rotation=50,color = 'black')
    plt.yticks(color='black')
    plt.legend()
    plt.title('Genome distribution of Oral',fontsize = 15)

    plt.xlabel('ANI',fontsize = 15,color = 'black')
    plt.ylabel('Num',fontsize = 15,color = 'black')
    plt.savefig('CAMI_II_distribution_Oral_1.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()


    plt.plot(x,y[15], label='2017.12.04_18.56.22_sample_15',color='#e6ab02')
    plt.plot(x, y[16], label='2017.12.04_18.56.22_sample_16',color='#e7298a')
    plt.plot(x, y[17], label='2017.12.04_18.56.22_sample_17',color='#7570b3')
    plt.plot(x, y[18], label='2017.12.04_18.56.22_sample_18',color='#d95f02')
    plt.plot(x, y[19], label='2017.12.04_18.56.22_sample_19',color='#1b9e77')
    plt.xticks(rotation=50,color = 'black')
    plt.yticks(color='black')
    plt.legend()
    plt.title('Genome distribution of Oral',fontsize = 15)

    plt.xlabel('ANI',fontsize = 15,color = 'black')
    plt.ylabel('Num',fontsize = 15,color = 'black')
    plt.savefig('CAMI_II_distribution_Oral_2.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

def get_hq_strain_ANI(dataset = 'skin'):
    if dataset == 'skin':
        data_index = [1,13,14,15,16,17,18,19,20,28]
    else:
        data_index = [6,7,8,13,14,15,16,17,18,19]

    distribution = get_similarity_distribution(dataset)
    genome_path = os.path.join('Results/Simulated/CAMI_II/{0}/{0}_{1}'.format(dataset,data_index[0]), 'genome')
    method_list = {}

    for root, dirs, files in os.walk(genome_path, topdown=False):
        for name in dirs:
            method_list[os.path.join(root, name).split('/')[-1]] = {90: [], 95: [], 97: [], 98: [], 99: [],
                                                                          99.5: [], 99.9: [], 100: []}

    for temp in data_index:
        method_path_list = []
        genome_path = os.path.join('Results/Simulated/CAMI_II/{0}/{0}_{1}'.format(dataset, temp), 'genome')
        for root, dirs, files in os.walk(genome_path, topdown=False):
            for name in dirs:
                method_path_list.append(os.path.join(root, name))

        for method_path in method_path_list:
            metric = pd.read_csv(os.path.join(method_path, 'metrics_per_bin.tsv'), sep='\t')
            com_90_pur_95 = metric[(metric['Completeness (bp)'].astype(float) > float(0.9)) & (
                    metric['Purity (bp)'].astype(float) >= float(0.95))]
            strain_list = com_90_pur_95['Most abundant genome'].values.tolist()
            for strain in strain_list:
                for value in [90, 95, 97, 98, 99, 99.5, 99.9, 100]:
                    if strain in distribution[temp][value]:
                        method_list[method_path.split('/')[-1]][value].append(strain)

    value_list = [90, 95, 97, 98, 99, 99.5, 99.9, 100]
    value_array = np.zeros(shape=(8, 3))

    for i in range(len(value_list)):
        value_array[i][0] = len(set(method_list['Metabat2_single'][value_list[i]]))
        value_array[i][1] = len(set(method_list['VAMB'][value_list[i]]))
        value_array[i][2] = len(set(method_list['SemiBin'][value_list[i]]))

    return value_array

def plot_bar_strain_simiarity():
    ### skin
    value_array = get_hq_strain_ANI()
    subset = pd.DataFrame(value_array,columns = ['Metabat2','VAMB','SemiBin'], index=['0~90','90~95','95~97','97~98','98~99','99~99.5','99.5~99.9','99.9~100'])
    print(subset)
    ax = subset.plot(kind='bar',width = 0.6, color = ['#ec7014', '#7570b3', '#1b9e77'])
    ax.set_yticks(ticks=[0,20,40,60,80])
    ax.set_yticklabels(labels=[0,20,40,60,80],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['0~90','90~95','95~97','97~98','98~99','99~99.5','99.5~99.9','99.9~100'], fontsize=15,color = 'black',rotation = 50)
    ax.set_ylabel('High quality genomes', fontsize=15,color = 'black')
    ax.set_title('Skin', fontsize=15, alpha=1.0)
    plt.savefig('CAMI_II_distribution_Skin_bar.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    ### Oral
    value_array = get_hq_strain_ANI(dataset= 'oral')
    subset = pd.DataFrame(value_array,columns = ['Metabat2','VAMB','SemiBin'], index=['0~90','90~95','95~97','97~98','98~99','99~99.5','99.5~99.9','99.9~100'])
    print(subset)
    ax = subset.plot(kind='bar',width = 0.6, color = ['#ec7014', '#7570b3', '#1b9e77'])
    ax.set_yticks(ticks=[0,30,60,90])
    ax.set_yticklabels(labels=[0,30,60,90],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['0~90','90~95','95~97','97~98','98~99','99~99.5','99.5~99.9','99.9~100'], fontsize=15,color = 'black',rotation = 50)
    ax.set_ylabel('High quality genomes', fontsize=15,color = 'black')
    ax.set_title('Oral', fontsize=15, alpha=1.0)
    plt.savefig('CAMI_II_distribution_Oral_bar.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

def plot_generalization():
    result = get_hq_taxi('skin')
    skin_SemiBin_strain = len(result['SemiBin']['strain'])
    skin_SemiBin_species = len(result['SemiBin']['species'])
    skin_SemiBin_genus = len(result['SemiBin']['genus'])

    skin_NoSemi_strain = len(result['NoSemi']['strain'])
    skin_NoSemi_species = len(result['NoSemi']['species'])
    skin_NoSemi_genus = len(result['NoSemi']['genus'])

    skin_SemiBin_m_strain = len(result['SemiBin_m']['strain'])
    skin_SemiBin_m_species = len(result['SemiBin_m']['species'])
    skin_SemiBin_m_genus = len(result['SemiBin_m']['genus'])

    skin_SemiBin_c_strain = len(result['SemiBin_c']['strain'])
    skin_SemiBin_c_species = len(result['SemiBin_c']['species'])
    skin_SemiBin_c_genus = len(result['SemiBin_c']['genus'])

    skin_SemiBin_mc_strain = len(result['SemiBin_mc']['strain'])
    skin_SemiBin_mc_species = len(result['SemiBin_mc']['species'])
    skin_SemiBin_mc_genus = len(result['SemiBin_mc']['genus'])
    print('Skin')
    print('Semibin compared to NoSemi {:.2%} improvement'.format( (skin_SemiBin_strain - skin_NoSemi_strain)/ skin_NoSemi_strain))
    print('Semibin compared to SemiBin_m {:.2%} improvement'.format( (skin_SemiBin_strain - skin_SemiBin_m_strain)/ skin_SemiBin_m_strain))
    print('Semibin compared to SemiBin_c {:.2%} improvement'.format( (skin_SemiBin_strain - skin_SemiBin_c_strain)/ skin_SemiBin_c_strain))
    print('Semibin compared to SemiBin_mc {:.2%} improvement'.format((skin_SemiBin_strain - skin_SemiBin_mc_strain) / skin_SemiBin_mc_strain))
    print('average {:.2%} improvement'.format(
        (((skin_SemiBin_strain - skin_SemiBin_m_strain)/ skin_SemiBin_m_strain) +
        ((skin_SemiBin_strain - skin_SemiBin_c_strain)/ skin_SemiBin_c_strain) +
        ((skin_SemiBin_strain - skin_SemiBin_mc_strain) / skin_SemiBin_mc_strain))/3))

    result = get_hq_taxi('oral')
    oral_SemiBin_strain = len(result['SemiBin']['strain'])
    oral_SemiBin_species = len(result['SemiBin']['species'])
    oral_SemiBin_genus = len(result['SemiBin']['genus'])

    oral_NoSemi_strain = len(result['NoSemi']['strain'])
    oral_NoSemi_species = len(result['NoSemi']['species'])
    oral_NoSemi_genus = len(result['NoSemi']['genus'])

    oral_SemiBin_m_strain = len(result['SemiBin_m']['strain'])
    oral_SemiBin_m_species = len(result['SemiBin_m']['species'])
    oral_SemiBin_m_genus = len(result['SemiBin_m']['genus'])

    oral_SemiBin_c_strain = len(result['SemiBin_c']['strain'])
    oral_SemiBin_c_species = len(result['SemiBin_c']['species'])
    oral_SemiBin_c_genus = len(result['SemiBin_c']['genus'])

    oral_SemiBin_mc_strain = len(result['SemiBin_mc']['strain'])
    oral_SemiBin_mc_species = len(result['SemiBin_mc']['species'])
    oral_SemiBin_mc_genus = len(result['SemiBin_mc']['genus'])

    print('Oral')
    print('Semibin compared to NoSemi {:.2%} improvement'.format( (oral_SemiBin_strain - oral_NoSemi_strain)/ oral_NoSemi_strain))
    print('Semibin compared to SemiBin_m {:.2%} improvement'.format( (oral_SemiBin_strain - oral_SemiBin_m_strain)/ oral_SemiBin_m_strain))
    print('Semibin compared to SemiBin_c {:.2%} improvement'.format( (oral_SemiBin_strain - oral_SemiBin_c_strain)/ oral_SemiBin_c_strain))
    print('Semibin compared to SemiBin_mc {:.2%} improvement'.format((oral_SemiBin_strain - oral_SemiBin_mc_strain) / oral_SemiBin_mc_strain))
    print('average {:.2%} improvement'.format(
        (((oral_SemiBin_strain - oral_SemiBin_m_strain)/ oral_SemiBin_m_strain) +
        ((oral_SemiBin_strain - oral_SemiBin_c_strain)/ oral_SemiBin_c_strain) +
        ((oral_SemiBin_strain - oral_SemiBin_mc_strain) / oral_SemiBin_mc_strain))/3))

    subset = pd.DataFrame(np.array([[skin_NoSemi_strain,skin_SemiBin_m_strain,skin_SemiBin_c_strain,skin_SemiBin_mc_strain,skin_SemiBin_strain ],[oral_NoSemi_strain,oral_SemiBin_m_strain,oral_SemiBin_c_strain,oral_SemiBin_mc_strain,oral_SemiBin_strain]]),columns = ['No_semi','SemiBin_m','SemiBin_c','SemiBin_mc','SemiBin'], index=['Skin','Oral'])
    print(subset)
    ax = subset.plot(kind='bar',width = 0.6,color = ['#e6ab02','#e7298a', '#7570b3', '#d95f02', '#1b9e77'])
    ax.set_yticks(ticks=[0,20,40,60,80,100,120,140])
    ax.set_yticklabels(labels=[0,20,40,60,80,100,120,140],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality strains', fontsize=15,color = 'black')
    ax.set_title('strain', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_com_strain_generalization.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    subset = pd.DataFrame(np.array([[skin_NoSemi_species,skin_SemiBin_m_species,skin_SemiBin_c_species,skin_SemiBin_mc_species,skin_SemiBin_species ],[oral_NoSemi_species,oral_SemiBin_m_species,oral_SemiBin_c_species,oral_SemiBin_mc_species,oral_SemiBin_species]]),columns = ['No_semi','SemiBin_m','SemiBin_c','SemiBin_mc','SemiBin'], index=['Skin','Oral'])
    print(subset)
    ax = subset.plot(kind='bar',width = 0.6,color = ['#e6ab02','#e7298a', '#7570b3', '#d95f02', '#1b9e77'])
    ax.set_yticks(ticks=[0,15,30,45,60,75,90,105])
    ax.set_yticklabels(labels=[0,15,30,45,60,75,90,105],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality species', fontsize=15,color = 'black')
    ax.set_title('species', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_com_species_generalization.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    subset = pd.DataFrame(np.array([[skin_NoSemi_genus,skin_SemiBin_m_genus,skin_SemiBin_c_genus,skin_SemiBin_mc_genus,skin_SemiBin_genus ],[oral_NoSemi_genus,oral_SemiBin_m_genus,oral_SemiBin_c_genus,oral_SemiBin_mc_genus,oral_SemiBin_genus]]),columns = ['No_semi','SemiBin_m','SemiBin_c','SemiBin_mc','SemiBin'], index=['Skin','Oral'])
    print(subset)
    ax = subset.plot(kind='bar',width = 0.6,color = ['#e6ab02','#e7298a', '#7570b3', '#d95f02', '#1b9e77'])
    ax.set_yticks(ticks=[0,10,20,30,40,50,60])
    ax.set_yticklabels(labels=[0,10,20,30,40,50,60],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality genus', fontsize=15,color = 'black')
    ax.set_title('genus', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_com_genus_generalization.pdf', dpi=300, bbox_inches='tight')
    plt.close()

def plot_comparsion_Metabat2():
    result = get_hq_taxi('skin')

    Metabat2_single_strain_skin = len(result['Metabat2']['strain'])
    Metabat2_single_species_skin = len(result['Metabat2']['species'])
    Metabat2_single_genus_skin = len(result['Metabat2']['genus'])

    Metabat2_multi_strain_skin = len(result['Metabat2_multi']['strain'])
    Metabat2_multi_species_skin = len(result['Metabat2_multi']['species'])
    Metabat2_multi_genus_skin = len(result['Metabat2_multi']['genus'])

    result = get_hq_taxi('oral')

    Metabat2_single_strain_oral = len(result['Metabat2']['strain'])
    Metabat2_single_species_oral = len(result['Metabat2']['species'])
    Metabat2_single_genus_oral = len(result['Metabat2']['genus'])

    Metabat2_multi_strain_oral = len(result['Metabat2_multi']['strain'])
    Metabat2_multi_species_oral = len(result['Metabat2_multi']['species'])
    Metabat2_multi_genus_oral = len(result['Metabat2_multi']['genus'])

    print(Metabat2_single_strain_skin, Metabat2_single_species_skin, Metabat2_single_genus_skin)
    print(Metabat2_multi_strain_skin, Metabat2_multi_species_skin, Metabat2_multi_genus_skin)

    print(Metabat2_single_strain_oral, Metabat2_single_species_oral, Metabat2_single_genus_oral)
    print(Metabat2_multi_strain_oral, Metabat2_multi_species_oral, Metabat2_multi_genus_oral)

    subset = pd.DataFrame(np.array([[Metabat2_single_strain_skin,Metabat2_multi_strain_skin],[Metabat2_single_strain_oral,Metabat2_multi_strain_oral]]),columns = ['Metabat2_single','Metabat2_multi'], index=['Skin','Oral'])
    print(subset)
    ax = subset.plot(kind='bar',width = 0.8,color = ['#ec7014', '#fec44f'],figsize=(3,4))
    ax.set_yticks(ticks=[0,20,40,60,80])
    ax.set_yticklabels(labels=[0,20,40,60,80],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality strains', fontsize=15,color = 'black')
    ax.set_title('strain', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_com_strain_Metabat2.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    subset = pd.DataFrame(np.array([[Metabat2_single_species_skin,Metabat2_multi_species_skin],[Metabat2_single_species_oral,Metabat2_multi_species_oral]]),columns = ['Metabat2_single','Metabat2_multi'], index=['Skin','Oral'])
    print(subset)
    ax = subset.plot(kind='bar',width= 0.8,legend = False,color = ['#ec7014', '#fec44f'],figsize=(3,4))
    ax.set_yticks(ticks=[0,15,30,45,60,75])
    ax.set_yticklabels(labels=[0,15,30,45,60,75],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    ax.set_title('species', fontsize=20, alpha=1.0,color = 'black')
    ax.set_ylabel('High quality species', fontsize=15,color = 'black')
    plt.savefig('CAMI_II_com_species_Metabat2.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    subset = pd.DataFrame(np.array([[Metabat2_single_genus_skin, Metabat2_multi_genus_skin],[Metabat2_single_genus_oral, Metabat2_multi_genus_oral]]),columns = ['Metabat2_single','Metabat2_multi'], index=['Skin','Oral'])
    print(subset)
    ax = subset.plot(kind='bar',width = 0.8,legend = False,color = ['#ec7014', '#fec44f'],figsize=(3,4))
    ax.set_yticks(ticks=[0,10,20,30,40,50])
    ax.set_yticklabels(labels=[0,10,20,30,40,50],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    ax.set_title('genus', fontsize=20, alpha=1.0,color = 'black')
    ax.legend(['Metabat2_single','Metabat2_multi'],fontsize=10)
    ax.set_ylabel('High quality genera', fontsize=15,color = 'black')
    plt.savefig('CAMI_II_com_genus_Metabat2.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

def get_result(dataset='skin', input = None):
    if dataset == 'skin':
        data_index = [1, 13, 14, 15, 16, 17, 18, 19, 20, 28]
    else:
        data_index = [6, 7, 8, 13, 14, 15, 16, 17, 18, 19]
    genome_path = os.path.join('{0}/{1}/{1}_{2}'.format(input, dataset,data_index[0]),'genome')
    method_list = {}
    species_list = {}
    genus_list = {}

    for root, dirs, files in os.walk(genome_path, topdown=False):
        for name in dirs:
            method_list[name] = []
            species_list[name] = []
            genus_list[name] = []

    for temp in data_index:
        taxi = pd.read_csv(
            'Results/Simulated/CAMI_II/{0}/taxonomic_profile.txt'.format(dataset, temp), sep='\t', skiprows=3, dtype={'@@TAXID': str, 'TAXPATH': str})
        taxi_genus = taxi[taxi['RANK'] == 'genus']['@@TAXID'].values.tolist()
        taxi_species = taxi[taxi['RANK'] == 'species']['@@TAXID'].values.tolist()
        method_path_list = []
        genome_path = os.path.join('{0}/{1}/{1}_{2}'.format(input, dataset, temp), 'genome')

        for root, dirs, files in os.walk(genome_path, topdown=False):
            for name in dirs:
                method_path_list.append(os.path.join(root, name))

        for method_path in method_path_list:
            metric = pd.read_csv(os.path.join(method_path, 'metrics_per_bin.tsv'), sep='\t')
            com_90_pur_95 = metric[(metric['Completeness (bp)'].astype(float) > float(0.9)) & (metric['Purity (bp)'].astype(float) >= float(0.95))]
            strain_list = com_90_pur_95['Most abundant genome'].values.tolist()
            method_list[method_path.split('/')[-1]].extend(strain_list)

            for temp_strain in strain_list:
                if temp_strain in taxi['_CAMI_GENOMEID'].values.tolist():
                    taxi_split = taxi[taxi['_CAMI_GENOMEID'] == temp_strain]['TAXPATH'].values[0].split('|')
                    if taxi_split[-2] in taxi_species:
                        species_list[method_path.split('/')[-1]].append(taxi_split[-2])
                    else:
                        print(method_path, 'error1')
                    if taxi_split[-3] in taxi_genus:
                        genus_list[method_path.split('/')[-1]].append(taxi_split[-3])
                    else:
                        print(method_path, 'error2')
    return method_list, species_list, genus_list

def plot_recluster():
    def get_recluster(dataset, input = None):
        method_list, species_list, genus_list = get_result(dataset, input)
        result = {'No_recluster': {
            'strain': list(set(method_list['SemiBin_no_relcluster'])),
            'species': list(set(species_list['SemiBin_no_relcluster'])),
            'genus': list(set(genus_list['SemiBin_no_relcluster']))},
                  'SemiBin': {'strain': list(set(method_list['SemiBin'])),
                              'species': list(set(species_list['SemiBin'])),
                              'genus': list(set(genus_list['SemiBin']))}, }
        return result
    skin_result = get_recluster('skin', 'updated_results/effect_reclustering')
    skin_strain_no_recluster = len(skin_result['No_recluster']['strain'])
    skin_species_no_recluster = len(skin_result['No_recluster']['species'])
    skin_genus_no_recluster = len(skin_result['No_recluster']['genus'])

    skin_strain_SemiBin = len(skin_result['SemiBin']['strain'])
    skin_species_SemiBin = len(skin_result['SemiBin']['species'])
    skin_genus_SemiBin = len(skin_result['SemiBin']['genus'])
    print(skin_strain_no_recluster, skin_species_no_recluster, skin_genus_no_recluster)
    print(skin_strain_SemiBin, skin_species_SemiBin, skin_genus_SemiBin)

    oral_result = get_recluster('oral', 'updated_results/effect_reclustering')
    oral_strain_no_recluster = len(oral_result['No_recluster']['strain'])
    oral_species_no_recluster = len(oral_result['No_recluster']['species'])
    oral_genus_no_recluster = len(oral_result['No_recluster']['genus'])

    oral_strain_SemiBin = len(oral_result['SemiBin']['strain'])
    oral_species_SemiBin = len(oral_result['SemiBin']['species'])
    oral_genus_SemiBin = len(oral_result['SemiBin']['genus'])
    print(oral_strain_no_recluster, oral_species_no_recluster, oral_genus_no_recluster)
    print(oral_strain_SemiBin, oral_species_SemiBin, oral_genus_SemiBin)

    line_width = 1

    plt.plot(['genus', 'species', 'strain'],
             [skin_genus_no_recluster, skin_species_no_recluster, skin_strain_no_recluster],
             label='No recluster', color='#ec7014', linewidth=line_width,
             marker='o', )
    plt.plot(['genus', 'species', 'strain'],
             [skin_genus_SemiBin, skin_species_SemiBin, skin_strain_SemiBin],
             label='SemiBin',
             color='#1b9e77', linewidth=line_width, marker='o', )
    plt.legend()
    plt.title('Skin', fontsize=15, alpha=1.0)
    plt.savefig('CAMI_II_skin_recluster.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    plt.plot(['genus', 'species', 'strain'],
             [oral_genus_no_recluster, oral_species_no_recluster, oral_strain_no_recluster],
             label='No recluster', color='#ec7014', linewidth=line_width,
             marker='o', )
    plt.plot(['genus', 'species', 'strain'],
             [oral_genus_SemiBin, oral_species_SemiBin, oral_strain_SemiBin],
             label='SemiBin',
             color='#1b9e77', linewidth=line_width, marker='o', )
    plt.legend()
    plt.title('Oral', fontsize=15, alpha=1.0)
    plt.savefig('CAMI_II_oral_recluster.pdf', dpi=300, bbox_inches='tight')
    plt.close()

def plot_cluster_alternative():
    def get_cluster(dataset='skin', input = None):
        method_list, species_list, genus_list = get_result(dataset, input)
        result = {'lp': {'strain': list(set(method_list['lp'])),'species': list(set(species_list['lp'])),'genus': list(set(genus_list['lp']))},
                  'leiden': {'strain': list(set(method_list['leiden'])),'species': list(set(species_list['leiden'])),'genus': list(set(genus_list['leiden']))},
                  'multi_level': {'strain': list(set(method_list['multi_level'])),'species': list(set(species_list['multi_level'])),'genus': list(set(genus_list['multi_level']))},
                  'infomap': {'strain': list(set(method_list['infomap'])),'species': list(set(species_list['infomap'])),'genus': list(set(genus_list['infomap']))}, }
        return result
    skin_result = get_cluster('skin', 'updated_results/effect_clustering/cluster')
    skin_strain_lp = len(skin_result['lp']['strain'])
    skin_species_lp = len(skin_result['lp']['species'])
    skin_genus_lp = len(skin_result['lp']['genus'])

    skin_strain_leiden = len(skin_result['leiden']['strain'])
    skin_species_leiden = len(skin_result['leiden']['species'])
    skin_genus_leiden = len(skin_result['leiden']['genus'])

    skin_strain_multi_level = len(skin_result['multi_level']['strain'])
    skin_species_multi_level = len(skin_result['multi_level']['species'])
    skin_genus_multi_level = len(skin_result['multi_level']['genus'])

    skin_strain_infomap = len(skin_result['infomap']['strain'])
    skin_species_infomap = len(skin_result['infomap']['species'])
    skin_genus_infomap = len(skin_result['infomap']['genus'])

    print(skin_strain_lp, skin_species_lp, skin_genus_lp)
    print(skin_strain_leiden, skin_species_leiden, skin_genus_leiden)
    print(skin_strain_multi_level, skin_species_multi_level, skin_genus_multi_level)
    print(skin_strain_infomap, skin_species_infomap, skin_genus_infomap)

    oral_result = get_cluster('oral', 'updated_results/effect_clustering/cluster')
    oral_strain_lp = len(oral_result['lp']['strain'])
    oral_species_lp = len(oral_result['lp']['species'])
    oral_genus_lp = len(oral_result['lp']['genus'])

    oral_strain_leiden = len(oral_result['leiden']['strain'])
    oral_species_leiden = len(oral_result['leiden']['species'])
    oral_genus_leiden = len(oral_result['leiden']['genus'])

    oral_strain_multi_level = len(oral_result['multi_level']['strain'])
    oral_species_multi_level = len(oral_result['multi_level']['species'])
    oral_genus_multi_level = len(oral_result['multi_level']['genus'])

    oral_strain_infomap = len(oral_result['infomap']['strain'])
    oral_species_infomap = len(oral_result['infomap']['species'])
    oral_genus_infomap = len(oral_result['infomap']['genus'])
    
    print(oral_strain_lp, oral_species_lp, oral_genus_lp)
    print(oral_strain_leiden, oral_species_leiden, oral_genus_leiden)
    print(oral_strain_multi_level, oral_species_multi_level, oral_genus_multi_level)
    print(oral_strain_infomap, oral_species_infomap, oral_genus_infomap)

    line_width = 1

    plt.plot(['Genus', 'Species', 'Strain'],
             [skin_genus_lp, skin_species_lp,
              skin_strain_lp],
             label='Label propagation', color='#ec7014', linewidth=line_width,
             marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [skin_genus_multi_level, skin_species_multi_level, skin_strain_multi_level],
             label='Louvain',
             color='#7570b3', linewidth=line_width, marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [skin_genus_leiden, skin_species_leiden, skin_strain_leiden],
             label='Leiden',
             color='#e7298a', linewidth=line_width, marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [skin_genus_infomap, skin_species_infomap, skin_strain_infomap],
             label='Infomap',
             color='#1b9e77', linewidth=line_width, marker='o', )
    plt.legend()
    # plt.xticks(['Strain', 'Species', 'Genus'])
    plt.title('Skin', fontsize=15, alpha=1.0)
    plt.savefig('CAMI_II_skin_cluster_alternative.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    plt.plot(['Genus', 'Species', 'Strain'],
             [oral_genus_lp, oral_species_lp, oral_strain_lp],
             label='Label propagation', color='#ec7014', linewidth=line_width,
             marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [oral_genus_multi_level, oral_species_multi_level, oral_strain_multi_level],
             label='Louvain',
             color='#7570b3', linewidth=line_width, marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [oral_genus_leiden, oral_species_leiden, oral_strain_leiden],
             label='Leiden',
             color='#e7298a', linewidth=line_width, marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [oral_genus_infomap, oral_species_infomap, oral_strain_infomap],
             label='Infomap',
             color='#1b9e77', linewidth=line_width, marker='o', )
    plt.legend()
    plt.title('Oral', fontsize=15, alpha=1.0)
    # plt.xticks(['Strain', 'Species', 'Genus'])
    plt.savefig('CAMI_II_oral_cluster_alternative.pdf', dpi=300, bbox_inches='tight')
    plt.close()

def plot_max_edge():
    def get_cluster(dataset='skin', input=None):
        method_list, species_list, genus_list = get_result(dataset, input)
        result = {'200': {'strain': list(set(method_list['200'])),
                         'species': list(set(species_list['200'])),
                         'genus': list(set(genus_list['200']))},
                  '500': {'strain': list(set(method_list['500'])),
                             'species': list(set(species_list['500'])),
                             'genus': list(set(genus_list['500']))},
                  '1000': {'strain': list(set(method_list['1000'])),
                      'species': list(set(species_list['1000'])),
                      'genus': list(set(genus_list['1000']))}, }
        return result

    skin_result = get_cluster('skin', 'updated_results/effect_max_edge')
    skin_strain_200 = len(skin_result['200']['strain'])
    skin_species_200 = len(skin_result['200']['species'])
    skin_genus_200 = len(skin_result['200']['genus'])

    skin_strain_500 = len(skin_result['500']['strain'])
    skin_species_500 = len(skin_result['500']['species'])
    skin_genus_500 = len(skin_result['500']['genus'])

    skin_strain_1000 = len(skin_result['1000']['strain'])
    skin_species_1000 = len(skin_result['1000']['species'])
    skin_genus_1000 = len(skin_result['1000']['genus'])

    print(skin_strain_200, skin_species_200, skin_genus_200)
    print(skin_strain_500, skin_species_500, skin_genus_500)
    print(skin_strain_1000, skin_species_1000, skin_genus_1000)


    oral_result = get_cluster('oral', 'updated_results/effect_max_edge')
    oral_strain_200 = len(oral_result['200']['strain'])
    oral_species_200 = len(oral_result['200']['species'])
    oral_genus_200 = len(oral_result['200']['genus'])

    oral_strain_500 = len(oral_result['500']['strain'])
    oral_species_500 = len(oral_result['500']['species'])
    oral_genus_500 = len(oral_result['500']['genus'])

    oral_strain_1000 = len(oral_result['1000']['strain'])
    oral_species_1000 = len(oral_result['1000']['species'])
    oral_genus_1000 = len(oral_result['1000']['genus'])

    print(oral_strain_200, oral_species_200, oral_genus_200)
    print(oral_strain_500, oral_species_500, oral_genus_500)
    print(oral_strain_1000, oral_species_1000, oral_genus_1000)

    line_width = 1

    plt.plot(['Genus', 'Species', 'Strain'],
             [skin_genus_1000, skin_species_1000,
              skin_strain_1000],
             label='1000', color='#ec7014', linewidth=line_width,
             marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [skin_genus_500, skin_species_500, skin_strain_500],
             label='500',
             color='#7570b3', linewidth=line_width, marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [skin_genus_200, skin_species_200, skin_strain_200],
             label='200',
             color='#1b9e77', linewidth=line_width, marker='o', )
    plt.legend()
    # plt.xticks(['Strain', 'Species', 'Genus'])
    plt.title('Skin', fontsize=15, alpha=1.0)
    plt.savefig('CAMI_II_skin_max_edges.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    plt.plot(['Genus', 'Species', 'Strain'],
             [oral_genus_1000, oral_species_1000, oral_strain_1000],
             label='1000', color='#ec7014', linewidth=line_width,
             marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [oral_genus_500, oral_species_500, oral_strain_500],
             label='500',
             color='#7570b3', linewidth=line_width, marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [oral_genus_200, oral_species_200, oral_strain_200],
             label='200',
             color='#1b9e77', linewidth=line_width, marker='o', )
    plt.legend()
    plt.title('Oral', fontsize=15, alpha=1.0)
    # plt.xticks(['Strain', 'Species', 'Genus'])
    plt.savefig('CAMI_II_oral_max_edges.pdf', dpi=300, bbox_inches='tight')
    plt.close()

def plot_recluster_alternative():
    def get_recluster(dataset='skin', input=None):
        method_list, species_list, genus_list = get_result(dataset, input)
        result = {'kmeans': {'strain': list(set(method_list['kmeans'])),
                         'species': list(set(species_list['kmeans'])),
                         'genus': list(set(genus_list['kmeans']))},
                  'spec': {'strain': list(set(method_list['spec'])),
                             'species': list(set(species_list['spec'])),
                             'genus': list(set(genus_list['spec']))},
                  'agg': {
                      'strain': list(set(method_list['agg'])),
                      'species': list(set(species_list['agg'])),
                      'genus': list(set(genus_list['agg']))},
                  'dbscan': {'strain': list(set(method_list['dbscan'])),
                              'species': list(set(species_list['dbscan'])),
                              'genus': list(set(genus_list['dbscan']))}, }
        return result

    skin_result = get_recluster('skin',
                              'updated_results/effect_clustering/recluster')
    skin_strain_kmeans = len(skin_result['kmeans']['strain'])
    skin_species_kmeans = len(skin_result['kmeans']['species'])
    skin_genus_kmeans = len(skin_result['kmeans']['genus'])

    skin_strain_spec = len(skin_result['spec']['strain'])
    skin_species_spec = len(skin_result['spec']['species'])
    skin_genus_spec = len(skin_result['spec']['genus'])

    skin_strain_agg = len(skin_result['agg']['strain'])
    skin_species_agg = len(skin_result['agg']['species'])
    skin_genus_agg = len(skin_result['agg']['genus'])

    skin_strain_dbscan = len(skin_result['dbscan']['strain'])
    skin_species_dbscan = len(skin_result['dbscan']['species'])
    skin_genus_dbscan = len(skin_result['dbscan']['genus'])

    print(skin_strain_kmeans, skin_species_kmeans, skin_genus_kmeans)
    print(skin_strain_spec, skin_species_spec, skin_genus_spec)
    print(skin_strain_agg, skin_species_agg,
          skin_genus_agg)
    print(skin_strain_dbscan, skin_species_dbscan, skin_genus_dbscan)

    oral_result = get_recluster('oral',
                              'updated_results/effect_clustering/recluster')
    oral_strain_kmeans = len(oral_result['kmeans']['strain'])
    oral_species_kmeans = len(oral_result['kmeans']['species'])
    oral_genus_kmeans = len(oral_result['kmeans']['genus'])

    oral_strain_spec = len(oral_result['spec']['strain'])
    oral_species_spec = len(oral_result['spec']['species'])
    oral_genus_spec = len(oral_result['spec']['genus'])

    oral_strain_agg = len(oral_result['agg']['strain'])
    oral_species_agg = len(oral_result['agg']['species'])
    oral_genus_agg = len(oral_result['agg']['genus'])

    oral_strain_dbscan = len(oral_result['dbscan']['strain'])
    oral_species_dbscan = len(oral_result['dbscan']['species'])
    oral_genus_dbscan = len(oral_result['dbscan']['genus'])

    print(oral_strain_kmeans, oral_species_kmeans, oral_genus_kmeans)
    print(oral_strain_spec, oral_species_spec, oral_genus_spec)
    print(oral_strain_agg, oral_species_agg,
          oral_genus_agg)
    print(oral_strain_dbscan, oral_species_dbscan, oral_genus_dbscan)

    line_width = 1

    plt.plot(['Genus', 'Species', 'Strain'],
             [skin_genus_dbscan, skin_species_dbscan,
              skin_strain_dbscan],
             label='DBSCAN', color='#ec7014', linewidth=line_width,
             marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [skin_genus_agg, skin_species_agg, skin_strain_agg],
             label='Agglomerative',
             color='#7570b3', linewidth=line_width, marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [skin_genus_spec, skin_species_spec, skin_strain_spec],
             label='Spectral',
             color='#e7298a', linewidth=line_width, marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [skin_genus_kmeans, skin_species_kmeans, skin_strain_kmeans],
             label='KMeans',
             color='#1b9e77', linewidth=line_width, marker='o', )
    plt.legend()
    # plt.xticks(['Strain', 'Species', 'Genus'])
    plt.title('Skin', fontsize=15, alpha=1.0)
    plt.savefig('CAMI_II_skin_recluster_alternative.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    plt.plot(['Genus', 'Species', 'Strain'],
             [oral_genus_dbscan, oral_species_dbscan, oral_strain_dbscan],
             label='DBSCAN', color='#ec7014', linewidth=line_width,
             marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [oral_genus_agg, oral_species_agg, oral_strain_agg],
             label='Agglomerative',
             color='#7570b3', linewidth=line_width, marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [oral_genus_spec, oral_species_spec, oral_strain_spec],
             label='Spectral',
             color='#e7298a', linewidth=line_width, marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [oral_genus_kmeans, oral_species_kmeans, oral_strain_kmeans],
             label='KMeans',
             color='#1b9e77', linewidth=line_width, marker='o', )
    plt.legend()
    plt.title('Oral', fontsize=15, alpha=1.0)
    # plt.xticks(['Strain', 'Species', 'Genus'])
    plt.savefig('CAMI_II_oral_recluster_alternative.pdf', dpi=300, bbox_inches='tight')
    plt.close()

def plot_embeddings():
    def get_recluster(dataset='skin', input=None):
        method_list, species_list, genus_list = get_result(dataset, input)
        result = {'hidden1': {'strain': list(set(method_list['hidden1'])),
                         'species': list(set(species_list['hidden1'])),
                         'genus': list(set(genus_list['hidden1']))},
                  'hidden2': {'strain': list(set(method_list['hidden2'])),
                             'species': list(set(species_list['hidden2'])),
                             'genus': list(set(genus_list['hidden2']))},
                  'output': {
                      'strain': list(set(method_list['output'])),
                      'species': list(set(species_list['output'])),
                      'genus': list(set(genus_list['output']))}, }
        return result

    skin_result = get_recluster('skin', 'updated_results/effect_embeddings')
    skin_strain_hidden1 = len(skin_result['hidden1']['strain'])
    skin_species_hidden1 = len(skin_result['hidden1']['species'])
    skin_genus_hidden1 = len(skin_result['hidden1']['genus'])

    skin_strain_hidden2 = len(skin_result['hidden2']['strain'])
    skin_species_hidden2 = len(skin_result['hidden2']['species'])
    skin_genus_hidden2 = len(skin_result['hidden2']['genus'])

    skin_strain_output = len(skin_result['output']['strain'])
    skin_species_output = len(skin_result['output']['species'])
    skin_genus_output = len(skin_result['output']['genus'])

    print(skin_strain_hidden1, skin_species_hidden1, skin_genus_hidden1)
    print(skin_strain_hidden2, skin_species_hidden2, skin_genus_hidden2)
    print(skin_strain_output, skin_species_output,
          skin_genus_output)

    oral_result = get_recluster('oral', 'updated_results/effect_embeddings')
    oral_strain_hidden1 = len(oral_result['hidden1']['strain'])
    oral_species_hidden1 = len(oral_result['hidden1']['species'])
    oral_genus_hidden1 = len(oral_result['hidden1']['genus'])

    oral_strain_hidden2 = len(oral_result['hidden2']['strain'])
    oral_species_hidden2 = len(oral_result['hidden2']['species'])
    oral_genus_hidden2 = len(oral_result['hidden2']['genus'])

    oral_strain_output = len(oral_result['output']['strain'])
    oral_species_output = len(oral_result['output']['species'])
    oral_genus_output = len(oral_result['output']['genus'])

    print(oral_strain_hidden1, oral_species_hidden1, oral_genus_hidden1)
    print(oral_strain_hidden2, oral_species_hidden2, oral_genus_hidden2)
    print(oral_strain_output, oral_species_output,
          oral_genus_output)

    line_width = 1

    plt.plot(['Genus', 'Species', 'Strain'],
             [skin_genus_hidden1, skin_species_hidden1,
              skin_strain_hidden1],
             label='Hidden1', color='#ec7014', linewidth=line_width,
             marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [skin_genus_hidden2, skin_species_hidden2, skin_strain_hidden2],
             label='Hidden2',
             color='#7570b3', linewidth=line_width, marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [skin_genus_output, skin_species_output, skin_strain_output],
             label='Output',
             color='#1b9e77', linewidth=line_width, marker='o', )
    plt.legend()
    # plt.xticks(['Strain', 'Species', 'Genus'])
    plt.title('Skin', fontsize=15, alpha=1.0)
    plt.savefig('CAMI_II_skin_embeddings.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    plt.plot(['Genus', 'Species', 'Strain'],
             [oral_genus_hidden1, oral_species_hidden1, oral_strain_hidden1],
             label='Hidden1', color='#ec7014', linewidth=line_width,
             marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [oral_genus_hidden2, oral_species_hidden2, oral_strain_hidden2],
             label='Hidden2',
             color='#7570b3', linewidth=line_width, marker='o', )
    plt.plot(['Genus', 'Species', 'Strain'],
             [oral_genus_output, oral_species_output, oral_strain_output],
             label='Output',
             color='#1b9e77', linewidth=line_width, marker='o', )
    plt.legend()
    plt.title('Oral', fontsize=15, alpha=1.0)
    # plt.xticks(['Strain', 'Species', 'Genus'])
    plt.savefig('CAMI_II_oral_embeddings.pdf', dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    ### bar plot in strain/species/genus level

    # plot_hq_num()
    #plot_comparsion_Metabat2()

    ### venn plot in strain/species/genus level
    #plot_venn_plot()

    ### plot similarity distribution of oral and skin
    #plot_similarity_distribution()

    ### bar plot with different similarities
    #plot_bar_strain_simiarity()

    plot_generalization()
    # plot_recluster()
    # plot_cluster_alternative()
    # plot_max_edge()
    # plot_recluster_alternative()
    # plot_embeddings()
