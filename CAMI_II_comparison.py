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
    genome_path = os.path.join('Results/Simulated/CAMI_II/{0}/amber_{0}_{1}'.format(dataset,data_index[0]), 'genome')
    method_list = {}
    species_list = {}
    genus_list = {}

    for root, dirs, files in os.walk(genome_path, topdown=False):
        for name in dirs:
            method_list[name] = []
            species_list[name] = []
            genus_list[name] = []

    for temp in data_index:
        taxi = pd.read_csv('Results/Simulated/CAMI_II/{0}/taxonomic_profile_{1}.txt'.format(dataset,temp),
                           sep='\t', skiprows=3, dtype={'@@TAXID': str, 'TAXPATH': str})
        taxi_genus = taxi[taxi['RANK'] == 'genus']['@@TAXID'].values.tolist()
        taxi_species = taxi[taxi['RANK'] == 'species']['@@TAXID'].values.tolist()
        method_path_list = []
        genome_path = os.path.join('Results/Simulated/CAMI_II/{0}/amber_{0}_{1}'.format(dataset,temp), 'genome')

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
                    else:
                        if taxi_split[-2] in taxi_genus:
                            genus_list[method_path.split('/')[-1]].append(taxi_split[-2])
                    if taxi_split[-3] in taxi_genus:
                        genus_list[method_path.split('/')[-1]].append(taxi_split[-3])

    result = {'Metabat2':{'strain':list(set(method_list['Metabat2'])),'species':list(set(species_list['Metabat2'])),'genus':list(set(genus_list['Metabat2']))}, 'VAMB': {'strain':list(set(method_list['VAMB'])),'species':list(set(species_list['VAMB'])),'genus':list(set(genus_list['VAMB']))}, 'SemiBin': {'strain':list(set(method_list['S3N2Bin'])),'species':list(set(species_list['S3N2Bin'])),'genus':list(set(genus_list['S3N2Bin']))}}

    return result

def get_generalization_hq_taxi(dataset = 'skin'):
    if dataset == 'skin':
        data_index = [1,13,14,15,16,17,18,19,20,28]
    else:
        data_index = [6,7,8,13,14,15,16,17,18,19]
    genome_path = os.path.join('Results/Simulated/CAMI_II/{0}/generalization/amber_{0}_{1}'.format(dataset,data_index[0]), 'genome')
    method_list = {}
    species_list = {}
    genus_list = {}

    for root, dirs, files in os.walk(genome_path, topdown=False):
        for name in dirs:
            method_list[name] = []
            species_list[name] = []
            genus_list[name] = []

    for temp in data_index:
        taxi = pd.read_csv('Results/Simulated/CAMI_II/{0}/taxonomic_profile_{1}.txt'.format(dataset,temp),
                           sep='\t', skiprows=3, dtype={'@@TAXID': str, 'TAXPATH': str})
        taxi_genus = taxi[taxi['RANK'] == 'genus']['@@TAXID'].values.tolist()
        taxi_species = taxi[taxi['RANK'] == 'species']['@@TAXID'].values.tolist()
        method_path_list = []
        genome_path = os.path.join('Results/Simulated/CAMI_II/{0}/generalization/amber_{0}_{1}'.format(dataset,temp), 'genome')
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
                    else:
                        if taxi_split[-2] in taxi_genus:
                            genus_list[method_path.split('/')[-1]].append(taxi_split[-2])
                    if taxi_split[-3] in taxi_genus:
                        genus_list[method_path.split('/')[-1]].append(taxi_split[-3])

    result = {'NoSemi':{'strain':list(set(method_list['NoSemi'])),'species':list(set(species_list['NoSemi'])),'genus':list(set(genus_list['NoSemi']))}, 'SemiBin_c': {'strain':list(set(method_list['S3N2Bin_c'])),'species':list(set(species_list['S3N2Bin_c'])),'genus':list(set(genus_list['S3N2Bin_c']))}, 'SemiBin_m': {'strain':list(set(method_list['S3N2Bin_m'])),'species':list(set(species_list['S3N2Bin_m'])),'genus':list(set(genus_list['S3N2Bin_m']))}, 'SemiBin_mc': {'strain':list(set(method_list['S3N2Bin_mc'])),'species':list(set(species_list['S3N2Bin_mc'])),'genus':list(set(genus_list['S3N2Bin_mc']))}}

    return result

def plot_comparsion_Metabat2():
    subset = pd.DataFrame(np.array([[66,74],[68,75]]),columns = ['Metabat2_single','Metabat2_multi'], index=['Skin','Oral'])
    ax = subset.plot(kind='bar',width = 0.8,color = ['#fdbf6f', '#ff7f00'],figsize=(3,4))
    ax.set_yticks(ticks=[0,20,40,60,80])
    ax.set_yticklabels(labels=[0,20,40,60,80],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality genomes', fontsize=15,color = 'black')
    ax.set_title('strain', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_com_strain_Metabat2.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    subset = pd.DataFrame(np.array([[63,65],[65,71]]),columns = ['Metabat2_single','Metabat2_multi'], index=['Skin','Oral'])
    ax = subset.plot(kind='bar',width= 0.8,legend = False,color = ['#fdbf6f', '#ff7f00'],figsize=(3,4))
    ax.set_yticks(ticks=[0,15,30,45,60,75])
    ax.set_yticklabels(labels=[0,15,30,45,60,75],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    #ax.set_ylabel('Num', fontsize=15,color = 'black')
    ax.set_title('species', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_com_species_Metabat2.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    subset = pd.DataFrame(np.array([[41, 43],[43, 47]]),columns = ['Metabat2_single','Metabat2_multi'], index=['Skin','Oral'])
    ax = subset.plot(kind='bar',width = 0.8,legend = False,color = ['#fdbf6f', '#ff7f00'],figsize=(3,4))
    ax.set_yticks(ticks=[0,10,20,30,40,50])
    ax.set_yticklabels(labels=[0,10,20,30,40,50],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    #ax.set_ylabel('Num', fontsize=15,color = 'black')
    ax.set_title('genus', fontsize=20, alpha=1.0,color = 'black')
    ax.legend(['Metabat2_single','Metabat2_multi'],fontsize=10)
    plt.savefig('CAMI_II_com_genus_Metabat2.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

def plot_bar_plot():
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

    subset = pd.DataFrame(np.array([[skin_strain_metabat2,skin_strain_vamb,skin_strain_SemiBin],[oral_strain_metabat2,oral_strain_vamb,oral_strain_SemiBin]]),columns = ['Metabat2','VAMB','SemiBin'], index=['Skin','Oral'])
    ax = subset.plot(kind='bar',width = 0.6,color = ['#fdbf6f', '#b2df8a', '#a6cee3'],figsize=(3,4))
    ax.set_yticks(ticks=[0,20,40,60,80,100,120,140])
    ax.set_yticklabels(labels=[0,20,40,60,80,100,120,140],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality genomes', fontsize=15,color = 'black')
    ax.set_title('strain', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_com_strain.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    subset = pd.DataFrame(np.array([[skin_species_metabat2,skin_species_vamb,skin_species_SemiBin],[oral_species_metabat2,oral_species_vamb,oral_species_SemiBin]]),columns = ['Metabat2','VAMB','SemiBin'], index=['Skin','Oral'])
    ax = subset.plot(kind='bar',width= 0.6,legend = False,color = ['#fdbf6f', '#b2df8a', '#a6cee3'],figsize=(3,4))
    ax.set_yticks(ticks=[0,15,30,45,60,75,90,105])
    ax.set_yticklabels(labels=[0,15,30,45,60,75,90,105],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    #ax.set_ylabel('Num', fontsize=15,color = 'black')
    ax.set_title('species', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_com_species.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    subset = pd.DataFrame(np.array([[skin_genus_metabat2,skin_genus_vamb,skin_genus_SemiBin],[oral_genus_metabat2,oral_genus_vamb,oral_genus_SemiBin]]),columns = ['Metabat2','VAMB','SemiBin'], index=['Skin','Oral'])
    ax = subset.plot(kind='bar',width = 0.6,legend = False,color = ['#fdbf6f', '#b2df8a', '#a6cee3'],figsize=(3,4))
    ax.set_yticks(ticks=[0,10,20,30,40,50,60])
    ax.set_yticklabels(labels=[0,10,20,30,40,50,60],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    #ax.set_ylabel('Num', fontsize=15,color = 'black')
    ax.set_title('genus', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_com_genus.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    plt.show()

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
    out = venn3_unweighted([set(skin_strain_metabat2),set(skin_strain_vamb),set(skin_strain_SemiBin)],set_labels=('Metabat2', 'VAMB','SemiBin'),normalize_to=5.0,set_colors=('#ff7f00', '#33a02c', '#1f78b4'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("strain(Skin)", fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_Strain_Skin.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(skin_species_metabat2),set(skin_species_vamb),set(skin_species_SemiBin)],set_labels=('Metabat2', 'VAMB','SemiBin'),normalize_to=5.0,set_colors=('#ff7f00', '#33a02c', '#1f78b4'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Skin)", fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_Species_Skin.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(skin_genus_metabat2),set(skin_genus_vamb),set(skin_genus_SemiBin)],set_labels=('Metabat2', 'VAMB','SemiBin'),normalize_to=5.0,set_colors=('#ff7f00', '#33a02c', '#1f78b4'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Skin)", fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_Genus_Skin.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(oral_strain_metabat2),set(oral_strain_vamb),set(oral_strain_SemiBin)],set_labels=('Metabat2', 'VAMB','SemiBin'),normalize_to=5.0,set_colors=('#ff7f00', '#33a02c', '#1f78b4'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("strain(Oral)", fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_Strain_Oral.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(oral_species_metabat2),set(oral_species_vamb),set(oral_species_SemiBin)],set_labels=('Metabat2', 'VAMB','SemiBin'),normalize_to=5.0,set_colors=('#ff7f00', '#33a02c', '#1f78b4'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Oral)", fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_Species_Oral.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(oral_genus_metabat2),set(oral_genus_vamb),set(oral_genus_SemiBin)],set_labels=('Metabat2', 'VAMB','SemiBin'),normalize_to=5.0,set_colors=('#ff7f00', '#33a02c', '#1f78b4'))
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

    plt.plot(x,y[1], label='2017.12.04_18.56.22_sample_1',color='#a6cee3')
    plt.plot(x, y[13], label='2017.12.04_18.56.22_sample_13',color='#b2df8a')
    plt.plot(x, y[14], label='2017.12.04_18.56.22_sample_14',color='#fb9a99')
    plt.plot(x, y[15], label='2017.12.04_18.56.22_sample_15',color='#fdbf6f')
    plt.plot(x, y[16], label='2017.12.04_18.56.22_sample_16',color='#1f78b4')
    plt.xticks(rotation=50,color = 'black')
    plt.yticks(color='black')
    plt.legend()
    plt.title('Genome distribution of Skin',fontsize = 15)

    plt.xlabel('ANI',fontsize = 15,color = 'black')
    plt.ylabel('Num',fontsize = 15,color = 'black')
    plt.savefig('CAMI_II_distribution_Skin_1.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()


    plt.plot(x,y[17], label='2017.12.04_18.56.22_sample_17',color='#a6cee3')
    plt.plot(x, y[18], label='2017.12.04_18.56.22_sample_18',color='#b2df8a')
    plt.plot(x, y[19], label='2017.12.04_18.56.22_sample_19',color='#fb9a99')
    plt.plot(x, y[20], label='2017.12.04_18.56.22_sample_20',color='#fdbf6f')
    plt.plot(x, y[28], label='2017.12.04_18.56.22_sample_28',color='#1f78b4')
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

    plt.plot(x,y[6], label='2017.12.04_18.56.22_sample_6',color='#a6cee3')
    plt.plot(x, y[7], label='2017.12.04_18.56.22_sample_7',color='#b2df8a')
    plt.plot(x, y[8], label='2017.12.04_18.56.22_sample_8',color='#fb9a99')
    plt.plot(x, y[13], label='2017.12.04_18.56.22_sample_13',color='#fdbf6f')
    plt.plot(x, y[14], label='2017.12.04_18.56.22_sample_14',color='#1f78b4')
    plt.xticks(rotation=50,color = 'black')
    plt.yticks(color='black')
    plt.legend()
    plt.title('Genome distribution of Oral',fontsize = 15)

    plt.xlabel('ANI',fontsize = 15,color = 'black')
    plt.ylabel('Num',fontsize = 15,color = 'black')
    plt.savefig('CAMI_II_distribution_Oral_1.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()


    plt.plot(x,y[15], label='2017.12.04_18.56.22_sample_15',color='#a6cee3')
    plt.plot(x, y[16], label='2017.12.04_18.56.22_sample_16',color='#b2df8a')
    plt.plot(x, y[17], label='2017.12.04_18.56.22_sample_17',color='#fb9a99')
    plt.plot(x, y[18], label='2017.12.04_18.56.22_sample_18',color='#fdbf6f')
    plt.plot(x, y[19], label='2017.12.04_18.56.22_sample_19',color='#1f78b4')
    plt.xticks(rotation=50,color = 'black')
    plt.yticks(color='black')
    plt.legend()
    plt.title('Genome distribution of Oral',fontsize = 15)

    plt.xlabel('ANI',fontsize = 15,color = 'black')
    plt.ylabel('Num',fontsize = 15,color = 'black')
    plt.savefig('CAMI_II_distribution_Oral_2.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

def get_hq_strain_similarity(dataset = 'skin'):
    if dataset == 'skin':
        data_index = [1,13,14,15,16,17,18,19,20,28]
    else:
        data_index = [6,7,8,13,14,15,16,17,18,19]

    distribution = get_similarity_distribution(dataset)
    genome_path = os.path.join('Results/Simulated/CAMI_II/{0}/amber_{0}_{1}'.format(dataset,data_index[0]), 'genome')
    method_list = {}

    for root, dirs, files in os.walk(genome_path, topdown=False):
        for name in dirs:
            method_list[os.path.join(root, name).split('/')[-1]] = {90: [], 95: [], 97: [], 98: [], 99: [],
                                                                          99.5: [], 99.9: [], 100: []}

    for temp in data_index:
        method_path_list = []
        genome_path = os.path.join('Results/Simulated/CAMI_II/{0}/amber_{0}_{1}'.format(dataset, temp), 'genome')
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
        value_array[i][0] = len(set(method_list['Metabat2'][value_list[i]]))
        value_array[i][1] = len(set(method_list['VAMB'][value_list[i]]))
        value_array[i][2] = len(set(method_list['S3N2Bin'][value_list[i]]))

    return value_array

def plot_bar_strain_simiarity():
    ### skin
    value_array = get_hq_strain_similarity()
    subset = pd.DataFrame(value_array,columns = ['Metabat2','VAMB','SemiBin'], index=['0~90','90~95','95~97','97~98','98~99','99~99.5','99.5~99.9','99.9~100'])
    ax = subset.plot(kind='bar',width = 0.6, color = ['#fdbf6f', '#b2df8a', '#a6cee3'])
    ax.set_yticks(ticks=[0,20,40,60,80])
    ax.set_yticklabels(labels=[0,20,40,60,80],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['0~90','90~95','95~97','97~98','98~99','99~99.5','99.5~99.9','99.9~100'], fontsize=15,color = 'black',rotation = 50)
    ax.set_ylabel('High quality genomes', fontsize=15,color = 'black')
    ax.set_title('Skin', fontsize=15, alpha=1.0)
    plt.savefig('CAMI_II_distribution_Skin_bar.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    ### Oral
    value_array = get_hq_strain_similarity(dataset= 'oral')

    subset = pd.DataFrame(value_array,columns = ['Metabat2','VAMB','SemiBin'], index=['0~90','90~95','95~97','97~98','98~99','99~99.5','99.5~99.9','99.9~100'])
    ax = subset.plot(kind='bar',width = 0.6, color = ['#fdbf6f', '#b2df8a', '#a6cee3'])
    ax.set_yticks(ticks=[0,30,60,90])
    ax.set_yticklabels(labels=[0,30,60,90],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['0~90','90~95','95~97','97~98','98~99','99~99.5','99.5~99.9','99.9~100'], fontsize=15,color = 'black',rotation = 50)
    ax.set_ylabel('High quality genomes', fontsize=15,color = 'black')
    ax.set_title('Oral', fontsize=15, alpha=1.0)
    plt.savefig('CAMI_II_distribution_Oral_bar.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()


def plot_abundance(dataset = 'skin',level = 'strain', output = None):
    """
    dataset: skin or oral
    level: strain or species
    """
    if dataset == 'skin':
        data_index = [1,13,14,15,16,17,18,19,20,28]
    else:
        data_index = [6,7,8,13,14,15,16,17,18,19]

    environment_dict = {}

    for index in data_index:
        taxi = pd.read_csv('Results/Simulated/CAMI_II/{0}/taxonomic_profile_{1}.txt'.format(dataset,index),
                           sep='\t', skiprows=3, dtype={'@@TAXID': str, 'TAXPATH': str})
        taxi_species = taxi[taxi['RANK'] == 'species']['@@TAXID'].values.tolist()
        abundance_dict = {}
        abundance = pd.read_csv('Results/Simulated/CAMI_II/{0}/abundance{1}.tsv'.format(dataset, index),sep='\t',header=None).values
        for temp in abundance:
            abundance_dict[temp[0]] = float(temp[1])
        metric = pd.read_csv(os.path.join('Results/Simulated/CAMI_II/{0}/amber_{0}_{1}/genome/Gold standard'.format(dataset, index), 'metrics_per_bin.tsv'), sep='\t')
        strain_list = metric['Most abundant genome'].values.tolist()
        if level == 'strain':
            for temp_strain in strain_list:
                if temp_strain in taxi['_CAMI_GENOMEID'].values.tolist():
                    if temp_strain not in environment_dict:
                        environment_dict[temp_strain] = abundance_dict[temp_strain]
                    else:
                        environment_dict[temp_strain] += abundance_dict[temp_strain]
        else:
            for temp_strain in strain_list:
                if temp_strain in taxi['_CAMI_GENOMEID'].values.tolist():
                    taxi_split = taxi[taxi['_CAMI_GENOMEID'] == temp_strain]['TAXPATH'].values[0].split('|')
                    if taxi_split[-2] in taxi_species:
                        if taxi_split[-2] not in environment_dict:
                            environment_dict[taxi_split[-2]] = abundance_dict[temp_strain]
                        else:
                            environment_dict[taxi_split[-2]] += abundance_dict[temp_strain]

    result = get_hq_taxi(dataset)
    if level == 'strain':
        strain_metabat2 = result['Metabat2']['strain']
        strain_vamb = result['VAMB']['strain']
        strain_SemiBin = result['SemiBin']['strain']
        abun_metabat2 = [environment_dict[strain] for strain in strain_metabat2]
        abun_vamb = [environment_dict[strain] for strain in strain_vamb]
        abun_SemiBin = [environment_dict[strain] for strain in strain_SemiBin]

    if level == 'species':
        species_metabat2 = result['Metabat2']['species']
        species_vamb = result['VAMB']['species']
        species_SemiBin = result['SemiBin']['species']
        abun_metabat2 = [environment_dict[species] for species in species_metabat2]
        abun_vamb = [environment_dict[species] for species in species_vamb]
        abun_SemiBin = [environment_dict[species] for species in species_SemiBin]

    SemiBin_abun = {'S3N2Bin': abun_SemiBin}
    SemiBin_abun = pd.DataFrame(SemiBin_abun)

    Metabat2_abun = {'Metabat2': abun_metabat2}
    Metabat2_abun = pd.DataFrame(Metabat2_abun)

    vamb_abun = {'VAMB': abun_vamb}
    vamb_abun = pd.DataFrame(vamb_abun)

    subset = pd.concat([SemiBin_abun, Metabat2_abun], axis=1)
    subset = pd.concat([subset, vamb_abun], axis=1)
    subset = subset.apply(np.log10)

    fig, ax = plt.subplots(nrows=1, ncols=1)
    subset.columns = [['SemiBin', 'Metabat2', 'VAMB']]
    column_new = ['Metabat2', 'VAMB', 'SemiBin', ]
    subset = subset[column_new]

    sns.violinplot(data=subset, width=.5, fliersize=0, scale="count", palette= ['#fdbf6f', '#b2df8a', '#a6cee3'])
    ax.set_yticks(ticks=[0, 1, 2, 3])
    ax.set_yticklabels(labels=[0, 1, 2, 3], fontsize=15, color='black')
    ax.set_xticklabels(
        labels=['Metabat2', 'VAMB', 'SemiBin',], fontsize=15, color='black')
    ax.set_ylabel('Abundance', fontsize=15, color='black')
    ax.set_title('{}'.format('{0}({1})'.format(dataset,level)), fontsize=15, alpha=1.0)

    plt.tight_layout()
    plt.savefig(output , dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

def plot_generalization():
    result = get_hq_taxi('skin')
    skin_SemiBin_strain = len(result['SemiBin']['strain'])
    skin_SemiBin_species = len(result['SemiBin']['species'])
    skin_SemiBin_genus = len(result['SemiBin']['genus'])

    result = get_generalization_hq_taxi('skin')
    skin_NoSemi_strain = len(result['NoSemi']['strain'])
    skin_NoSemi_species = len(result['NoSemi']['species'])
    skin_NoSemi_genus = len(result['NoSemi']['genus'])

    skin_SemiBin_m_strain = len(result['SemiBin_c']['strain'])
    skin_SemiBin_m_species = len(result['SemiBin_c']['species'])
    skin_SemiBin_m_genus = len(result['SemiBin_c']['genus'])

    skin_SemiBin_c_strain = len(result['SemiBin_m']['strain'])
    skin_SemiBin_c_species = len(result['SemiBin_m']['species'])
    skin_SemiBin_c_genus = len(result['SemiBin_m']['genus'])

    skin_SemiBin_mc_strain = len(result['SemiBin_mc']['strain'])
    skin_SemiBin_mc_species = len(result['SemiBin_mc']['species'])
    skin_SemiBin_mc_genus = len(result['SemiBin_mc']['genus'])

    result = get_hq_taxi('oral')
    oral_SemiBin_strain = len(result['SemiBin']['strain'])
    oral_SemiBin_species = len(result['SemiBin']['species'])
    oral_SemiBin_genus = len(result['SemiBin']['genus'])

    result = get_generalization_hq_taxi('oral')
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

    subset = pd.DataFrame(np.array([[skin_NoSemi_strain,skin_SemiBin_m_strain,skin_SemiBin_c_strain,skin_SemiBin_mc_strain,skin_SemiBin_strain ],[oral_NoSemi_strain,oral_SemiBin_m_strain,oral_SemiBin_c_strain,oral_SemiBin_mc_strain,oral_SemiBin_strain]]),columns = ['No_semi','SemiBin_m','SemiBin_c','SemiBin_mc','SemiBin'], index=['Skin','Oral'])
    print(subset)
    ax = subset.plot(kind='bar',width = 0.6,color = ['#fdbf6f','#fb9a99', '#b2df8a', '#a6cee3', '#1f78b4'])
    ax.set_yticks(ticks=[0,20,40,60,80,100,120,140])
    ax.set_yticklabels(labels=[0,20,40,60,80,100,120,140],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality strains', fontsize=15,color = 'black')
    ax.set_title('strain', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_com_strain_generalization.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    subset = pd.DataFrame(np.array([[skin_NoSemi_species,skin_SemiBin_m_species,skin_SemiBin_c_species,skin_SemiBin_mc_species,skin_SemiBin_species ],[oral_NoSemi_species,oral_SemiBin_m_species,oral_SemiBin_c_species,oral_SemiBin_mc_species,oral_SemiBin_species]]),columns = ['No_semi','SemiBin_m','SemiBin_c','SemiBin_mc','SemiBin'], index=['Skin','Oral'])
    print(subset)
    ax = subset.plot(kind='bar',width = 0.6,color = ['#fdbf6f','#fb9a99', '#b2df8a', '#a6cee3', '#1f78b4'])
    ax.set_yticks(ticks=[0,15,30,45,60,75,90,105])
    ax.set_yticklabels(labels=[0,15,30,45,60,75,90,105],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality species', fontsize=15,color = 'black')
    ax.set_title('species', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_com_species_generalization.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    subset = pd.DataFrame(np.array([[skin_NoSemi_genus,skin_SemiBin_m_genus,skin_SemiBin_c_genus,skin_SemiBin_mc_genus,skin_SemiBin_genus ],[oral_NoSemi_genus,oral_SemiBin_m_genus,oral_SemiBin_c_genus,oral_SemiBin_mc_genus,oral_SemiBin_genus]]),columns = ['No_semi','SemiBin_m','SemiBin_c','SemiBin_mc','SemiBin'], index=['Skin','Oral'])
    print(subset)
    ax = subset.plot(kind='bar',width = 0.6,color = ['#fdbf6f','#fb9a99', '#b2df8a', '#a6cee3', '#1f78b4'])
    ax.set_yticks(ticks=[0,10,20,30,40,50,60])
    ax.set_yticklabels(labels=[0,10,20,30,40,50,60],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality genus', fontsize=15,color = 'black')
    ax.set_title('genus', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('CAMI_II_com_genus_generalization.pdf', dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    ### bar plot in strain/species/genus level
    plot_bar_plot()
    plot_comparsion_Metabat2()

    ### venn plot in strain/species/genus level
    plot_venn_plot()

    ### plot similarity distribution of oral and skin
    plot_similarity_distribution()

    ### bar plot with different similarities
    plot_bar_strain_simiarity()

    # ### violinplot of the abundance
    plot_abundance(dataset='skin',level = 'species', output='CAMI_II_abundance_skin_species.pdf')
    plot_abundance(dataset='skin', level='strain', output='CAMI_II_abundance_skin_strain.pdf')
    plot_abundance(dataset='oral', level='species', output='CAMI_II_abundance_oral_species.pdf')
    plot_abundance(dataset='oral', level='strain', output='CAMI_II_abundance_oral_strain.pdf')

    plot_generalization()