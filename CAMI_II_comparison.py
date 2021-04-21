"""
This script is used to reproduce the plot of the CAMI II results(compared to other binners).
"""
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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

    result = {'Metabat2':{'strain':list(set(method_list['Metabat2'])),'species':list(set(species_list['Metabat2'])),'genus':list(set(genus_list['Metabat2']))}, 'VAMB': {'strain':list(set(method_list['VAMB'])),'species':list(set(species_list['VAMB'])),'genus':list(set(genus_list['VAMB']))}, 'S3N2Bin': {'strain':list(set(method_list['S3N2Bin'])),'species':list(set(species_list['S3N2Bin'])),'genus':list(set(genus_list['S3N2Bin']))}}

    return result

def plot_bar_plot():
    skin_result = get_hq_taxi()
    skin_strain_metabat2 = len(skin_result['Metabat2']['strain'])
    skin_species_metabat2 = len(skin_result['Metabat2']['species'])
    skin_genus_metabat2 = len(skin_result['Metabat2']['genus'])

    skin_strain_vamb= len(skin_result['VAMB']['strain'])
    skin_species_vamb = len(skin_result['VAMB']['species'])
    skin_genus_vamb = len(skin_result['VAMB']['genus'])

    skin_strain_s3n2bin = len(skin_result['S3N2Bin']['strain'])
    skin_species_s3n2bin = len(skin_result['S3N2Bin']['species'])
    skin_genus_s3n2bin = len(skin_result['S3N2Bin']['genus'])

    oral_result = get_hq_taxi(dataset='oral')
    oral_strain_metabat2 = len(oral_result['Metabat2']['strain'])
    oral_species_metabat2 = len(oral_result['Metabat2']['species'])
    oral_genus_metabat2 = len(oral_result['Metabat2']['genus'])

    oral_strain_vamb = len(oral_result['VAMB']['strain'])
    oral_species_vamb = len(oral_result['VAMB']['species'])
    oral_genus_vamb = len(oral_result['VAMB']['genus'])

    oral_strain_s3n2bin = len(oral_result['S3N2Bin']['strain'])
    oral_species_s3n2bin = len(oral_result['S3N2Bin']['species'])
    oral_genus_s3n2bin = len(oral_result['S3N2Bin']['genus'])

    subset = pd.DataFrame(np.array([[skin_strain_metabat2,skin_strain_vamb,skin_strain_s3n2bin],[oral_strain_metabat2,oral_strain_vamb,oral_strain_s3n2bin]]),columns = ['Metabat2','VAMB','$S^3N^2Bin$'], index=['Skin','Oral'])
    ax = subset.plot(kind='bar',width = 0.6)
    ax.set_yticks(ticks=[0,20,40,60,80,100,120,140])
    ax.set_yticklabels(labels=[0,20,40,60,80,100,120,140],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality genomes', fontsize=15,color = 'black')
    ax.set_title('strain', fontsize=20, alpha=1.0,color = 'black')

    subset = pd.DataFrame(np.array([[skin_species_metabat2,skin_species_vamb,skin_species_s3n2bin],[oral_species_metabat2,oral_species_vamb,oral_species_s3n2bin]]),columns = ['Metabat2','VAMB','$S^3N^2Bin$'], index=['Skin','Oral'])
    ax = subset.plot(kind='bar',width= 0.6,legend = False)
    ax.set_yticks(ticks=[0,15,30,45,60,75,90,105])
    ax.set_yticklabels(labels=[0,15,30,45,60,75,90,105],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('Num', fontsize=15,color = 'black')
    ax.set_title('species', fontsize=20, alpha=1.0,color = 'black')

    subset = pd.DataFrame(np.array([[skin_genus_metabat2,skin_genus_vamb,skin_genus_s3n2bin],[oral_genus_metabat2,oral_genus_vamb,oral_genus_s3n2bin]]),columns = ['Metabat2','VAMB','$S^3N^2Bin$'], index=['Skin','Oral'])
    ax = subset.plot(kind='bar',width = 0.6,legend = False)
    ax.set_yticks(ticks=[0,10,20,30,40,50,60])
    ax.set_yticklabels(labels=[0,10,20,30,40,50,60],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('Num', fontsize=15,color = 'black')
    ax.set_title('genus', fontsize=20, alpha=1.0,color = 'black')

    plt.show()

def plot_venn_plot():
    skin_result = get_hq_taxi()
    skin_strain_metabat2 = skin_result['Metabat2']['strain']
    skin_species_metabat2 = skin_result['Metabat2']['species']
    skin_genus_metabat2 = skin_result['Metabat2']['genus']

    skin_strain_vamb= skin_result['VAMB']['strain']
    skin_species_vamb = skin_result['VAMB']['species']
    skin_genus_vamb = skin_result['VAMB']['genus']

    skin_strain_s3n2bin = skin_result['S3N2Bin']['strain']
    skin_species_s3n2bin = skin_result['S3N2Bin']['species']
    skin_genus_s3n2bin = skin_result['S3N2Bin']['genus']

    oral_result = get_hq_taxi(dataset='oral')
    oral_strain_metabat2 = oral_result['Metabat2']['strain']
    oral_species_metabat2 = oral_result['Metabat2']['species']
    oral_genus_metabat2 = oral_result['Metabat2']['genus']

    oral_strain_vamb = oral_result['VAMB']['strain']
    oral_species_vamb = oral_result['VAMB']['species']
    oral_genus_vamb = oral_result['VAMB']['genus']

    oral_strain_s3n2bin = oral_result['S3N2Bin']['strain']
    oral_species_s3n2bin = oral_result['S3N2Bin']['species']
    oral_genus_s3n2bin = oral_result['S3N2Bin']['genus']

    from matplotlib_venn import venn3_unweighted
    out = venn3_unweighted([set(skin_strain_metabat2),set(skin_strain_vamb),set(skin_strain_s3n2bin)],set_labels=('Metabat2', 'VAMB','$S^3N^2Bin$'),normalize_to=5.0)
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("strain(Skin)", fontsize=20, alpha=1.0,color = 'black')
    plt.show()

    out = venn3_unweighted([set(skin_species_metabat2),set(skin_species_vamb),set(skin_species_s3n2bin)],set_labels=('Metabat2', 'VAMB','$S^3N^2Bin$'),normalize_to=5.0)
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Skin)", fontsize=20, alpha=1.0,color = 'black')
    plt.show()

    out = venn3_unweighted([set(skin_genus_metabat2),set(skin_genus_vamb),set(skin_genus_s3n2bin)],set_labels=('Metabat2', 'VAMB','$S^3N^2Bin$'),normalize_to=5.0)
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Skin)", fontsize=20, alpha=1.0,color = 'black')
    plt.show()

    out = venn3_unweighted([set(oral_strain_metabat2),set(oral_strain_vamb),set(oral_strain_s3n2bin)],set_labels=('Metabat2', 'VAMB','$S^3N^2Bin$'),normalize_to=5.0)
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("strain(Oral)", fontsize=20, alpha=1.0,color = 'black')
    plt.show()

    out = venn3_unweighted([set(oral_species_metabat2),set(oral_species_vamb),set(oral_species_s3n2bin)],set_labels=('Metabat2', 'VAMB','$S^3N^2Bin$'),normalize_to=5.0)
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Oral)", fontsize=20, alpha=1.0,color = 'black')
    plt.show()

    out = venn3_unweighted([set(oral_genus_metabat2),set(oral_genus_vamb),set(oral_genus_s3n2bin)],set_labels=('Metabat2', 'VAMB','$S^3N^2Bin$'),normalize_to=5.0)
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Oral)", fontsize=20, alpha=1.0,color = 'black')
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

    plt.plot(x,y[1], label='2017.12.04_18.56.22_sample_1')
    plt.plot(x, y[13], label='2017.12.04_18.56.22_sample_13')
    plt.plot(x, y[14], label='2017.12.04_18.56.22_sample_14')
    plt.plot(x, y[15], label='2017.12.04_18.56.22_sample_15')
    plt.plot(x, y[16], label='2017.12.04_18.56.22_sample_16')
    plt.xticks(rotation=50,color = 'black')
    plt.yticks(color='black')
    plt.legend()
    plt.title('Genome distribution of Skin',fontsize = 15)

    plt.xlabel('ANI',fontsize = 15,color = 'black')
    plt.ylabel('Num',fontsize = 15,color = 'black')
    plt.show()


    plt.plot(x,y[17], label='2017.12.04_18.56.22_sample_17')
    plt.plot(x, y[18], label='2017.12.04_18.56.22_sample_18')
    plt.plot(x, y[19], label='2017.12.04_18.56.22_sample_19')
    plt.plot(x, y[20], label='2017.12.04_18.56.22_sample_20')
    plt.plot(x, y[28], label='2017.12.04_18.56.22_sample_28')
    plt.xticks(rotation=50,color = 'black')
    plt.yticks(color='black')
    plt.legend()
    plt.title('Genome distribution of Skin',fontsize = 15)

    plt.xlabel('ANI',fontsize = 15,color = 'black')
    plt.ylabel('Num',fontsize = 15,color = 'black')
    plt.show()

    oral_distribution = get_similarity_distribution(dataset='oral')

    x = ['0~90','90~95','95~97','97~98','98~99','99~99.5','99.5~99.9','99.9~100']
    y = {}
    for sample in  [6,7,8,13,14,15,16,17,18,19]:
        y[sample] = []
        for value in oral_distribution[sample]:
            y[sample].append(len(oral_distribution[sample][value]))
    print(y)
    plt.plot(x,y[6], label='2017.12.04_18.56.22_sample_6')
    plt.plot(x, y[7], label='2017.12.04_18.56.22_sample_7')
    plt.plot(x, y[8], label='2017.12.04_18.56.22_sample_8')
    plt.plot(x, y[13], label='2017.12.04_18.56.22_sample_13')
    plt.plot(x, y[14], label='2017.12.04_18.56.22_sample_14')
    plt.xticks(rotation=50,color = 'black')
    plt.yticks(color='black')
    plt.legend()
    plt.title('Genome distribution of Oral',fontsize = 15)

    plt.xlabel('ANI',fontsize = 15,color = 'black')
    plt.ylabel('Num',fontsize = 15,color = 'black')
    plt.show()


    plt.plot(x,y[15], label='2017.12.04_18.56.22_sample_15')
    plt.plot(x, y[16], label='2017.12.04_18.56.22_sample_16')
    plt.plot(x, y[17], label='2017.12.04_18.56.22_sample_17')
    plt.plot(x, y[18], label='2017.12.04_18.56.22_sample_18')
    plt.plot(x, y[19], label='2017.12.04_18.56.22_sample_19')
    plt.xticks(rotation=50,color = 'black')
    plt.yticks(color='black')
    plt.legend()
    plt.title('Genome distribution of Oral',fontsize = 15)

    plt.xlabel('ANI',fontsize = 15,color = 'black')
    plt.ylabel('Num',fontsize = 15,color = 'black')
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

    return  value_array

def plot_bar_strain_simiarity():
    ### skin
    value_array = get_hq_strain_similarity()
    subset = pd.DataFrame(value_array,columns = ['Metabat2','VAMB','$S^3N^2Bin$'], index=['0~90','90~95','95~97','97~98','98~99','99~99.5','99.5~99.9','99.9~100'])
    ax = subset.plot(kind='bar',width = 0.6)
    ax.set_yticks(ticks=[0,20,40,60,80])
    ax.set_yticklabels(labels=[0,20,40,60,80],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['0~90','90~95','95~97','97~98','98~99','99~99.5','99.5~99.9','99.9~100'], fontsize=15,color = 'black',rotation = 50)
    ax.set_ylabel('High quality genomes', fontsize=15,color = 'black')
    ax.set_title('Skin', fontsize=15, alpha=1.0)
    plt.show()

    ### Oral
    value_array = get_hq_strain_similarity(dataset= 'oral')

    subset = pd.DataFrame(value_array,columns = ['Metabat2','VAMB','$S^3N^2Bin$'], index=['0~90','90~95','95~97','97~98','98~99','99~99.5','99.5~99.9','99.9~100'])
    ax = subset.plot(kind='bar',width = 0.6)
    ax.set_yticks(ticks=[0,30,60,90])
    ax.set_yticklabels(labels=[0,30,60,90],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['0~90','90~95','95~97','97~98','98~99','99~99.5','99.5~99.9','99.9~100'], fontsize=15,color = 'black',rotation = 50)
    ax.set_ylabel('High quality genomes', fontsize=15,color = 'black')
    ax.set_title('Oral', fontsize=15, alpha=1.0)
    plt.show()


def plot_abundance(dataset = 'skin',level = 'strain'):
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
        strain_s3n2bin = result['S3N2Bin']['strain']
        abun_metabat2 = [environment_dict[strain] for strain in strain_metabat2]
        abun_vamb = [environment_dict[strain] for strain in strain_vamb]
        abun_s3n2bin = [environment_dict[strain] for strain in strain_s3n2bin]

    if level == 'species':
        species_metabat2 = result['Metabat2']['species']
        species_vamb = result['VAMB']['species']
        species_s3n2bin = result['S3N2Bin']['species']
        abun_metabat2 = [environment_dict[species] for species in species_metabat2]
        abun_vamb = [environment_dict[species] for species in species_vamb]
        abun_s3n2bin = [environment_dict[species] for species in species_s3n2bin]

    S3N2Bin_abun = {'S3N2Bin': abun_s3n2bin}
    S3N2Bin_abun = pd.DataFrame(S3N2Bin_abun)

    Metabat2_abun = {'Metabat2': abun_metabat2}
    Metabat2_abun = pd.DataFrame(Metabat2_abun)

    vamb_abun = {'VAMB': abun_vamb}
    vamb_abun = pd.DataFrame(vamb_abun)

    subset = pd.concat([S3N2Bin_abun, Metabat2_abun], axis=1)
    subset = pd.concat([subset, vamb_abun], axis=1)
    subset = subset.apply(np.log10)

    fig, ax = plt.subplots(nrows=1, ncols=1)
    subset.columns = [['S3N2Bin', 'Metabat2', 'VAMB']]
    column_new = ['Metabat2', 'VAMB', 'S3N2Bin', ]
    subset = subset[column_new]

    sns.violinplot(data=subset, width=.5, fliersize=0, scale="count")
    ax.set_yticks(ticks=[0, 1, 2, 3])
    ax.set_yticklabels(labels=[0, 1, 2, 3], fontsize=15, color='black')
    ax.set_xticklabels(
        labels=['Metabat2', 'VAMB', '$S^3N^2Bin$',], fontsize=15, color='black')
    ax.set_ylabel('Abundance', fontsize=15, color='black')
    ax.set_title('{}'.format('{0}({1})'.format(dataset,level)), fontsize=15, alpha=1.0)

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    ### bar plot in strain/species/genus level
    plot_bar_plot()

    ### venn plot in strain/species/genus level
    plot_venn_plot()

    ### plot similarity distribution of oral and skin
    plot_similarity_distribution()

    ### bar plot with different similarities
    plot_bar_strain_simiarity()

    ### violinplot of the abundance
    plot_abundance(dataset='skin',level = 'species')
    plot_abundance(dataset='skin', level='strain')
    plot_abundance(dataset='oral', level='species')
    plot_abundance(dataset='oral', level='strain')