

def plot_per_sample_comparison(dataset = 'dog', output = None):
    """
    dataset: dog, human, tara
    """
    SemiBin_result = get_result(dataset, method='SemiBin_pretrain')
    Metabat2_result = get_result(dataset, method='Metabat2')

    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    colormap = ['#084594','#2171b5','#4292c6','#6baed6','#9ecae1','#c6dbef','#deebf7']
    newcmp = LinearSegmentedColormap.from_list('cmps', colormap)
    result = {}

    num_SemiBin = []
    num_Metabat2 = []

    for sample in SemiBin_result:
        num_SemiBin.append(len(SemiBin_result[sample]['high quality']))
        num_Metabat2.append(len(Metabat2_result[sample]['high quality']))
    print(wilcoxon(num_Metabat2, num_SemiBin))
    for sample in SemiBin_result:
        if (len(SemiBin_result[sample]['high quality']), len(Metabat2_result[sample]['high quality'])) not in result:
            result[(len(SemiBin_result[sample]['high quality']), len(Metabat2_result[sample]['high quality']))] = 1
        else:
            result[(len(SemiBin_result[sample]['high quality']), len(Metabat2_result[sample]['high quality']))] += 1

    if dataset != 'tara':
        data = np.zeros(shape=(len(result), 3))

        for i, temp in enumerate(result):
            data[i][0] = temp[0]
            data[i][1] = temp[1]
            data[i][2] = result[temp]

        data = pd.DataFrame(data).astype(int)

        data.columns = ['x', 'y', 'num']

        plt.scatter(data.x, data.y,
                    c=data.num, s=(data.num ** 2) * 60, cmap=newcmp)

        plt.colorbar(shrink=0.5)
        plt.xlabel("SemiBin")
        plt.ylabel("Metabat2")
        if dataset == 'dog':
            plt.plot([0, 35], [0, 35], c='black')
            plt.title('Dog gut'.format(dataset), fontsize=15)
        else:
            plt.plot([0, 40], [0, 40], c='black')
            plt.title('Human gut'.format(dataset), fontsize=15)
        plt.show()
        plt.savefig(output, dpi=300, bbox_inches='tight')
        plt.close()


    else:
        data = []
        for i, temp in enumerate(result):
            if temp[0] == 0 and temp[1] == 0:
                continue
            data.append([temp[0], temp[1], result[temp]])
        data = np.array(data)
        data = pd.DataFrame(data).astype(int)

        data.columns = ['x', 'y', 'num']

        plt.scatter(data.x, data.y,
                    c=data.num, s=(data.num) * 60, cmap=newcmp)

        plt.colorbar()
        plt.xlabel("SemiBin")
        plt.ylabel("Metabat2")
        plt.plot([0, 25], [0, 25], c='black')
        plt.title('Ocean', fontsize=15)
        plt.show()
        plt.savefig(output, dpi=300, bbox_inches='tight')
        plt.close()


def plot_high_quality_comparison():
    num_dog_maxbin2_single = get_num_high_quality()
    num_dog_vamb_single = get_num_high_quality(method='VAMB')
    num_dog_metabat2_single = get_num_high_quality(method='Metabat2')
    num_dog_semibin_single = get_num_high_quality(method='S3N2Bin')
    num_dog_semibin_pretrain_single = get_num_high_quality(method='SemiBin_pretrain')

    num_human_maxbin2_single = get_num_high_quality(dataset='human')
    num_human_vamb_single = get_num_high_quality(dataset='human', method='VAMB')
    num_human_metabat2_single = get_num_high_quality(dataset='human', method='Metabat2')
    num_human_semibin_single = get_num_high_quality(dataset='human', method='S3N2Bin')
    num_human_semibin_pretrain_single = get_num_high_quality(dataset='human', method='SemiBin_pretrain')

    num_tara_maxbin2_single = get_num_high_quality(dataset='tara')
    num_tara_vamb_single = get_num_high_quality(dataset='tara', method='VAMB')
    num_tara_metabat2_single = get_num_high_quality(dataset='tara', method='Metabat2')
    num_tara_semibin_single = get_num_high_quality(dataset='tara', method='S3N2Bin')
    num_tara_semibin_pretrain_single = get_num_high_quality(dataset='tara', method='SemiBin_pretrain')

    num_dog_vamb_mulit = get_num_high_quality(method='VAMB',binning_mode='multi_sample')
    num_dog_semibin_multi = get_num_high_quality(method='S3N2Bin',binning_mode='multi_sample')

    num_human_vamb_multi = get_num_high_quality(dataset='human', method='VAMB',binning_mode='multi_sample')
    num_human_semibin_multi = get_num_high_quality(dataset='human', method='S3N2Bin',binning_mode='multi_sample')

    num_tara_vamb_multi = get_num_high_quality(dataset='tara', method='VAMB', binning_mode='multi_sample')
    num_tara_semibin_multi = get_num_high_quality(dataset='tara', method='S3N2Bin', binning_mode='multi_sample')

    # print(num_dog_maxbin2_single, num_dog_vamb_single, num_dog_metabat2_single, num_dog_semibin_single,
    #       num_dog_semibin_pretrain_single)
    # print(num_human_maxbin2_single, num_human_vamb_single, num_human_metabat2_single, num_human_semibin_single,
    #       num_human_semibin_pretrain_single)
    # print(num_tara_maxbin2_single, num_tara_vamb_single, num_tara_metabat2_single, num_tara_semibin_single,
    #       num_tara_semibin_pretrain_single)
    #
    # print(num_dog_vamb_mulit,num_dog_semibin_multi)
    # print(num_human_vamb_multi, num_human_semibin_multi)
    # print(num_tara_vamb_multi,num_tara_semibin_multi)

    subset = pd.DataFrame(np.array([[num_dog_maxbin2_single,num_dog_vamb_single,num_dog_metabat2_single,num_dog_semibin_single, num_dog_semibin_pretrain_single,num_dog_vamb_mulit,num_dog_semibin_multi]]),columns = ['Maxbin2','VAMB(single)','Metabat2','SemiBin(single)','SemiBin(pre-train)','VAMB(multi)', 'SemiBin(multi)'], index=['Dog gut'])
    ax = subset.plot(kind='bar',figsize=(2,4),legend = False, color=['#fb9a99','#b2df8a','#fdbf6f','#a6cee3','#1d91c0','#33a02c','#1f78b4'])
    ax.set_yticks(ticks=[0,500,1000,1500,2000,2500,3000,3500])
    ax.set_yticklabels(labels=[0,500,1000,1500,2000,2500,3000,3500],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Dog gut'], fontsize=15,color = 'black',rotation = 360)
    plt.savefig('Real_dog_hq.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()
    print(subset)
    subset = pd.DataFrame(np.array([[num_human_maxbin2_single,num_human_vamb_single,num_human_metabat2_single,num_human_semibin_single, num_human_semibin_pretrain_single,num_human_vamb_multi,num_human_semibin_multi]]),columns = ['Maxbin2','VAMB(single)','Metabat2','SemiBin(single)','SemiBin(pre-train)','VAMB(multi)', 'SemiBin(multi)'], index=['Human gut'])
    ax = subset.plot(kind='bar',figsize=(2,4),legend = False, color=['#fb9a99','#b2df8a','#fdbf6f','#a6cee3','#1d91c0','#33a02c','#1f78b4'])
    ax.set_yticks(ticks=[0,300,600,900,1200,1500])
    ax.set_yticklabels(labels=[0,300,600,900,1200,1500],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Human gut'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality bins', fontsize=15,color = 'black')
    plt.savefig('Real_human_hq.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()
    print(subset)
    subset = pd.DataFrame(np.array([[num_tara_maxbin2_single,num_tara_vamb_single,num_tara_metabat2_single,num_tara_semibin_single,num_tara_semibin_pretrain_single, num_tara_vamb_multi,num_tara_semibin_multi]]),columns = ['Maxbin2','VAMB(single)','Metabat2','SemiBin(single)','SemiBin(pre-train)','VAMB(multi)', 'SemiBin(multi)'], index=['Tara'])
    ax = subset.plot(kind='bar',figsize=(2,4),legend = False, color=['#fb9a99','#b2df8a','#fdbf6f','#a6cee3','#1d91c0','#33a02c','#1f78b4'])
    ax.set_yticks(ticks=[0,100,200,300,400,500])
    ax.set_yticklabels(labels=[0,100,200,300,400,500],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Ocean'], fontsize=15,color = 'black',rotation = 360)
    ax.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0., fontsize = 8)
    plt.savefig('Real_tara_hq.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()
    print(subset)




    # subset = pd.DataFrame(np.array([[skin_strain_metabat2,skin_strain_vamb,skin_strain_SemiBin],[oral_strain_metabat2,oral_strain_vamb,oral_strain_SemiBin]]),columns = ['Metabat2','VAMB','SemiBin'], index=['Skin','Oral'])
    # ax = subset.plot(kind='bar',width = 0.6,color = ['#fdbf6f', '#b2df8a', '#a6cee3'],figsize=(3,4))
    # ax.set_yticks(ticks=[0,20,40,60,80,100,120,140])
    # ax.set_yticklabels(labels=[0,20,40,60,80,100,120,140],fontsize=12,color = 'black')
    # ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    # ax.set_ylabel('High quality genomes', fontsize=15,color = 'black')
    # ax.set_title('strain', fontsize=20, alpha=1.0,color = 'black')
    # plt.savefig('CAMI_II_com_strain.pdf', dpi=300, bbox_inches='tight')
    # plt.close()
    #
    # subset = pd.DataFrame(np.array([[skin_species_metabat2,skin_species_vamb,skin_species_SemiBin],[oral_species_metabat2,oral_species_vamb,oral_species_SemiBin]]),columns = ['Metabat2','VAMB','SemiBin'], index=['Skin','Oral'])
    # ax = subset.plot(kind='bar',width= 0.6,legend = False,color = ['#fdbf6f', '#b2df8a', '#a6cee3'],figsize=(3,4))
    # ax.set_yticks(ticks=[0,15,30,45,60,75,90,105])
    # ax.set_yticklabels(labels=[0,15,30,45,60,75,90,105],fontsize=12,color = 'black')
    # ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    # #ax.set_ylabel('Num', fontsize=15,color = 'black')
    # ax.set_title('species', fontsize=20, alpha=1.0,color = 'black')
    # plt.savefig('CAMI_II_com_species.pdf', dpi=300, bbox_inches='tight')
    # plt.close()
    #
    # subset = pd.DataFrame(np.array([[skin_genus_metabat2,skin_genus_vamb,skin_genus_SemiBin],[oral_genus_metabat2,oral_genus_vamb,oral_genus_SemiBin]]),columns = ['Metabat2','VAMB','SemiBin'], index=['Skin','Oral'])
    # ax = subset.plot(kind='bar',width = 0.6,legend = False,color = ['#fdbf6f', '#b2df8a', '#a6cee3'],figsize=(3,4))
    # ax.set_yticks(ticks=[0,10,20,30,40,50,60])
    # ax.set_yticklabels(labels=[0,10,20,30,40,50,60],fontsize=12,color = 'black')
    # ax.set_xticklabels(labels=['Skin','Oral'], fontsize=15,color = 'black',rotation = 360)
    # #ax.set_ylabel('Num', fontsize=15,color = 'black')
    # ax.set_title('genus', fontsize=20, alpha=1.0,color = 'black')
    # plt.savefig('CAMI_II_com_genus.pdf', dpi=300, bbox_inches='tight')
    # plt.close()

    #plt.show()