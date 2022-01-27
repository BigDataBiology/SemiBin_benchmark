import pandas as pd
import numpy as np

human_list = ['CCMD75147712ST', 'CCMD66848156ST', 'CCMD65222621ST', 'CCMD35801800ST', 'CCMD56948710ST',
              'CCMD30627121ST', 'CCMD95676152ST', 'CCMD76409700ST', 'CCMD30626189ST', 'CCMD99929634ST']

dog_list = ['SAMN06172516', 'SAMN06172518', 'SAMN06172405', 'SAMN06172415', 'SAMN06172478', 'SAMN06172403',
            'SAMN06172439', 'SAMN06172442', 'SAMN06172504', 'SAMEA103957794']

ocean_list = ['TARA_148_SRF_0.22-3', 'TARA_128_SRF_0.22-3', 'TARA_066_SRF_lt-0.22', 'TARA_031_SRF_0.22-1.6',
              'TARA_036_SRF_0.22-1.6', 'TARA_048_SRF_0.22-1.6', 'TARA_078_SRF_0.45-0.8', 'TARA_023_SRF_0.22-1.6',
              'TARA_068_SRF_lt-0.22', 'TARA_018_SRF_lt-0.22']


def process_results(log_file):
    data = []

    with open(log_file) as fh:
        row = {}
        for line in fh:
            if line.strip().startswith("Command exited with non-zero"):
                continue

            name, value = line.strip().split(": ")
            value = value.strip('"')

            if name == "Elapsed (wall clock) time (h:mm:ss or m:ss)":
                if '.' in value:
                    value = value.split(":")
                    value[-1] = round(float(value[-1]))
                else:
                    value = value.split(":")

                if len(value) == 2:
                    seconds = int(value[0]) * 60 + int(value[1])
                elif len(value) == 3:
                    seconds = int(value[0]) * 3600 + int(value[1]) * 60 + int(value[2])
                else:
                    raise Exception("Unknown time format {}".format(value))
                value = seconds
                row['Time(s)'] = value

            if name == 'Maximum resident set size (kbytes)':
                row['Memory(Mb)'] = int(value)/1024

    return row

def process_result(environment, mode):
    generate_data = {'time':[], 'Mem':[]}
    generate_cannot = {'time':[], 'Mem':[]}
    train_cpu = {'time':[], 'Mem':[]}
    train_gpu = {'time':[], 'Mem':[]}
    bin = {'time':[], 'Mem':[]}

    if environment == 'human':
        sample_list = human_list
    if environment == 'dog':
        sample_list = dog_list
    if environment == 'ocean':
        sample_list = ocean_list

    for temp in sample_list:
        # generate_data
        if mode == 'single':
            result = process_results('updated_results/log_SemiBin_0.5/{1}/{0}/{2}/generate_data_single.txt'.format(mode,environment,temp))
            generate_data['time'].append(result['Time(s)'])
            generate_data['Mem'].append(result['Memory(Mb)'])

        if mode == 'single':
        # generate_cannot
            result = process_results('Results/computation_log/{0}/{1}/{2}/taxonomy_log.txt'.format(mode,environment,temp))
            generate_cannot['time'].append(result['Time(s)'])
            generate_cannot['Mem'].append(result['Memory(Mb)'])

        # train cpu
        result = process_results('Results/computation_log/{0}/{1}/{2}/train_log.txt'.format(mode,environment,temp))
        train_cpu['time'].append(result['Time(s)'])
        train_cpu['Mem'].append(result['Memory(Mb)'])

        # train gpu
        result = process_results('Results/computation_log/{0}/{1}/{2}/train_log_gpu.txt'.format(mode,environment,temp))
        train_gpu['time'].append(result['Time(s)'])
        train_gpu['Mem'].append(result['Memory(Mb)'])

        # bin
        if mode == 'single':
            result = process_results('updated_results/log_SemiBin_0.5/{1}/{0}/{2}/bin_single.txt'.format(mode,environment,temp))
            bin['time'].append(result['Time(s)'])
            bin['Mem'].append(result['Memory(Mb)'])
        else:
            result = process_results('updated_results/log_SemiBin_0.5/{1}/{0}/{2}/bin_multi.txt'.format(mode,environment,temp))
            bin['time'].append(result['Time(s)'])
            bin['Mem'].append(result['Memory(Mb)'])

    if mode == 'multi':
        result = process_results('updated_results/log_SemiBin_0.5/{1}/{0}/generate_data_multi.txt'.format(mode,environment))
        generate_data['time'].append(result['Time(s)'])
        generate_data['Mem'].append(result['Memory(Mb)'])

    print('{0} {1} Time(minutes) and mem(MB)'.format(environment, mode))
    print('Generating data:', np.mean(generate_data['time'])/60, np.mean(generate_data['Mem']))
    print(generate_data)
    if mode == 'single':
        print('Generating cannot:', np.mean(generate_cannot['time'])/60, np.mean(generate_cannot['Mem']))
        print(generate_cannot)
    print('train cpu:', np.mean(train_cpu['time'])/60, np.mean(train_cpu['Mem']))
    print(train_cpu)
    print('train gpu:', np.mean(train_gpu['time'])/60, np.mean(train_gpu['Mem']))
    print(train_gpu)
    print('bin:', np.mean(bin['time'])/60, np.mean(bin['Mem']))
    print(bin)

    return generate_data, generate_cannot, train_cpu, train_gpu, bin


def sum_SemiBin(environment, mode):
    generate_data = {'time':0, 'Mem':[]}
    generate_cannot = {'time':0, 'Mem':[]}
    train_cpu = {'time':0, 'Mem':[]}
    train_gpu = {'time':0, 'Mem':[]}
    bin = {'time':0, 'Mem':[]}

    if environment == 'human':
        sample_list = human_list
    if environment == 'dog':
        sample_list = dog_list
    if environment == 'ocean':
        sample_list = ocean_list

    for temp in sample_list:
        # generate_data
        if mode == 'single':
            result = process_results('updated_results/log_SemiBin_0.5/{1}/{0}/{2}/generate_data_single.txt'.format(mode,environment,temp))
            generate_data['time'] += result['Time(s)'] / 60
            generate_data['Mem'].append(result['Memory(Mb)'])
        # if mode == 'single':
        # generate_cannot
        result = process_results('Results/computation_log/{0}/{1}/{2}/taxonomy_log.txt'.format('single',environment,temp))
        generate_cannot['time'] += result['Time(s)'] / 60
        generate_cannot['Mem'].append(result['Memory(Mb)'])

        # train cpu
        result = process_results('Results/computation_log/{0}/{1}/{2}/train_log.txt'.format(mode,environment,temp))
        train_cpu['time'] += result['Time(s)'] / 60
        train_cpu['Mem'].append(result['Memory(Mb)'])

        # train gpu
        result = process_results('Results/computation_log/{0}/{1}/{2}/train_log_gpu.txt'.format(mode,environment,temp))
        train_gpu['time'] += result['Time(s)'] / 60
        train_gpu['Mem'].append(result['Memory(Mb)'])

        # bin
        if mode == 'single':
            result = process_results('updated_results/log_SemiBin_0.5/{1}/{0}/{2}/bin_single.txt'.format(mode,environment,temp))
            bin['time'] += result['Time(s)'] / 60
            bin['Mem'].append(result['Memory(Mb)'])
        else:
            result = process_results('updated_results/log_SemiBin_0.5/{1}/{0}/{2}/bin_multi.txt'.format(mode,environment,temp))
            bin['time'] += result['Time(s)'] / 60
            bin['Mem'].append(result['Memory(Mb)'])

    if mode == 'multi':
        result = process_results('updated_results/log_SemiBin_0.5/{1}/{0}/generate_data_multi.txt'.format(mode,environment))
        generate_data['time'] += result['Time(s)'] / 60
        generate_data['Mem'].append(result['Memory(Mb)'])

    generate_data['Mem'] = np.max(generate_data['Mem'])

    generate_cannot['Mem'] = np.max(generate_cannot['Mem'])
    train_cpu['Mem'] = np.max(train_cpu['Mem'])
    train_gpu['Mem'] = np.max(train_gpu['Mem'])
    bin['Mem'] = np.max(bin['Mem'])

    print('{0} {1} Time(minutes) and mem(MB)'.format(environment, mode))
    if mode == 'multi':
        print('SemiBin(CPU), multi-sample', np.sum((generate_data['time'], generate_cannot['time'], train_cpu['time'], bin['time'])), np.max((generate_data['Mem'], generate_cannot['Mem'], train_cpu['Mem'],  bin['Mem'])))
        print('SemiBin(GPU), multi-sample', np.sum((generate_data['time'], generate_cannot['time'], train_gpu['time'], bin['time'])), np.max((generate_data['Mem'], generate_cannot['Mem'],
                     train_gpu['Mem'], bin['Mem'])))
    else:
        print('SemiBin(CPU), single-sample', np.sum((generate_data['time'], generate_cannot['time'], train_cpu['time'], bin['time'])), np.max((generate_data['Mem'], generate_cannot['Mem'], train_cpu['Mem'],  bin['Mem'])))
        print('SemiBin(GPU), single-sample', np.sum((generate_data['time'], generate_cannot['time'], train_gpu['time'], bin['time'])), np.max((generate_data['Mem'], generate_cannot['Mem'], train_gpu['Mem'], bin['Mem'])))
        print('SemiBin(pretrain), single-sample', np.sum((generate_data['time'], bin['time'])), np.max((generate_data['Mem'], bin['Mem'])))

def sum_Maxbin2(environment):
    bin = {'time':0, 'Mem':[]}

    if environment == 'human':
        sample_list = human_list
    if environment == 'dog':
        sample_list = dog_list
    if environment == 'ocean':
        sample_list = ocean_list

    for temp in sample_list:
        result = process_results('updated_results/running_log/{0}/single/{1}/maxbin2_log.txt'.format(environment, temp))
        print(result)
        bin['time'] += result['Time(s)'] / 60
        bin['Mem'].append(result['Memory(Mb)'])

    bin['Mem'] = np.max(bin['Mem'])
    print(bin)

def sum_Metabat2(environment):
    bin = {'time':0, 'Mem':[]}

    if environment == 'human':
        sample_list = human_list
    if environment == 'dog':
        sample_list = dog_list
    if environment == 'ocean':
        sample_list = ocean_list

    for temp in sample_list:
        result = process_results('updated_results/running_log/{0}/single/{1}/metabat2_log.txt'.format(environment, temp))
        print(result)
        bin['time'] += result['Time(s)'] / 60
        bin['Mem'].append(result['Memory(Mb)'])

    bin['Mem'] = np.max(bin['Mem'])
    print(bin)

def sum_Vamb(environment):
    bin = {'time':0, 'Mem':[]}

    if environment == 'human':
        sample_list = human_list
    if environment == 'dog':
        sample_list = dog_list
    if environment == 'ocean':
        sample_list = ocean_list

    for temp in sample_list:
        result = process_results('updated_results/running_log/{0}/single/{1}/vamb_cpu_single.txt'.format(environment, temp))
        print(result)
        bin['time'] += result['Time(s)'] / 60
        bin['Mem'].append(result['Memory(Mb)'])

    bin['Mem'] = np.max(bin['Mem'])
    print('single-sample(CPU)', bin)
    print('\n\n')

    bin = {'time':0, 'Mem':[]}

    for temp in sample_list:
        result = process_results('updated_results/running_log/{0}/single/{1}/vamb_gpu_single.txt'.format(environment, temp))
        print(result)
        bin['time'] += result['Time(s)'] / 60
        bin['Mem'].append(result['Memory(Mb)'])

    bin['Mem'] = np.max(bin['Mem'])
    print('single-sample(GPU)', bin)
    print('\n\n')

    bin = {'time':0, 'Mem':[]}

    result = process_results('updated_results/running_log/{0}/multi/vamb_cpu_multi.txt'.format(environment))
    print(result)
    bin['time'] += result['Time(s)'] / 60
    bin['Mem'] = result['Memory(Mb)']
    print('multi-sample(CPU)', bin)
    print('\n\n')

    bin = {'time':0, 'Mem':[]}

    result = process_results('updated_results/running_log/{0}/multi/vamb_gpu_multi.txt'.format(environment))
    print(result)
    bin['time'] += result['Time(s)'] / 60
    bin['Mem'] = result['Memory(Mb)']
    print('multi-sample(GPU)', bin)

if __name__ == '__main__':
    # sum_SemiBin('dog', 'single')
    # sum_SemiBin('dog', 'multi')
    # sum_SemiBin('ocean', 'single')
    # sum_SemiBin('ocean', 'multi')
    # sum_SemiBin('human', 'single')
    sum_SemiBin('human', 'multi')

    # sum_Vamb('human')
    # sum_Vamb('dog')
    # sum_Vamb('ocean')

    # process_result('dog','single')
    # process_result('dog', 'multi')
    # process_result('human', 'single')
    # process_result('human', 'multi')
    # process_result('ocean', 'single')
    # process_result('ocean', 'multi')

    # sum_Maxbin2('ocean')

    # sum_Metabat2('ocean')
