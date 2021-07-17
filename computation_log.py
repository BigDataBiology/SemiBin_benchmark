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
            result = process_results('Results/computation_log/{0}/{1}/{2}/generate_data_log.txt'.format(mode,environment,temp))
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
        result = process_results('Results/computation_log/{0}/{1}/{2}/bin_log.txt'.format(mode,environment,temp))
        bin['time'].append(result['Time(s)'])
        bin['Mem'].append(result['Memory(Mb)'])

    if mode == 'multi':
        result = process_results('Results/computation_log/{0}/{1}/generate_data_log.txt'.format(mode,environment))
        generate_data['time'].append(result['Time(s)'])
        generate_data['Mem'].append(result['Memory(Mb)'])

    print('{0} {1} Time(minutes) and mem(MB)'.format(environment, mode))
    print('Generating data:', np.mean(generate_data['time'])/60, np.mean(generate_data['Mem']))
    if mode == 'single':
        print('Generating cannot:', np.mean(generate_cannot['time'])/60, np.mean(generate_cannot['Mem']))
    print('train cpu:', np.mean(train_cpu['time'])/60, np.mean(train_cpu['Mem']))
    print('train gpu:', np.mean(train_gpu['time'])/60, np.mean(train_gpu['Mem']))
    print('bin:', bin['time'], np.mean(bin['time'])/60, np.mean(bin['Mem']))

    return  generate_data, generate_cannot, train_cpu, train_gpu, bin

if __name__ == '__main__':
    generate_data, generate_cannot, train_cpu, train_gpu, bin = process_result('human', 'multi')

