#! usr/bin/env python
import sys

'''Program to find anomolous loci in given sequences and outputs list of candidate Horizontally Transferred sequences based on variation in Dinucleotide Frequency profiles, using a custon vhi square based formula.
Inputs in order: name of input file(str),number of bootstraps required(int), quiry window size(int), sliding value of window(int), multiple of standard deviations of mean of chi square value beyond which deviations of dinucleotide signatures will be considered relevant(int), name of output file(str)
References:
Tsirigos, A., & Rigoutsos, I. (2005). A new computational method for the detection of horizontal gene transfer events. Nucleic acids research, 33(3), 922-933.
Baran, R. H., & Ko, H. (2008). Detecting horizontally transferred and essential genes based on dinucleotide relative abundance. DNA research, 15(5), 267-276.
Von Wintersdorff, C. J., Penders, J., Van Niekerk, J. M., Mills, N. D., Majumder, S., Van Alphen, L. B., ... & Wolffs, P. F. (2016). Dissemination of antimicrobial resistance in microbial ecosystems through horizontal gene transfer. Frontiers in microbiology, 7, 173.
Ravenhall, M., Škunca, N., Lassalle, F., & Dessimoz, C. (2015). Inferring horizontal gene transfer. PLoS computational biology, 11(5), e1004095.
Baran, R. H., & Ko, H. (2008). Detecting horizontally transferred and essential genes based on dinucleotide relative abundance. DNA research, 15(5), 267-276.'''

def main():
    input_file = sys.argv[1]
    strap_ref = int(sys.argv[2])
    window_size = int(sys.argv[3])
    slide_value = int(sys.argv[4])
    variation_value = int(sys.argv[5])
    output_file = sys.argv[6]
    m9(input_file, strap_ref, window_size, slide_value, variation_value, output_file)

def read_sequence(input_file):
    file = open(input_file, 'r')
    seqs = file.readlines()
    ref = ''
    for t in range(0, len(seqs)):       #reads sequence
        if seqs[t][0] != '>':
            ref = ref + seqs[t].strip()
    return ref

def bootstrap_frequencies(ref, strap_ref, mono_list, dinuc_list):
    import random
    import numpy
    sample_dinuc = numpy.zeros((strap_ref, len(dinuc_list)))
    mono_arr = numpy.zeros((strap_ref, len(mono_list)))
    t = 0
    while t < strap_ref:
        value = random.randrange(len(ref) - 50001)               #Bootstraps dinucleotide fequencies using dictionaries of nucleotides
        sample = ref[value:value+50001]
        if sample.count('N') / len(sample) >= 0.1:
            continue
        temp_dit = {}
        mono_dit = {}
        for q in sample:
            if q == 'N':
                continue
            if q not in mono_dit and len(q.strip()) == 1:
                mono_dit[q] = 1
            elif q in mono_dit and len(q.strip()) == 1:
                mono_dit[q] = mono_dit[q] + 1
        for p in range(len(mono_list)):
            mono_arr[t, p] = mono_dit[mono_list[p]]                      #forms dictionaries of nucleotides
        for a in range(len(sample) - 1):
            dinuc = sample[a:a+2]
            if dinuc.count('N') > 0:
                continue
            if dinuc not in temp_dit and len(dinuc.strip()) == 2:
                temp_dit[dinuc] = 1
            elif dinuc in temp_dit and len(dinuc.strip()) == 2:
                temp_dit[dinuc] = temp_dit[dinuc] + 1
        for x in range(len(dinuc_list)):
            sample_dinuc[t, x] = temp_dit[dinuc_list[x]]
        t = t + 1
    return mono_arr, sample_dinuc

def calculate_final_frequencies(mono_arr, sample_dinuc, mono_list, dinuc_list):
    import numpy
    final_sample_mono = numpy.zeros((1, len(mono_list)))
    final_sample_dinuc = numpy.zeros((1, len(dinuc_list)))
    for a in range(len(mono_list)):
        final_sample_mono[0, a] = numpy.mean(mono_arr[:, a])         #forms arrays containing bootstrapped values
    mono_dit2 = {mono_list[p]: final_sample_mono[0, p] for p in range(len(mono_list))}
    for x in range(len(dinuc_list)):
        final_sample_dinuc[0, x] = numpy.mean(sample_dinuc[:, x])
        for t in range(len(dinuc_list[x])):
            final_sample_dinuc[0, x] = final_sample_dinuc[0, x] / mono_dit2[dinuc_list[x][t]]     #finds Dinucleotide freqency
    return final_sample_dinuc

def sliding_window_analysis(ref, window_size, slide_value, final_sample_dinuc, mono_list, dinuc_list):
    import numpy
    initial = 0
    final_list = []
    initial_list = []                                                 
    while initial < len(ref) - window_size:                           #queries the sequence using a sliding window
        test_sample = ref[initial:initial+window_size]              #Form test window for chi square calculation
        if test_sample.count('N') / len(test_sample) > 0.1:           #If more than 10% of the test window is Ns, discard the window and move 50 nucleotides ahead
            initial = initial + 50
            continue
        test_arr = numpy.zeros((1, len(dinuc_list)))                 #Initialise array for dinucleatide percentage storage for test window
        mono_dit = {}
        for q in test_sample:                                     #Make a mononucleotide dictionary
            if q not in mono_dit:
                mono_dit[q] = 1
            elif q in mono_dit:
                mono_dit[q] = mono_dit[q] + 1
        if 'N' in mono_dit:                                          #Remove uncertain bases from mononucleotide dictionary
            mono_dit.pop('N')
        temp_dit = {}
        for a in range(len(test_sample)):                                #forms dictionaries of dinucleotides
            dinuc = test_sample[a:a+2]
            if dinuc.count('N') > 0:                                 #Remove uncertain bases from dinucleotide dictionary
                continue
            if dinuc not in temp_dit and len(dinuc.strip()) == 2:
                temp_dit[dinuc] = 1
            elif dinuc in temp_dit and len(dinuc.strip()) == 2:
                temp_dit[dinuc] = temp_dit[dinuc] + 1
        for x in range(len(dinuc_list)):                                  #Fill an array with dinucleotide counts                    
            test_arr[0, x] = temp_dit[dinuc_list[x]]
            for p in range(len(dinuc_list[x])):                           #Divide dinuc counts by mono counts to find dinuc frequency
                test_arr[0, x] = test_arr[0, x] / mono_dit[dinuc_list[x][p]]
        final1 = 0
        for gra in range(len(dinuc_list)):
            temp1 = ((final_sample_dinuc[0, gra] - test_arr[0, gra]) ** 2)              #calculates a chi2 value for each window
            temp1 = temp1 / final_sample_dinuc[0, gra]
            final1 = final1 + temp1
        final_list.append(final1)
        initial_list.append(initial)
        initial = initial + slide_value
    return final_list, initial_list

def plot_scores(initial_list, final_list, output_file):
    from matplotlib import pyplot
    import numpy
    pyplot.plot(numpy.array(initial_list), numpy.array(final_list))              #plots calculated values of each window as observed at specific locations
    pyplot.minorticks_on()
    pyplot.xlabel('window start position')
    pyplot.ylabel('Signature Distinction Score')
    pyplot.savefig(output_file)

def detect_anomalies(final_list, initial_list, variation_value):
    import numpy
    final_arr = numpy.array(final_list)
    final_mean = numpy.mean(final_arr)
    final_std = numpy.std(final_arr)
    var_list = []
    for t in range(len(final_list)):
        if final_list[t] >= (final_mean + (variation_value * final_std)):              #finds relevant loci based on score variation from mean
            var_list.append(initial_list[t])
    print(var_list)

def m9(input_file, strap_ref, window_size, slide_value, variation_value, output_file):
    mono_list = ['T', 'C', 'G', 'A']
    dinuc_list = ['AT','AG','AC','AA','GC','GG','GA','GT','CA','CC','CT','CG','TA','TT','TG','TC']
    ref = read_sequence(input_file)
    mono_arr, sample_dinuc = bootstrap_frequencies(ref, strap_ref, mono_list, dinuc_list)
    final_sample_dinuc = calculate_final_frequencies(mono_arr, sample_dinuc, mono_list, dinuc_list)
    final_list, initial_list = sliding_window_analysis(ref, window_size, slide_value, final_sample_dinuc, mono_list, dinuc_list)
    plot_scores(initial_list, final_list, output_file)
    detect_anomalies(final_list, initial_list, variation_value)

main()
