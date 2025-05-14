# Horizon.py
A program made to discover possibly Horizontallly transferred sequences in a query sequence using Dinucleotide frequency signatures

Horizontally Transferred insert detector


This documentation encompasses the workings and usage of the accompanying program.


Said program is made to find anomalous DNA loci in a given sequence, which can range from a sufficiently long gene cluster to an entire genome.


Inputs in order: 


Name of input file(str): Name of file containing input sequence. Can be fasta or txt


Number of bootstraps required(int):Dinucleotide signature of reference will be bootstrapped for provided number of times


Query window size(int): window size for examining the sequence


Sliding value of window(int): value over which the sliding window moves across the genome.


Significant deviation value(int): multiple of standard deviations from mean of obtained dinucleotide frequencies,  beyond which dinucleotide signatures will be considered relevantly varying from the host.


Name of output file(str): output for chi square graph, can be jpeg ,eps or other relevant image formats.


Outputs a list of relevant window start locations based on variation_value entered, full locus=start position+window_length, along with an image file containing a graph of calculated frequency variation values, larger peaks denote higher chances of being a horizontally transferred sequence.


Description


Evolution in sequences is more likely to affect larger sequence length than smaller ones, and as such, the 16 Dinucleotide frequencies within genomes have been found to be relatively conserved and distinct within several species, to the extent that they can be utilised as genomic signatures.


Moreover, these signatures have been found to be conserved for any 50kb subsample from the genome.


This provides a possible means to track Horizontal Gene transfer between organisms, as the foreign DNA would possess a distinct signature of Dinucleotide frequencies than the rest of the host genome.


This program queries substrings in a given DNA sequence for their Dinucleotide signature and compares them to the signature of the host sequence, outputting a list of start sites of substrings that present distinct signatures. 


The signature of the full sequence is calculated by subsampling 50kb regions at random and finding bootstrapped Dinucleotide frequencies.


Dinucleotide frequencies are defined as:


DF=f(KA)/f(K)*f(A)


Where KA is a Dinucleotide, f(KA) is the frequency of the Dinucleotide in a given sequence, and f(K) and f(A) represent frequencies of the individual mononucleotides forming the selected Dinucleotide.


The sequence is then queried using a sliding window. To retain detection capability for small anomalies and maintaining enough resolution to classify non anomalous regions, a window must be chosen based on expected length of anomalous inserts. The sliding length should also be selected based on required resolution and number of measurement points.


Dinucleotide frequencies are calculated for each window and then compared to the reference bootstrapped value for the whole sequence using Chi square values:(Observed value-Expected value)**2/Expected value. Larger Chi Square values denote a larger deviation from the reference.


As Chi square distributions approach bell curves at higher degrees of freedom, the loci deviating from mean by length larger than a selected value times standard deviation(> mean+(std*variation_value)) will considered as possible candidates for Horizontally Transferred Genes, and are output by the program based on the input variation_value.


The program also provides a graph of the calculated Frequency deviation scores at each queried locus, to point out areas of interest.

Input files providede are: actinomycetes_microcystis.txt and actinoplanes_mod.txt


References


* Tsirigos, A., & Rigoutsos, I. (2005). A new computational method for the detection of horizontal gene transfer events. Nucleic acids research, 33(3), 922-933.
* Baran, R. H., & Ko, H. (2008). Detecting horizontally transferred and essential genes based on dinucleotide relative abundance. DNA research, 15(5), 267-276.
* Von Wintersdorff, C. J., Penders, J., Van Niekerk, J. M., Mills, N. D., Majumder, S., Van Alphen, L. B., ... & Wolffs, P. F. (2016). Dissemination of antimicrobial resistance in microbial ecosystems through horizontal gene transfer. Frontiers in microbiology, 7, 173.
* Ravenhall, M., Škunca, N., Lassalle, F., & Dessimoz, C. (2015). Inferring horizontal gene transfer. PLoS computational biology, 11(5), e1004095.
* Baran, R. H., & Ko, H. (2008). Detecting horizontally transferred and essential genes based on dinucleotide relative abundance. DNA research, 15(5), 267-276.
* Williams Jr, C. A. (1950). The choice of the number and width of classes for the chi-square test of goodness of fit. Journal of the American Statistical Association, 45(249), 77-86.
