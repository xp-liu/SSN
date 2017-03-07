# SSN
This is source code for constructing SSN in python, and the “SciPy” package is required for this program.

The "construct_single_network.py" can be used to construct the single-sample network (SSN), and the parameter is shown below:

usage: python construct_single_network.py -pvalue=threshold_p_value | -threshold=threshold_value -background=background_network_file -ref=reference_sample_file  -sample=sample_data_file -out=results_output_fold

Options and arguments:

-pvalue : set the threshold of p-value [0..1], if the -pvalue set 1, all edges will be outputted to the SSN

-threshold : set the threshold value of the absolute value of deltaPCC [0..2]

-background : background network to calculate the deltaPCC of edges based on the network

-ref : the expression profile of reference samples

-sample : the expression profile for the sample to be constructed the SSN

-out : the directory to store the SSN



There is a simple example below: 

The file "reference_samples.txt" is the expression profile of reference samples.

The file "background.txt" is the background network for calculating SSN, and the user can use the co-expression network, existing regulation network, or a sub-network in here.

The file "10_samples_data.txt" is a 10 samples expression profile, the SSNs will be constructed for the 10 samples profile.


To construct the SSN for the 10 samples, and the results will be put in "target" directory:

python construct_single_network.py -ref=reference_samples.txt -background=background.txt -sample=10_samples_data.txt -pvalue=0.01 -out=target

or 

python construct_single_network.py -ref=reference_samples.txt -background=background.txt -sample=10_samples_data.txt -threshold=0.5 -out=target


