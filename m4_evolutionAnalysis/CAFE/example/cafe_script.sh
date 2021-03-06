#!shell
date

#specify data file, p-value threshold, # of threads to use, and log file
load -i example_data.tab -p 0.01 -t 10 -l log.txt

#the phylogenetic tree structure with branch lengths
tree (((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:93)

#search for 2 parameter model
lambda -s -t (((2,2)1,(1,1)1)1,1)

#specify the global lambda to for generating simulated data
lambda -l 0.0017

#generate 10 simulated data sets
genfamily rndtree/rnd -t 10

#estimate lambdas and compare likelihoods of global lambda and 2-parameter models
lhtest -d rndtree -l 0.0017 -t (((2,2)1,(1,1)1)1,1) -o lh2.out

date
