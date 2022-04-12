#!shell
date
version

#specify data file, p-value threshold, # of threads to use, and log file
load -i example_data.tab -p 0.01 -t 10 -l log.txt

#the phylogenetic tree structure with branch lengths
tree (((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:93)

#search for 2 parameter model
lambdamu -s -t (((2,2)1,(1,1)1)1,1)

# generate a report
report report.txt

