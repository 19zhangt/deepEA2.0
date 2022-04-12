
kmc_dir=/storage_server/home/zhangt/freePEA/scripts/KMC/bin

input=~/a2z/m6A_13spp/workflow/data_12ss/ppa/01cleandata/ppa_Input1_R1.fastq.gz
ip=~/a2z/m6A_13spp/workflow/data_12ss/ppa/01cleandata/ppa_IP1_R1.fastq.gz

${kmc_dir}/kmc -t10 -k51 -ci20 -cs1000000 ${input} Input1_R1_k51.res .
${kmc_dir}/kmc -t10 -k51 -ci20 -cs1000000 ${ip} IP1_R1_k51.res .

${kmc_dir}/kmc_tools simple Input1_R1_k51.res -ci10 IP1_R1_k51.res -ci10 \
intersect Input_union -ocleft intersect IP_union -ocright

${kmc_dir}/kmc_dump Input_union Input_union 
${kmc_dir}/kmc_dump IP_union IP_union

# sort -nk2 -r 0db1_union_db2  > 0db1_union_db2.freq
