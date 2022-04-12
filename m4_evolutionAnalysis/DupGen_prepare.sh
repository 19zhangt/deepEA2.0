## background sequences
# cat genome/Aco/Aco.pep.fa genome/Spo/Spo.pep.fa genome/Atr/Atr.pep.fa genome/Mpo/Mpo.pep.fa > ${filedir}/Background/four_background.fa
# cat genome/Aco/Aco_DupGen.gff genome/Spo/Spo_DupGen.gff genome/Atr/Atr_DupGen.gff genome/Mpo/Mpo_DupGen.gff > ${filedir}/Background/four_background.gff
# makeblastdb -in ${filedir}/Background/four_background.fa -dbtype nucl -parse_seqids -out ${filedir}/Background/four_background_index
# diamond makedb --in ${filedir}/Background/four_background.fa -d ${filedir}/Background/pep_four_background_index

export PATH=/home/zhangt/a2z/m6A_13spp/tools/DupGen_finder:$PATH

genomedir=/distorage_server/home/zhangt/genomes
scriptsdir=/home/zhangt/a2z/m6A_13spp/tools
filedir=/home/zhangt/a2z/m6A_13spp/associated_data/Duplication_Genes


for species in pvu #ath gar ghi pvu gma sly sbi zma ata tdi tae osa ppa
do

    mkdir -p ${filedir}/${species}
    # blast with itself
    awk -v tmpn=${species} '{OFS="\t";print tmpn"_"$1,$4,$2,$3}' \
    ${genomedir}/${species^}/Annotation/${species^}.exons.gff3.gene.bed >${filedir}/${species}/${species}.gff

    diamond makedb --in ${genomedir}/${species^}/Pep/${species^}.longest.pep.rename.fa -d ${genomedir}/${species^}/Pep/${species^}.longest.pep.rename_index
    diamond blastp -d ${genomedir}/${species^}/Pep/${species^}.longest.pep.rename_index -q ${genomedir}/${species^}/Pep/${species^}.longest.pep.rename.fa -o ${filedir}/${species}/${species}_tmp.blast
    python ${scriptsdir}/DupGen_prepare.py ${filedir}/${species}/${species}_tmp.blast ${filedir}/${species}/${species}.blast
    
    # blast with background
    diamond blastp -d ${filedir}/Background/pep_four_background_index -q ${genomedir}/${species^}/Pep/${species^}.longest.pep.rename.fa -o ${filedir}/${species}/${species}_BG_tmp.blast
    python ${scriptsdir}/DupGen_prepare.py ${filedir}/${species}/${species}_BG_tmp.blast ${filedir}/${species}/${species}_BG.blast
    
    cat ${filedir}/${species}/${species}.gff ${filedir}/Background/four_background.gff > ${filedir}/${species}/${species}_BG.gff

    # makeblastdb -in pep_cds/cds/${species}.fa -dbtype nucl -parse_seqids -out pep_cds/cds/${species}_index
    # blastn -query pep_cds/cds/${species}.fa -db pep_cds/cds/${species}_index -out ${filedir}/${species}/${species}.blast -evalue 1e-5 -num_threads 10 -outfmt 6
    # blastn -query pep_cds/cds/${species}.fa -db ${filedir}/Background/four_background_index -out ${filedir}/${species}/${species}_BG.blast -evalue 1e-5 -num_threads 10 -outfmt 6
    # mkdir ${filedir}/${species}/pep_raw_DupGen_results    
    # ${scriptsdir}/DupGen_finder/DupGen_finder.pl -i ${filedir}/${species} -t ${species} -c BG -o ${filedir}/${species}/pep_raw_DupGen_results -m 20 -e 1e-10


    ${scriptsdir}/DupGen_finder/DupGen_finder-unique.pl -i ${filedir}/${species} -t ${species} -c BG -o ${filedir}/${species}/pep_DupGen_results -m 20 -e 1e-10
done
