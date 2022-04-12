# Snakemake file for MeRIP-seq analysis
# conda activate metaPlants
# ~/miniconda3/bin/snakemake -s m6A_running_12ss.smk --configfile run.yaml -np

# fastp v0.20.1
# hisat2 v2.2.1
# featureCounts v2.0.1

DATA_DIR     =  config["data_dir"]
SPECIES_DIR   =  config["species_dir"]

## GFF3 to GTF
GFF3_FILE   = '{0}/{1}/Annotation/{1}.exons.gff3'.format(SPECIES_DIR, config["species_name"])
GTF_FILE    = '{0}/{1}/Annotation/{1}.exons.gtf'.format(SPECIES_DIR, config["species_name"])
cDNA_FILE   = '{0}/{1}/cDNA/{1}.cDNA.fa'.format(SPECIES_DIR, config["species_name"])
GENOME_FILE = '{0}/{1}/Genome/{1}.fa'.format(SPECIES_DIR, config["species_name"])
GENOME_SIZE = '{0}/{1}/Genome/{1}.genome.size'.format(SPECIES_DIR, config["species_name"])
HISAT_INDEX = '{0}/{1}/Genome/hisat2_index'.format(SPECIES_DIR, config["species_name"])

if not os.path.exists( GTF_FILE ):
    from bioinfokit.analys import gff
    gff.gff_to_gtf(file= GFF3_FILE)
    os.system('mv {} {}'.format(os.path.basename(GTF_FILE), GTF_FILE))
else:
    print("GTF file exists")

if not os.path.exists( GFF3_FILE + ".gene.bed" ) or os.path.getsize( GFF3_FILE + ".gene.bed" ) == 0:
    os.system('awk -F"\t" \'$3~/gene/{{OFS="\t";split($9,a,"=|;");gsub(/gene:/, "", a[2]);print $1,$4-1,$5,a[2],$6,$7,$3}}\' '+ GFF3_FILE + ' >' + GFF3_FILE + ".gene.bed")

if not os.path.exists( cDNA_FILE ):
    os.system('gffread -g {} -w {} {}'.format(GENOME_FILE, cDNA_FILE, GFF3_FILE))

if not os.path.exists( GENOME_SIZE ) or os.path.getsize(GENOME_SIZE) == 0:
    os.system('cut -f 1-2 {}.fai > {}'.format(GENOME_FILE, GENOME_SIZE))

if not os.path.exists( HISAT_INDEX ):
    os.system("mkdir -p " + HISAT_INDEX)
    os.system("extract_splice_sites.py {0} > {0}.genome.ss".format(GTF_FILE))
    os.system("extract_exons.py {0} > {0}.genome.exon".format(GTF_FILE))
    os.system("hisat2-build -p 20 --ss {0}.genome.ss --exon {0}.genome.exon {1} {2}/Genome".format(
        GTF_FILE, GENOME_FILE, HISAT_INDEX))


## workflow
rule zero:
    input:
        # expand(DATA_DIR + "00rawdata/{sample}.{treat}.{rep}.strand_info.txt", sample=config["samples"], rep = config["replicates"], treat = ["input", "ip"]),
        # expand(DATA_DIR + "final_output/read_density.png"),
        expand(DATA_DIR + "08fasta/{sample}_peaks_T.fa", sample=config["samples"])

if config["se_or_pe"] == "SE":
    rule fastp_se:
        input: DATA_DIR + "00rawdata/{sample}.{treat}.{rep}.fastq.gz"
        output: 
            reads = DATA_DIR + "01cleandata/{sample}.{treat}.{rep}.clean.fastq.gz",
            json = DATA_DIR + "01QC/fastp/{sample}.{treat}.{rep}_fastp.json",
            html = DATA_DIR + "01QC/fastp/{sample}.{treat}.{rep}_fastp.html"
        threads: 1
        params:
            datadir = DATA_DIR,
            mapt = 16,
            adapter=config["adapter1"]
        shell:
            """
            phred_value=`bash ./workflow/scripts/base_phred.sh {input}`
            if [[ ${{phred_value}} == 64 ]];then
                if [[ {params.adapter} == "auto" ]];then
                    fastp -w {params.mapt} --phred64 -5 -3 -r -n 0 -e 20 -i {input} -o {output.reads} \
                    -j {output.json} -h {output.html} --trim_poly_x --poly_x_min_len 15 \
                    -R "{wildcards.sample}.{wildcards.treat}.{wildcards.rep}"
                else
                    fastp -w {params.mapt} --phred64 -5 -3 -r -n 0 -e 20 -i {input} -o {output.reads} \ 
                    --adapter_sequence {params.adapter} \
                    -j {output.json} -h {output.html} --trim_poly_x --poly_x_min_len 15 \
                    -R "{wildcards.sample}.{wildcards.treat}.{wildcards.rep}"
                fi
            else
                if [[ {params.adapter} == "auto" ]];then
                    fastp -w {params.mapt} -5 -3 -r -n 0 -e 20 -i {input} -o {output.reads} \
                    -j {output.json} -h {output.html} --trim_poly_x --poly_x_min_len 15 \
                    -R "{wildcards.sample}.{wildcards.treat}.{wildcards.rep}"
                else
                    fastp -w {params.mapt} -5 -3 -r -n 0 -e 20 -i {input} -o {output.reads} \ 
                    --adapter_sequence {params.adapter} \
                    -j {output.json} -h {output.html} --trim_poly_x --poly_x_min_len 15 \
                    -R "{wildcards.sample}.{wildcards.treat}.{wildcards.rep}"
                fi
            fi
            """
elif config["se_or_pe"] == "PE":
    rule fastp_pe:
        input:
            reads1 = DATA_DIR + "00rawdata/{sample}.{treat}.{rep}_1.fastq.gz",
            reads2 = DATA_DIR + "00rawdata/{sample}.{treat}.{rep}_2.fastq.gz"
        output:
            creads1 = DATA_DIR + "01cleandata/{sample}.{treat}.{rep}_1.clean.fastq.gz",
            creads2 = DATA_DIR + "01cleandata/{sample}.{treat}.{rep}_2.clean.fastq.gz",
            json = DATA_DIR + "01QC/fastp/{sample}.{treat}.{rep}_fastp.json",
            html = DATA_DIR + "01QC/fastp/{sample}.{treat}.{rep}_fastp.html"
        threads: 1
        params:
            datadir = DATA_DIR,
            mapt = 16,
            adapter1=config["adapter1"],
            adapter2=config["adapter2"]
        shell:
            """
            phred_value=`bash ./workflow/scripts/base_phred.sh {input.reads1}`
            if [[ ${{phred_value}} == 64 ]];then
                if [[ {params.adapter1} == "auto" ]];then
                    fastp -w {params.mapt} --phred64 -i {input.reads1} -o {output.creads1} \
                    -I {input.reads2} -O {output.creads2} \
                    -j {output.json} -h {output.html} -5 -3 -r -n 0 -e 20 --trim_poly_x --poly_x_min_len 15 \
                    --detect_adapter_for_pe \
                    -R "{wildcards.sample}.{wildcards.treat}.{wildcards.rep}"
                else
                    fastp -w {params.mapt} --phred64 -i {input.reads1} -o {output.creads1} \
                    -I {input.reads2} -O {output.creads2} \
                    -j {output.json} -h {output.html} -5 -3 -r -n 0 -e 20 --trim_poly_x --poly_x_min_len 15 \
                    --adapter_sequence={params.adapter1} --adapter_sequence_r2={params.adapter2} \
                    -R "{wildcards.sample}.{wildcards.treat}.{wildcards.rep}"
                fi
            else
                if [[ {params.adapter1} == "auto" ]];then
                    fastp -w {params.mapt} -i {input.reads1} -o {output.creads1} \
                    -I {input.reads2} -O {output.creads2} \
                    -j {output.json} -h {output.html} -5 -3 -r -n 0 -e 20 --trim_poly_x --poly_x_min_len 15 \
                    --detect_adapter_for_pe \
                    -R "{wildcards.sample}.{wildcards.treat}.{wildcards.rep}"
                else
                    fastp -w {params.mapt} -i {input.reads1} -o {output.creads1} \
                    -I {input.reads2} -O {output.creads2} \
                    -j {output.json} -h {output.html} -5 -3 -r -n 0 -e 20 --trim_poly_x --poly_x_min_len 15 \
                    --adapter_sequence={params.adapter1} --adapter_sequence_r2={params.adapter2} \
                    -R "{wildcards.sample}.{wildcards.treat}.{wildcards.rep}"
                fi
            fi
            """
else:
    sys.exit()


if config["se_or_pe"] == "SE":
    rule fastp_se:
        input: DATA_DIR + "00rawdata/{sample}.{treat}.{rep}.fastq.gz"
        output: 
            reads = DATA_DIR + "01cleandata/{sample}.{treat}.{rep}.clean.fastq.gz",
            json = DATA_DIR + "01QC/fastp/{sample}.{treat}.{rep}_fastp.json",
            html = DATA_DIR + "01QC/fastp/{sample}.{treat}.{rep}_fastp.html"
        threads: 1
        params:
            datadir = DATA_DIR,
            mapt = 16,
            adapter=config["adapter1"]
        shell:
            """
            phred_value=`bash ./workflow/scripts/base_phred.sh {input}`
            if [[ ${{phred_value}} == 64 ]];then
                if [[ {params.adapter} == "auto" ]];then
                    fastp -w {params.mapt} --phred64 -5 -3 -r -n 0 -e 20 -i {input} -o {output.reads} \
                    -j {output.json} -h {output.html} --trim_poly_x --poly_x_min_len 15 \
                    -R "{wildcards.sample}.{wildcards.treat}.{wildcards.rep}" \
                    -U --umi_loc=per_read --umi_len=9
                else
                    fastp -w {params.mapt} --phred64 -5 -3 -r -n 0 -e 20 -i {input} -o {output.reads} \ 
                    --adapter_sequence {params.adapter} \
                    -j {output.json} -h {output.html} --trim_poly_x --poly_x_min_len 15 \
                    -R "{wildcards.sample}.{wildcards.treat}.{wildcards.rep}" \
                    -U --umi_loc=per_read --umi_len=9
                fi
            else
                if [[ {params.adapter} == "auto" ]];then
                    fastp -w {params.mapt} -5 -3 -r -n 0 -e 20 -i {input} -o {output.reads} \
                    -j {output.json} -h {output.html} --trim_poly_x --poly_x_min_len 15 \
                    -R "{wildcards.sample}.{wildcards.treat}.{wildcards.rep}" \
                    -U --umi_loc=per_read --umi_len=9
                else
                    fastp -w {params.mapt} -5 -3 -r -n 0 -e 20 -i {input} -o {output.reads} \ 
                    --adapter_sequence {params.adapter} \
                    -j {output.json} -h {output.html} --trim_poly_x --poly_x_min_len 15 \
                    -R "{wildcards.sample}.{wildcards.treat}.{wildcards.rep}" \
                    -U --umi_loc=per_read --umi_len=9
                fi
            fi
            """
elif config["se_or_pe"] == "PE":
    rule fastp_pe:
        input:
            reads1 = DATA_DIR + "00rawdata/{sample}.{treat}.{rep}_1.fastq.gz",
            reads2 = DATA_DIR + "00rawdata/{sample}.{treat}.{rep}_2.fastq.gz"
        output:
            creads1 = DATA_DIR + "01cleandata/{sample}.{treat}.{rep}_1.clean.fastq.gz",
            creads2 = DATA_DIR + "01cleandata/{sample}.{treat}.{rep}_2.clean.fastq.gz",
            json = DATA_DIR + "01QC/fastp/{sample}.{treat}.{rep}_fastp.json",
            html = DATA_DIR + "01QC/fastp/{sample}.{treat}.{rep}_fastp.html"
        threads: 1
        params:
            datadir = DATA_DIR,
            mapt = 16,
            adapter1=config["adapter1"],
            adapter2=config["adapter2"]
        shell:
            """
            phred_value=`bash ./workflow/scripts/base_phred.sh {input.reads1}`
            if [[ ${{phred_value}} == 64 ]];then
                if [[ {params.adapter1} == "auto" ]];then
                    fastp -w {params.mapt} --phred64 -i {input.reads1} -o {output.creads1} \
                    -I {input.reads2} -O {output.creads2} \
                    -j {output.json} -h {output.html} -5 -3 -r -n 0 -e 20 --trim_poly_x --poly_x_min_len 15 \
                    --detect_adapter_for_pe \
                    -R "{wildcards.sample}.{wildcards.treat}.{wildcards.rep}" \
                    -U --umi_loc=per_read --umi_len=25
                else
                    fastp -w {params.mapt} --phred64 -i {input.reads1} -o {output.creads1} \
                    -I {input.reads2} -O {output.creads2} \
                    -j {output.json} -h {output.html} -5 -3 -r -n 0 -e 20 --trim_poly_x --poly_x_min_len 15 \
                    --adapter_sequence={params.adapter1} --adapter_sequence_r2={params.adapter2} \
                    -R "{wildcards.sample}.{wildcards.treat}.{wildcards.rep}" \
                    -U --umi_loc=per_read --umi_len=25
                fi
            else
                if [[ {params.adapter1} == "auto" ]];then
                    fastp -w {params.mapt} -i {input.reads1} -o {output.creads1} \
                    -I {input.reads2} -O {output.creads2} \
                    -j {output.json} -h {output.html} -5 -3 -r -n 0 -e 20 --trim_poly_x --poly_x_min_len 15 \
                    --detect_adapter_for_pe \
                    -R "{wildcards.sample}.{wildcards.treat}.{wildcards.rep}" \
                    -U --umi_loc=per_read --umi_len=25
                else
                    fastp -w {params.mapt} -i {input.reads1} -o {output.creads1} \
                    -I {input.reads2} -O {output.creads2} \
                    -j {output.json} -h {output.html} -5 -3 -r -n 0 -e 20 --trim_poly_x --poly_x_min_len 15 \
                    --adapter_sequence={params.adapter1} --adapter_sequence_r2={params.adapter2} \
                    -R "{wildcards.sample}.{wildcards.treat}.{wildcards.rep}" \
                    -U --umi_loc=per_read --umi_len=25
                fi
            fi
            """
else:
    sys.exit()

rule multiQC:
    input: expand(DATA_DIR + "01QC/fastp/{sample}.{treat}.{rep}_fastp.json", sample=config["samples"], rep = config["replicates"], treat = ["input", "ip"])
    output: DATA_DIR + "01QC/multiqc/multiQC.html"
    params:
        datadir = DATA_DIR
    threads: 1
    shell:
        """
        multiqc {params.datadir}01QC/fastp -o {params.datadir}01QC/multiqc -n multiQC
        """

if config["se_or_pe"] == "SE":
    rule strand_analyze:
        input:
            creads = DATA_DIR + "01cleandata/{sample}.{treat}.{rep}.clean.fastq.gz",
            gene_bed = GFF3_FILE + ".gene.bed",
            # miltiqc = DATA_DIR + "01QC/multiqc/multiQC.html"
        output: DATA_DIR + "00rawdata/{sample}.{treat}.{rep}.strand_info.txt"
        params:
            hisat2_index = HISAT_INDEX,
            bamidx = config["bamidxtype"],
            mapt=10
        shell:
            """
            mkdir -p {output}_tmp
            workflow/scripts/extract_fq.sh {input.creads} {output}_tmp/tmp.fastq
            hisat2 -k 20 -p {params.mapt} -x {params.hisat2_index}/Genome -U {output}_tmp/tmp.fastq | \
            samtools view -F 4 -Sb - | samtools sort -T {output}_tmp/tmp -o {output}_tmp/tmp.bam -
            if [[ {params.bamidx} == "large" ]];then
                samtools index -c {output}_tmp/tmp.bam && cp {output}_tmp/tmp.bam.csi {output}_tmp/tmp.bam.bai
            else
                samtools index {output}_tmp/tmp.bam
            fi
            workflow/scripts/SE_infer_experiment.py -i {output}_tmp/tmp.bam -r {input.gene_bed} > {output}
            rm -r {output}_tmp
            """
elif config["se_or_pe"] == "PE":
    rule strand_analyze:
        input:
            creads1 = DATA_DIR + "01cleandata/{sample}.{treat}.{rep}_1.clean.fastq.gz",
            creads2 = DATA_DIR + "01cleandata/{sample}.{treat}.{rep}_2.clean.fastq.gz",
            gene_bed = GFF3_FILE + ".gene.bed",
            # miltiqc = DATA_DIR + "01QC/multiqc/multiQC.html"
        output: DATA_DIR + "00rawdata/{sample}.{treat}.{rep}.strand_info.txt"
        threads: 1
        params:
            hisat2_index = HISAT_INDEX,
            bamidx = config["bamidxtype"],
            mapt=10
        shell:
            """
            mkdir -p {output}_tmp
            ./workflow/scripts/base_extract_fq.sh {input.creads1} {output}_tmp/tmp_1.fastq
            ./workflow/scripts/base_extract_fq.sh {input.creads2} {output}_tmp/tmp_2.fastq
            hisat2 --no-spliced-alignment -k 20 -p {params.mapt} -x {params.hisat2_index}/Genome \
            -1 {output}_tmp/tmp_1.fastq -2 {output}_tmp/tmp_2.fastq | \
            samtools view -F 4 -Sb - | samtools sort -T {output}_tmp/tmp -o {output}_tmp/tmp.bam -
            if [[ {params.bamidx} == "large" ]];then
                samtools index -c {output}_tmp/tmp.bam && cp {output}_tmp/tmp.bam.csi {output}_tmp/tmp.bam.bai
            else
                samtools index {output}_tmp/tmp.bam
            fi
            python ./workflow/scripts/base_infer_strand.py -i {output}_tmp/tmp.bam -r {input.gene_bed} > {output}
            rm -r {output}_tmp
            """
else:
    sys.exit()

if config["se_or_pe"] == "SE":
    rule hisat2_align:
        input: 
            creads = DATA_DIR + "01cleandata/{sample}.{treat}.{rep}.clean.fastq.gz",
            strand = DATA_DIR + "00rawdata/{sample}.input.{rep}.strand_info.txt",
            genome_file = GENOME_FILE
        output:             
            bam = DATA_DIR + "02hisat2/{sample}.{treat}.{rep}.bam",
            bed = DATA_DIR + "02hisat2/{sample}.{treat}.{rep}.bed"
        params:
            datadir = DATA_DIR,
            hisat2_index = HISAT_INDEX,
            minval = 10,
            maxval = config["intron_max"],
            bamidx = config["bamidxtype"],
            mapt = 20,
            umilog = DATA_DIR + "02hisat2/{sample}.{treat}.{rep}.umi.log"
        threads: 1
        shell:
            """
            strand_info=`tail -n 1 {input.strand}`
            TMPOUTDIR={params.datadir}02hisat2 && mkdir -p ${{TMPOUTDIR}}/alignment_log
            XNAME=${{TMPOUTDIR}}/output_{wildcards.sample}.{wildcards.treat}.{wildcards.rep}
            YNAME=${{TMPOUTDIR}}/alignment_log/{wildcards.sample}.{wildcards.treat}.{wildcards.rep}
            if [[ $strand_info != "unstranded" ]];then
                hisat2 -p {params.mapt} -k 5 --max-intronlen {params.maxval} \
                --rna-strandness ${{strand_info}} \
                --summary-file ${{YNAME}}_hisat2.txt \
                -x {params.hisat2_index}/Genome -U {input.creads} -S ${{XNAME}}.sam
            else
                hisat2 -p {params.mapt} -k 5 --max-intronlen {params.maxval} \
                --summary-file ${{YNAME}}_hisat2.txt \
                -x {params.hisat2_index}/Genome -U {input.creads} -S ${{XNAME}}.sam
            fi
            awk '$0~/^@/||$0~/NH:i:1/' ${{XNAME}}.sam | \
            samtools view -b -@ {params.mapt} --reference {input.genome_file} -F 1796 -q 30 -o ${{XNAME}}.bam -
            samtools sort -T ${{XNAME}}_tmp --output-fmt BAM --reference {input.genome_file} -@ {params.mapt} -o {output.bam}_tmp.bam ${{XNAME}}.bam 
            if [[ {params.bamidx} == "large" ]];then
                samtools index -c -@ 10 {output.bam}_tmp.bam && cp {output.bam}_tmp.bam.csi {output.bam}_tmp.bam.bai
            else
                samtools index -@ 10 {output.bam}_tmp.bam
            fi
            bam_stat.py -i {output.bam}_tmp.bam > ${{YNAME}}_unique.txt
            mkdir -p {output.bam}_umitmp
            umi_tools dedup --temp-dir={output.bam}_umitmp --umi-separator=":" --stdin={output.bam}_tmp.bam --log={params.umilog} > {output.bam}
            rm -r {output.bam}_umitmp
            bam_stat.py -i {output.bam} > ${{YNAME}}_umi_re.txt
            bedtools bamtobed -i {output.bam} > {output.bed}
            """
elif config["se_or_pe"] == "PE":
    rule hisat2_align:
        input: 
            creads1 = DATA_DIR + "01cleandata/{sample}.{treat}.{rep}_1.clean.fastq.gz",
            creads2 = DATA_DIR + "01cleandata/{sample}.{treat}.{rep}_2.clean.fastq.gz",
            strand = DATA_DIR + "00rawdata/{sample}.input.{rep}.strand_info.txt",
            genome_file = GENOME_FILE
        output: 
            bam = DATA_DIR + "02hisat2/{sample}.{treat}.{rep}.bam",
            bed = DATA_DIR + "02hisat2/{sample}.{treat}.{rep}.bed"
        params:
            datadir = DATA_DIR,
            hisat2_index = HISAT_INDEX,
            minval=10,
            maxval=config["intron_max"],
            bamidx = config["bamidxtype"],
            mapt=20,
            umilog = DATA_DIR + "02hisat2/{sample}.{treat}.{rep}.umi.log"
        threads: 1
        shell:
            """
            strand_info=`tail -n 1 {input.strand}`
            TMPOUTDIR={params.datadir}02hisat2 && mkdir -p ${{TMPOUTDIR}}/alignment_log
            XNAME=${{TMPOUTDIR}}/output_{wildcards.sample}.{wildcards.treat}.{wildcards.rep}
            YNAME=${{TMPOUTDIR}}/alignment_log/{wildcards.sample}.{wildcards.treat}.{wildcards.rep}
            if [[ $strand_info != "unstranded" ]];then
                hisat2 -p {params.mapt} -k 5 --max-intronlen {params.maxval} \
                --rna-strandness ${{strand_info}} \
                --summary-file ${{YNAME}}_hisat2.txt \
                -x {params.hisat2_index}/Genome -1 {input.creads1} -2 {input.creads2} \
                -S ${{XNAME}}.sam
            else
                hisat2 -p {params.mapt} -k 5 --max-intronlen {params.maxval} \
                --summary-file ${{YNAME}}_hisat2.txt \
                -x {params.hisat2_index}/Genome -1 {input.creads1} -2 {input.creads2} \
                -S ${{XNAME}}.sam
            fi
            samtools view -@ 10 -bS -o ${{XNAME}}.raw.bam ${{XNAME}}.sam
            awk '$0~/^@/||$0~/NH:i:1/' ${{XNAME}}.sam | \
            samtools view -b -@ {params.mapt} --reference {input.genome_file} -f 3 -F 1804 -q 30 -o ${{XNAME}}.bam -            
            samtools sort -T ${{XNAME}}_tmp --output-fmt BAM --reference {input.genome_file} -@ {params.mapt} -o {output.bam}_tmp.bam ${{XNAME}}.bam
            if [[ {params.bamidx} == "large" ]];then
                samtools index -c -@ 10 {output.bam}_tmp.bam && cp {output.bam}_tmp.bam.csi {output.bam}_tmp.bam.bai
            else
                samtools index -@ 10 {output.bam}_tmp.bam
            fi
            bam_stat.py -i {output.bam}_tmp.bam > ${{YNAME}}_unique.txt
            mkdir -p {output.bam}_umitmp
            umi_tools dedup --temp-dir={output.bam}_umitmp --paired --umi-separator=":" --stdin={output.bam}_tmp.bam --log={params.umilog} > {output.bam}
            rm -r {output.bam}_umitmp
            bam_stat.py -i {output.bam} > ${{YNAME}}_umi_re.txt
            samtools sort -T ${{XNAME}}_nsort --output-fmt BAM -@ {params.mapt} -n -o ${{XNAME}}_nsort.bam {output.bam}
            bedtools bamtobed -bedpe -mate1 -i ${{XNAME}}_nsort.bam > ${{XNAME}}_tmp.bed
            if [[ $strand_info == "RF" ]];then
                awk '{{a=i++;if($2>=$5){{min=$5}}else{{min=$2}};
                if($3>$6){{max=$3}}else{{max=$6}};
                OFS="\t";print $1,min,max,"read_"a+1,".", $10}}' ${{XNAME}}_tmp.bed | \
                sort -k 1,1 -k 2,2n - > {output.bed}
            else
                awk '{{a=i++;if($2>=$5){{min=$5}}else{{min=$2}};
                if($3>$6){{max=$3}}else{{max=$6}};
                OFS="\t";print $1,min,max,"read_"a+1,".", $9}}' ${{XNAME}}_tmp.bed | \
                sort -k 1,1 -k 2,2n - > {output.bed}
            fi
            if [[ $strand_info != "unstranded" ]];then
                bash workflow/scripts/bam_split.sh {output.bam} $strand_info {params.mapt} {params.bamidx}
            fi
            rm ${{XNAME}}.sam ${{XNAME}}.bam {output.bam}_tmp.bam {output.bam}_tmp.bam.* ${{XNAME}}_nsort.bam ${{XNAME}}_tmp.bed
            """
else:
    sys.exit()

# samtools view -@ {params.mapt} -bS -o ${{XNAME}}.raw.bam ${{XNAME}}.sam
# rm {output.bam}_nsort.bam {output.bam}_tmp.bed ${{XNAME}}_tmp.bam ${{XNAME}}.sam

rule read_density:
    input:
        IP = expand(DATA_DIR + "02hisat2/{sample}.ip.{rep}.bam", sample=config["samples"], rep = config["replicates"]),
        Input =  expand(DATA_DIR + "02hisat2/{sample}.input.{rep}.bam", sample=config["samples"], rep = config["replicates"]),
        gtf = '{0}/{1}/Annotation/{1}.exons_tmp.gtf'.format(SPECIES_DIR, config["species_name"])
    output: DATA_DIR + "final_output/read_density.png"
    threads: 1
    shell:
        """
        set -x
        /usr/bin/Rscript workflow/scripts/read_distribution_using_trumpet.R {input.gtf} {input.IP} {input.Input} {output}
        """

rule peak_calling:
    input:
        IP = DATA_DIR + "02hisat2/{sample}.ip.{rep}.bed",
        Input = DATA_DIR + "02hisat2/{sample}.input.{rep}.bed",
        strand = DATA_DIR + "00rawdata/{sample}.input.{rep}.strand_info.txt",
        genome = GENOME_SIZE,
        read_density = DATA_DIR + "final_output/read_density.png"
    output: DATA_DIR + "03calling/{sample}.{rep}.bed"
    threads: 1
    shell:
        """
        strand_info=`tail -n 1 {input.strand}`
        if [[ $strand_info != "unstranded" ]] ;then
            awk '$6=="+"' {input.IP} > {input.IP}.fwd.bed
            awk '$6=="-"' {input.IP} > {input.IP}.rev.bed
            awk '$6=="+"' {input.Input} > {input.Input}.fwd.bed
            awk '$6=="-"' {input.Input} > {input.Input}.rev.bed
            ip_fwd=`wc -l {input.IP}.fwd.bed | cut -f 1 -d " "`
            ip_rev=`wc -l {input.IP}.rev.bed | cut -f 1 -d " "`
            input_fwd=`wc -l {input.Input}.fwd.bed | cut -f 1 -d " "`
            input_rev=`wc -l {input.Input}.rev.bed | cut -f 1 -d " "`
            ip_count=`expr ${{ip_fwd}} + ${{ip_rev}}`
            input_count=`expr ${{input_fwd}} + ${{input_rev}}`

            for kk in fwd rev
            do
                Rscript workflow/scripts/PEA_peak_calling.R \
                -i {input.IP}.${{kk}}.bed \
                -n {input.Input}.${{kk}}.bed \
                -g {input.genome} \
                -p ${{ip_count}} -q ${{input_count}} \
                -o {output}.${{kk}}.txt
                sed "s/ //g" {output}.${{kk}}.txt > {output}.${{kk}}_tmp.txt
                awk -F"\t" '{{if($2<=0){{$2=1}};OFS="\t";print $0}}' {output}.${{kk}}_tmp.txt >{output}.${{kk}}.txt
            done
            awk -F"\t" '{{OFS="\t";print $1,$2-1,$3, ".", ".", "+"}}' {output}.fwd.txt > {output}.fwd.bed
            awk -F"\t" '{{OFS="\t";print $1,$2-1,$3, ".", ".", "-"}}' {output}.rev.txt > {output}.rev.bed
            cat {output}.fwd.bed {output}.rev.bed | sort -k 1,1 -k 2,2n > {output}
        else
            ip_count=`wc -l {input.IP} | cut -f 1 -d " "`
            input_count=`wc -l {input.Input} | cut -f 1 -d " "`
            Rscript workflow/scripts/PEA_peak_calling.R \
            -i {input.IP} \
            -n {input.Input} \
            -g {input.genome} \
            -p ${{ip_count}} -q ${{input_count}} \
            -o {output}.txt
            sed "s/ //g" {output}.txt > {output}_tmp.txt
            awk -F"\t" '{{if($2<=0){{$2=1}};OFS="\t";print $0}}' {output}_tmp.txt >{output}.txt
            awk -F"\t" '{{OFS="\t";print $1,$2-1,$3, ".", ".", "."}}' {output}.txt > {output}
        fi
        """

rule stringtie:
    input:
        bam = DATA_DIR + "02hisat2/{sample}.input.{rep}.bam",
        gtf = GTF_FILE,
        strand = DATA_DIR + "00rawdata/{sample}.input.{rep}.strand_info.txt"
    output:
        gtf_out = DATA_DIR + "04gene/{sample}.{rep}.gtf",
        gene_exp = DATA_DIR + "04gene/{sample}.{rep}.exp.txt",
        count_out = DATA_DIR + "04gene/{sample}.{rep}.count.txt"
    params:
        tinshell = 20
    threads: 2
    shell: 
        """
        strand_info=`tail -n 1 {input.strand}`
        if [[ ${{strand_info}} == "RF" ]] || [[ ${{strand_info}} == "R" ]];then
            stringtie -e -p {params.tinshell} --rf -G {input.gtf} -j 10 -o {output.gtf_out} \
            -A {output.gene_exp} {input.bam}
        elif [[ ${{strand_info}} == "FR" ]] || [[ ${{strand_info}} == "F" ]];then
            stringtie -e -p {params.tinshell} --fr -G {input.gtf} -j 10 -o {output.gtf_out} \
            -A {output.gene_exp} {input.bam}
        else
            stringtie -e -p {params.tinshell} -G {input.gtf} -j 10 -o {output.gtf_out} \
            -A {output.gene_exp} {input.bam}
        fi
        featureCounts -T {params.tinshell} -p -t exon -g gene_id -a {input.gtf} -o {output.count_out} {input.bam}
        """

rule add_gene_strand:
    input:
        peak = rules.peak_calling.output,
        expGene = rules.stringtie.output.gene_exp,
        gene = GFF3_FILE + ".gene.bed",
        strand = DATA_DIR + "00rawdata/{sample}.input.{rep}.strand_info.txt"
    output: DATA_DIR + "05filter/{sample}.{rep}.bed"
    threads: 1
    shell:
        """
        sed -i "s/gene://" {input.expGene}
        strand_info=`tail -n 1 {input.strand}`
        if [[ $strand_info != "unstranded" ]];then
            awk -F"\t" 'NR>1&&$9>=1{{OFS="\t";print $3,$5-1,$6,$1,$9,$4}}' {input.expGene} | \
            awk -F"\t" 'NR==FNR{{a[$4]=$4;next}}{{if(a[$4]){{print $0}}}}' {input.gene} - | \
            bedtools intersect -nonamecheck -wo -s -a {input.peak}.fwd.bed -b - | \
            awk -F"\t" '$13>=100{{OFS="\t";print $1,$2,$3,$10,$11,$12}}' -  > {output}.fwd.bed
            awk -F"\t" 'NR>1&&$9>=1{{OFS="\t";print $3,$5-1,$6,$1,$9,$4}}' {input.expGene} | \
            awk -F"\t" 'NR==FNR{{a[$4]=$4;next}}{{if(a[$4]){{print $0}}}}' {input.gene} - | \
            bedtools intersect -nonamecheck -wo -s -a {input.peak}.rev.bed -b - | \
            awk -F"\t" '$13>=100{{OFS="\t";print $1,$2,$3,$10,$11,$12}}' -  > {output}.rev.bed
            cat {output}.fwd.bed {output}.rev.bed >{output}
        else
            awk -F"\t" 'NR>1&&$9>=1{{OFS="\t";print $3,$5-1,$6,$1,$9,$4}}' {input.expGene} | \
            awk -F"\t" 'NR==FNR{{a[$4]=$4;next}}{{if(a[$4]){{print $0}}}}' {input.gene} - | \
            bedtools intersect -nonamecheck -wo -a {input.peak} -b - | \
            awk -F"\t" '$13>=100{{OFS="\t";print $1,$2,$3,$10,$11,$12}}' -  > {output}
        fi
        """

rule rep_merge:
    input:
        peak = expand(DATA_DIR + "05filter/{sample}.{rep}.bed", sample=config["samples"], rep = config["replicates"]),
        strand = DATA_DIR + "00rawdata/{sample}.input.rep1.strand_info.txt"
    output: DATA_DIR + "05filter/{sample}.bed"
    params:
        SAPNUM = len(config['replicates']),
        datadir = DATA_DIR + "05filter/{sample}"
    shell:
        """
        strand_info=`tail -n 1 {input.strand}`
        peakfix={params.datadir}
        if [[ $strand_info != "unstranded" ]];then
            for kk in fwd rev
            do
                if [[ {params.SAPNUM} == 3 ]];then
                    /usr/bin/Rscript workflow/scripts/samples_peak_merge.R ${{peakfix}}.rep1.bed.${{kk}}.bed \
                    ${{peakfix}}.rep2.bed.${{kk}}.bed ${{peakfix}}.rep3.bed.${{kk}}.bed \
                    ${{peakfix}}.bed.${{kk}}_merge.bed
                 elif [[ {params.SAPNUM} == 2 ]];then
                    /usr/bin/Rscript workflow/scripts/samples_peak_merge.R ${{peakfix}}.rep1.bed.${{kk}}.bed \
                    ${{peakfix}}.rep2.bed.${{kk}}.bed ${{peakfix}}.bed.${{kk}}_merge.bed
                fi
                if [[ ${{kk}} == "fwd" ]];then
                    awk -F"\t" '{{OFS="\t";$6="+";print $0}}' ${{peakfix}}.bed.fwd_merge.bed >{output}.fwd_merge_mod.bed
                else
                    awk -F"\t" '{{OFS="\t";$6="-";print $0}}' ${{peakfix}}.bed.rev_merge.bed >{output}.rev_merge_mod.bed
                fi
            done
            cat {output}.fwd_merge_mod.bed {output}.rev_merge_mod.bed > {output}
        else
            if [[ {params.SAPNUM} == 3 ]];then
                /usr/bin/Rscript workflow/scripts/samples_peak_merge.R ${{peakfix}}.rep1.bed \
                ${{peakfix}}.rep2.bed ${{peakfix}}.rep3.bed {output}
            elif [[ {params.SAPNUM} == 2 ]];then
                /usr/bin/Rscript workflow/scripts/samples_peak_merge.R ${{peakfix}}.rep1.bed ${{peakfix}}.rep2.bed {output}
            fi
        fi
        """

rule readd_strand:
    input:
        peak = DATA_DIR + "05filter/{sample}.bed",
        gene = GFF3_FILE + ".gene.bed",
        strand = DATA_DIR + "00rawdata/{sample}.input.rep1.strand_info.txt"
    output: DATA_DIR + "06peak/{sample}.bed"
    params:
        SAPNUM = len(config['replicates']),
        datadir = DATA_DIR + "05filter/{sample}"
    shell:
        """
        strand_info=`tail -n 1 {input.strand}`
        peakfix={params.datadir}
        if [[ $strand_info != "unstranded" ]];then
            for kk in fwd rev
            do
                if [[ {params.SAPNUM} == 3 ]];then
                    cat ${{peakfix}}.rep1.bed.${{kk}}.bed ${{peakfix}}.rep2.bed.${{kk}}.bed \
                    ${{peakfix}}.rep3.bed.${{kk}}.bed > ${{peakfix}}.${{kk}}.bed
                elif [[ {params.SAPNUM} == 2 ]];then
                    cat ${{peakfix}}.rep1.bed.${{kk}}.bed ${{peakfix}}.rep2.bed.${{kk}}.bed > ${{peakfix}}.${{kk}}.bed
                fi
                cut -f 4 ${{peakfix}}.${{kk}}.bed | sort -u | \
                awk -F"\t" 'NR==FNR{{a[$1]=$1;next}}{{if(a[$4]){{OFS="\t";print $0}}}}' - {input.gene} |
                bedtools intersect -nonamecheck -wo -s -a {input.peak}.${{kk}}_merge_mod.bed -b - | \
                awk -F"\t" '$13>=100{{OFS="\t";print $1,$2,$3,$10,$11,$12}}' - | sort -k1,1 -k2,2n - \
                > {output}.${{kk}}.bed
            done
            cat {output}.fwd.bed {output}.rev.bed > {output}
        else
            if [[ {params.SAPNUM} == 3 ]];then
                cat ${{peakfix}}.rep1.bed ${{peakfix}}.rep2.bed \
                ${{peakfix}}.rep3.bed > ${{peakfix}}
            elif [[ {params.SAPNUM} == 2 ]];then
                cat ${{peakfix}}.rep1.bed ${{peakfix}}.rep2.bed > ${{peakfix}}
            fi
            cut -f 4 ${{peakfix}} | sort -u | \
            awk -F"\t" 'NR==FNR{{a[$1]=$1;next}}{{if(a[$4]){{OFS="\t";print $0}}}}' - {input.gene} |
            bedtools intersect -nonamecheck -wo -a {input.peak} -b - | \
            awk -F"\t" '$13>=100{{OFS="\t";print $1,$2,$3,$10,$11,$12}}' - | sort -k1,1 -k2,2n - > {output}
        fi
        """

rule diffbind_prepare:
    input:
        peak = DATA_DIR + "06peak/{sample}.bed",
        bam = expand(DATA_DIR + "02hisat2/{sample}.{treat}.{rep}.bam", sample=config["samples"], rep = config["replicates"], treat = ["input", "ip"]),
        strand = DATA_DIR + "00rawdata/{sample}.input.rep1.strand_info.txt"
    output: 
        difftable = DATA_DIR + "07diffbind/{sample}.DiffBind_conf.csv",
        diffout = DATA_DIR + "07diffbind/{sample}_diffbind.txt"
    params:
        SAPNUM = len(config['replicates']),
        datadir = DATA_DIR + "02hisat2/{sample}"
    shell:
        """
        peakfix={params.datadir}
        strand_info=`tail -n 1 {input.strand}`
        if [[ $strand_info != "unstranded" ]];then
            for kk in fwd rev
            do
                cut -f 1-3 {input.peak}.${{kk}}.bed > {input.peak}_cut.bed.${{kk}}.bed 
                if [[ {params.SAPNUM} == 3 ]];then
echo -e "SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
Sample1,Leaf,Leaf1,IP1,1,${{peakfix}}.ip.rep1.bam.${{kk}}.bam,Input1,${{peakfix}}.input.rep1.bam.${{kk}}.bam,{input.peak}_cut.bed.${{kk}}.bed,bed
Sample2,Leaf,Leaf2,IP2,2,${{peakfix}}.ip.rep2.bam.${{kk}}.bam,Input2,${{peakfix}}.input.rep2.bam.${{kk}}.bam,{input.peak}_cut.bed.${{kk}}.bed,bed
Sample3,Leaf,Leaf3,IP3,3,${{peakfix}}.ip.rep3.bam.${{kk}}.bam,Input3,${{peakfix}}.input.rep3.bam.${{kk}}.bam,{input.peak}_cut.bed.${{kk}}.bed,bed
" > {output.difftable}.${{kk}}.txt
                elif [[ {params.SAPNUM} == 2 ]];then
echo -e "SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
Sample1,Leaf,Leaf1,IP1,1,${{peakfix}}.ip.rep1.bam.${{kk}}.bam,Input1,${{peakfix}}.input.rep1.bam.${{kk}}.bam,{input.peak}_cut.bed.${{kk}}.bed,bed
Sample2,Leaf,Leaf2,IP2,2,${{peakfix}}.ip.rep2.bam.${{kk}}.bam,Input2,${{peakfix}}.input.rep2.bam.${{kk}}.bam,{input.peak}_cut.bed.${{kk}}.bed,bed
" > {output.difftable}.${{kk}}.txt
                else
echo "SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
Sample1,Leaf,Leaf1,IP1,1,${{peakfix}}.ip.rep1.bam.${{kk}}.bam,Input1,${{peakfix}}.input.rep1.bam.${{kk}}.bam,{input.peak}_cut.bed.${{kk}}.bed,bed" > {output.difftable}.${{kk}}.txt
                fi
                /usr/bin/Rscript workflow/scripts/diffbind2.R {output.difftable}.${{kk}}.txt {input.peak}_cut.bed.${{kk}}.bed {output.diffout}.${{kk}}.txt
            done
            echo "Finished" >{output.difftable}
            cat {output.diffout}.fwd.txt {output.diffout}.rev.txt >{output.diffout}
        else
            cut -f 1-3 {input.peak} > {input.peak}_cut.bed
            if [[ {params.SAPNUM} == 3 ]];then
echo -e "SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
Sample1,Leaf,Leaf1,IP1,1,${{peakfix}}.ip.rep1.bam,Input1,${{peakfix}}.input.rep1.bam,{input.peak}_cut.bed,bed
Sample2,Leaf,Leaf2,IP2,2,${{peakfix}}.ip.rep2.bam,Input2,${{peakfix}}.input.rep2.bam,{input.peak}_cut.bed,bed
Sample3,Leaf,Leaf3,IP3,3,${{peakfix}}.ip.rep3.bam,Input3,${{peakfix}}.input.rep3.bam,{input.peak}_cut.bed,bed
" > {output.difftable}
            elif [[ {params.SAPNUM} == 2 ]];then
echo -e "SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
Sample1,Leaf,Leaf1,IP1,1,${{peakfix}}.ip.rep1.bam,Input1,${{peakfix}}.input.rep1.bam,{input.peak}_cut.bed,bed
Sample2,Leaf,Leaf2,IP2,2,${{peakfix}}.ip.rep2.bam,Input2,${{peakfix}}.input.rep2.bam,{input.peak}_cut.bed,bed
" > {output.difftable}
            else
echo "SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
Sample1,Leaf,Leaf1,IP1,1,${{peakfix}}.ip.rep1.bam,Input1,${{peakfix}}.input.rep1.bam,{input.peak}_cut.bed,bed" > {output.difftable}
            fi
            /usr/bin/Rscript workflow/scripts/diffbind2.R {output.difftable} {input.peak}_cut.bed {output.diffout}
        fi
        """

rule diffbind_peak:
    input:
        peak = rules.diffbind_prepare.output.diffout,
        strand = DATA_DIR + "00rawdata/{sample}.input.rep1.strand_info.txt"
    output: DATA_DIR + "07diffbind/{sample}_gene_ratio.txt"
    params:
        SAPNUM = len(config['replicates']),
        datadir = DATA_DIR + "06peak/{sample}"
    shell:
        """
        strand_info=`tail -n 1 {input.strand}`
        peakfix={params.datadir}
        if [[ $strand_info != "unstranded" ]];then
            for kk in fwd rev
            do
                if [[ {params.SAPNUM} == 3 ]];then
                    awk -F"\t" 'NR==FNR{{a[$1"_"$2"_"$3]=$0;next}}{{if(a[$1"_"$2"_"$3]){{print($0"\t"a[$1"_"$2"_"$3])}}}}' \
                    {input.peak}.${{kk}}.txt ${{peakfix}}.bed.${{kk}}.bed | \
                    cut -f 1-6,10-21 | awk '{{OFS="\t";$2+=1;print $0}}' - >{output}.${{kk}}.txt
                elif [[ {params.SAPNUM} == 2 ]];then
                    awk -F"\t" 'NR==FNR{{a[$1"_"$2"_"$3]=$0;next}}{{if(a[$1"_"$2"_"$3]){{print($0"\t"a[$1"_"$2"_"$3])}}}}' \
                    {input.peak}.${{kk}}.txt ${{peakfix}}.bed.${{kk}}.bed | \
                    cut -f 1-6,10-17 | awk '{{OFS="\t";$2+=1;print $0}}' - >{output}.${{kk}}.txt
                fi
            done
            cat {output}.fwd.txt {output}.rev.txt > {output}
        else
            if [[ {params.SAPNUM} == 3 ]];then
                awk -F"\t" 'NR==FNR{{a[$1"_"$2"_"$3]=$0;next}}{{if(a[$1"_"$2"_"$3]){{print($0"\t"a[$1"_"$2"_"$3])}}}}' \
                {input.peak} ${{peakfix}}.bed | \
                cut -f 1-6,10-21 | awk '{{OFS="\t";$2+=1;print $0}}' - >{output}
            elif [[ {params.SAPNUM} == 2 ]];then
                awk -F"\t" 'NR==FNR{{a[$1"_"$2"_"$3]=$0;next}}{{if(a[$1"_"$2"_"$3]){{print($0"\t"a[$1"_"$2"_"$3])}}}}' \
                {input.peak} ${{peakfix}}.bed | \
                cut -f 1-6,10-17 | awk '{{OFS="\t";$2+=1;print $0}}' - >{output}
            fi
        fi
        """

rule get_fasta:
    input:
        bed = DATA_DIR + "07diffbind/{sample}_gene_ratio.txt",
        genome = GENOME_FILE
    output: DATA_DIR + "08fasta/{sample}_peaks_T.fa"
    threads: 1
    shell:
        """
        set -x
        tmpv=`dirname {output}`/{wildcards.sample}
        cut -f 1-6 {input.bed} | bedtools getfasta -fi {input.genome} -bed - -fo ${{tmpv}}.fa -s
        awk '{{if($1!~/>/){{$1=toupper($1)}};print $0}}' ${{tmpv}}.fa >{output}
        /usr/bin/Rscript workflow/scripts/m6A_motif_ratio.R {output} >${{tmpv}}_motif_ratio.txt
        bedtools shuffle -i {input.bed} -g {input.genome}.fai >${{tmpv}}_shuffle.bed
        bedtools getfasta -fi {input.genome} -bed ${{tmpv}}_shuffle.bed -fo ${{tmpv}}_shuffle.fa -s
        awk '{{if($1!~/>/){{$1=toupper($1)}};print $0}}' ${{tmpv}}_shuffle.fa >${{tmpv}}_shuffle_T.fa
        /usr/bin/Rscript workflow/scripts/m6A_motif_ratio.R ${{tmpv}}_shuffle_T.fa >${{tmpv}}_motif_ratio_shuffle.txt
        findMotifs.pl {output} fasta ${{tmpv}}_homer_motif -fasta ${{tmpv}}_shuffle_T.fa -len 5,6 -rna -p 20 -S 10
        findMotifs.pl {output} fasta ${{tmpv}}_RRACH -fasta ${{tmpv}}_shuffle_T.fa -norevopp -nomotif -mknown workflow/scripts/homer/Homer_36_RRACH_merge.motif
        findMotifs.pl {output} fasta ${{tmpv}}_TRTAY -fasta ${{tmpv}}_shuffle_T.fa -norevopp -nomotif -mknown workflow/scripts/homer/Homer_9_TRTAY_merge.motif
        """

rule bam_bw:
    input: DATA_DIR + "02hisat2/{sample}.{treat}.{rep}.bam"
    output: DATA_DIR + "final_output/{sample}.{treat}.{rep}.bw"
    params:
        bamidx = config["bamidxtype"]
    threads: 1
    shell:
        """
        if [[ {params.bamidx} == "large" ]];then
            samtools index -c -@ 10 {input} && cp {input}.csi {input}.bai
        else
            samtools index -@ 10 {input}
        fi
        bamCoverage -b {input} -o {output} --normalizeUsing RPKM -bs 20 -p 10
        """

rule corrplot:
    input: 
        difftable = DATA_DIR + "07diffbind/{sample}.DiffBind_conf.csv",
        strand = DATA_DIR + "00rawdata/{sample}.input.rep1.strand_info.txt"
    output: DATA_DIR + "07diffbind/{sample}.corrplot.png"
    params:
        data_dir = DATA_DIR,
        SAPNUM = len(config['replicates'])
    threads: 1
    shell: 
        """
        /usr/bin/Rscript workflow/scripts/m6A_relicate_corrplot.R {params.SAPNUM} {input.difftable} {output} `tail -n 1 {input.strand}`
        cp {output} {params.data_dir}/final_output/corrplot.png
        """