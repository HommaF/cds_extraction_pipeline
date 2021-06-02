configfile: "config.json"

rule all:
	input:
		expand("output_cds/{genome}_{sample}/one_miss", genome=config["ref_genome"], sample=config["vcf_phased"])

rule phase_homozoygous:
	input:
		"vcf_phased/{genome}_{sample}_phased.bcf"
	output:
		first_nine="vcf_phased/{genome}_{sample}_first_nine.vcf",
		headers="vcf_phased/{genome}_{sample}_headers.vcf",
		tenth="vcf_phased/{genome}_{sample}_tenth.vcf",
		vcf_phased="vcf_phased/{genome}_{sample}_hz_phased.bcf"
	threads: 1
	run:
		shell("bcftools view -Ov {input} | grep $'^SL' | cut -f1,2,3,4,5,6,7,8,9 > {output.first_nine}")
		shell("bcftools view -Ov {input} | grep $'^SL' | cut -f10 | sed 's/^1\/1/1\|1/g' > {output.tenth}")
		shell("bcftools view -Ov {input} | grep $'#' > {output.headers}")
		shell("paste {output.first_nine} {output.tenth} >> {output.headers}")
		shell("bcftools view -Ob {output.headers} > {output.vcf_phased}")
		shell("bcftools index {output.vcf_phased}")

rule identify_phased_genes:
	input:
		work="vcf_phased/{genome}_{sample}_hz_phased.bcf",
		rm1="vcf_phased/{genome}_{sample}_first_nine.vcf",
		rm2="vcf_phased/{genome}_{sample}_headers.vcf",
		rm3="vcf_phased/{genome}_{sample}_tenth.vcf"
	output:
		gen_phased="candidates/{genome}_{sample}_fully_phased.txt",
		single_block="candidates/{genome}_{sample}_single_block_phased.txt",
		one_miss="candidates/{genome}_{sample}_one_phased_missing.txt",
		no_var="candidates/{genome}_{sample}_no_variants.txt"
	params:
		names=expand("ref_genome/{genes_of_interest}", genes_of_interest=config["file_with_goi"]),
		cds=expand("ref_genome/{cds_reference}.gff", cds_reference=config["ref_genome_annot"]),
		sample="{sample}"

	threads: 2
	conda:
		"envs/cds_extraction.yaml"
	shell:
		"cat {params.names} | xargs -I@ -P 2 sh -c 'scripts/phased_genes.sh @ {params.cds} {input.work} {output.gen_phased} {output.single_block} {output.one_miss} {params.sample} {output.no_var}'; "
		"rm {input.rm1}; rm {input.rm2}; rm {input.rm3}"

rule locus_bed:
	input:
		fully_phased="candidates/{genome}_{sample}_fully_phased.txt",
		one_miss="candidates/{genome}_{sample}_one_phased_missing.txt"
	output:
		fully_phased = directory("tmp/{genome}_{sample}/fully_phased_locus_bed/"),
		one_miss = directory("tmp/{genome}_{sample}/one_miss_locus_bed/")
	params:
		mrna="ref_genome/ITAG4.0_mrna_models.gff"
	threads: 2
	shell:
		"mkdir {output.fully_phased}; "
		"cat {input.fully_phased} | xargs -I@ -P 2 sh -c 'scripts/awk_gene.sh @ {params.mrna} {output.fully_phased}'; "

		"mkdir {output.one_miss}; "
		"cat {input.one_miss} | xargs -I@ -P 2 sh -c 'scripts/awk_gene.sh @ {params.mrna} {output.one_miss}'"


rule phased_cds_bed:
	input:
		fully_phased="candidates/{genome}_{sample}_fully_phased.txt",
		one_miss="candidates/{genome}_{sample}_one_phased_missing.txt"
	output:
		fully_phased="vcf_phased/{genome}_{sample}_hz_fully_phased.bed",
		one_miss="vcf_phased/{genome}_{sample}_hz_one_missing.bed"
	params:
		cds=expand("ref_genome/{cds_reference}.gff", cds_reference=config["ref_genome_annot"]),
	threads: 1
	conda:
		"envs/cds_extraction.yaml"
	shell:
		"cat {input.fully_phased} | xargs -I@ sh -c 'grep @ {params.cds} | cut -f1,4,5 >> {output.fully_phased}'; "
		"cat {input.one_miss} | xargs -I@ sh -c 'grep @ {params.cds} | cut -f1,4,5 >> {output.one_miss}'"


rule phased_cds_intersect:
	input:
		bed_fully_phased="vcf_phased/{genome}_{sample}_hz_fully_phased.bed",
		bed_one_miss="vcf_phased/{genome}_{sample}_hz_one_missing.bed",
		bcf="vcf_phased/{genome}_{sample}_hz_phased.bcf"
	output:
		vcf_fully_phased="vcf_phased/{genome}_{sample}_hz_fully_phased.vcf",
		vcf_gz_fully_phased="vcf_phased/{genome}_{sample}_hz_fully_phased.vcf.gz",
		vcf_one_miss="vcf_phased/{genome}_{sample}_hz_one_missing.vcf",
		vcf_gz_one_miss="vcf_phased/{genome}_{sample}_hz_one_missing.vcf.gz"
	threads: 1
	run:
		shell("bcftools view -Ov {input.bcf} | grep $'#' - > {output.vcf_fully_phased}")
		shell("bcftools view -Oz {input.bcf} | bedtools intersect -a stdin -b {input.bed_fully_phased} >> {output.vcf_fully_phased}")
		shell("bcftools view -Oz {output.vcf_fully_phased} > {output.vcf_gz_fully_phased}")
		shell("bcftools index {output.vcf_gz_fully_phased}")

		shell("bcftools view -Ov {input.bcf} | grep $'#' - > {output.vcf_one_miss}")
		shell("bcftools view -Oz {input.bcf} | bedtools intersect -a stdin -b {input.bed_one_miss} >> {output.vcf_one_miss}")
		shell("bcftools view -Oz {output.vcf_one_miss} > {output.vcf_gz_one_miss}")
		shell("bcftools index {output.vcf_gz_one_miss}")

rule phased_genomes:
	input:
		vcf_gz_fully_phased="vcf_phased/{genome}_{sample}_hz_fully_phased.vcf.gz",
		vcf_gz_one_miss="vcf_phased/{genome}_{sample}_hz_one_missing.vcf.gz"

	output:
		h1_fully_phased="tmp/{genome}_{sample}/fully_phased_h1.fasta",
		h2_fully_phased="tmp/{genome}_{sample}/fully_phased_h2.fasta",
		
		h1_one_miss="tmp/{genome}_{sample}/one_miss_h1.fasta",
		h2_one_miss="tmp/{genome}_{sample}/one_miss_h2.fasta"

	params:
		h1=1,
		h2=2,
		ref_genome=expand("ref_genome/{genome}.fa", genome=config["ref_genome"])

	threads: 1
	conda:
		"envs/cds_extraction.yaml"
	shell:
		"cat {params.ref_genome} | bcftools consensus -H {params.h1} {input.vcf_gz_fully_phased} > {output.h1_fully_phased};"
		"cat {params.ref_genome} | bcftools consensus -H {params.h2} {input.vcf_gz_fully_phased} > {output.h2_fully_phased};"

		"cat {params.ref_genome} | bcftools consensus -H {params.h1} {input.vcf_gz_one_miss} > {output.h1_one_miss};"
		"cat {params.ref_genome} | bcftools consensus -H {params.h2} {input.vcf_gz_one_miss} > {output.h2_one_miss}"


rule index_genomes:
	input:
		h1_fully_phased="tmp/{genome}_{sample}/fully_phased_h1.fasta",
		h2_fully_phased="tmp/{genome}_{sample}/fully_phased_h2.fasta",
		h1_one_miss="tmp/{genome}_{sample}/one_miss_h1.fasta",
		h2_one_miss="tmp/{genome}_{sample}/one_miss_h2.fasta"

	output:
		h1_fully_phased="tmp/{genome}_{sample}/fully_phased_h1.fasta.fai",
		h2_fully_phased="tmp/{genome}_{sample}/fully_phased_h2.fasta.fai",
		h1_one_miss="tmp/{genome}_{sample}/one_miss_h1.fasta.fai",
		h2_one_miss="tmp/{genome}_{sample}/one_miss_h2.fasta.fai"

	threads: 1

	conda:
		"envs/cds_extraction.yaml"
	shell:
		"samtools faidx {input.h1_fully_phased}; touch {output.h1_fully_phased}; "
		"samtools faidx {input.h2_fully_phased}; touch {output.h2_fully_phased}; "
		"samtools faidx {input.h1_one_miss}; touch {output.h1_one_miss};"
		"samtools faidx {input.h2_one_miss}; touch {output.h2_one_miss}"


rule mrna_phased:
	input:
		fully_phased="candidates/{genome}_{sample}_fully_phased.txt",
		locus_fully_phased="tmp/{genome}_{sample}/fully_phased_locus_bed",
		genome_h1_fully_phased="tmp/{genome}_{sample}/fully_phased_h1.fasta",
		genome_h2_fully_phased="tmp/{genome}_{sample}/fully_phased_h2.fasta",
		index_h1_fully_phased="tmp/{genome}_{sample}/fully_phased_h1.fasta.fai",
		index_h2_fully_phased="tmp/{genome}_{sample}/fully_phased_h2.fasta.fai",

		one_miss="candidates/{genome}_{sample}_one_phased_missing.txt",
		locus_one_miss="tmp/{genome}_{sample}/one_miss_locus_bed",
		genome_h1_one_miss="tmp/{genome}_{sample}/one_miss_h1.fasta",
		genome_h2_one_miss="tmp/{genome}_{sample}/one_miss_h2.fasta",
		index_h1_one_miss="tmp/{genome}_{sample}/one_miss_h1.fasta.fai",
		index_h2_one_miss="tmp/{genome}_{sample}/one_miss_h2.fasta.fai"

	output:
		h1_fully_phased=directory("tmp/{genome}_{sample}/fully_phased_locus_fasta_h1"),
		h2_fully_phased=directory("tmp/{genome}_{sample}/fully_phased_locus_fasta_h2"),

		h1_one_miss=directory("tmp/{genome}_{sample}/one_miss_locus_fasta_h1"),
		h2_one_miss=directory("tmp/{genome}_{sample}/one_miss_locus_fasta_h2")

	params:	
		h1=1,
		h2=2

	threads: 2

        run:
		shell("mkdir {output.h1_fully_phased}; mkdir {output.h2_fully_phased}")
		shell("cat {input.fully_phased} | xargs -I@ -P 2 sh -c 'scripts/mrna_phased.sh @ {input.genome_h1_fully_phased} {input.locus_fully_phased} {output.h1_fully_phased} {params.h1}'")
		shell("cat {input.fully_phased} | xargs -I@ -P 2 sh -c 'scripts/mrna_phased.sh @ {input.genome_h2_fully_phased} {input.locus_fully_phased} {output.h2_fully_phased} {params.h2}'")

		shell("mkdir {output.h1_one_miss}; mkdir {output.h2_one_miss}")
		shell("cat {input.one_miss} | xargs -I@ -P 2 sh -c 'scripts/mrna_phased.sh @ {input.genome_h1_one_miss} {input.locus_one_miss} {output.h1_one_miss} {params.h1}'")
		shell("cat {input.one_miss} | xargs -I@ -P 2 sh -c 'scripts/mrna_phased.sh @ {input.genome_h2_one_miss} {input.locus_one_miss} {output.h2_one_miss} {params.h2}'")

rule exonerate_fully_phased:
	input:
		fully_phased="candidates/{genome}_{sample}_fully_phased.txt",
		h1_fully_phased="tmp/{genome}_{sample}/fully_phased_locus_fasta_h1",
		h2_fully_phased="tmp/{genome}_{sample}/fully_phased_locus_fasta_h2",

		one_miss="candidates/{genome}_{sample}_one_phased_missing.txt",
		h1_one_miss="tmp/{genome}_{sample}/one_miss_locus_fasta_h1",
		h2_one_miss="tmp/{genome}_{sample}/one_miss_locus_fasta_h2"

	output:
		tmp_fully_phased=directory("tmp/{genome}_{sample}/fully_phased_exonerate"),
		tmp_one_miss=directory("tmp/{genome}_{sample}/one_miss_exonerate"),

		result_fully_phased=directory("output_cds/{genome}_{sample}/fully_phased"),
		result_one_miss=directory("output_cds/{genome}_{sample}/one_miss")
	params:
		proteome="proteome/no_stop_codons",
		sample="{sample}",
		h1=1,
		h2=2,
#		tidy_up="rm *tmp.bcf"
	conda:
		"envs/cds_extraction.yaml"
	threads: 2
	shell:
		"mkdir {output.tmp_fully_phased}; mkdir {output.result_fully_phased};"
		"cat {input.fully_phased} | xargs -I@ -n 1 -P 2 sh -c 'scripts/exonerate.sh @ {input.h1_fully_phased} {params.proteome} {output.tmp_fully_phased} {output.result_fully_phased} {params.sample} {params.h1}';"
		"cat {input.fully_phased} | xargs -I@ -n 1 -P 2 sh -c 'scripts/exonerate.sh @ {input.h2_fully_phased} {params.proteome} {output.tmp_fully_phased} {output.result_fully_phased} {params.sample} {params.h2}';"

		"mkdir {output.tmp_one_miss}; mkdir {output.result_one_miss};"
		"cat {input.one_miss} | xargs -I@ -n 1 -P 2 sh -c 'scripts/exonerate.sh @ {input.h1_one_miss} {params.proteome} {output.tmp_one_miss} {output.result_one_miss} {params.sample} {params.h1}';"
		"cat {input.one_miss} | xargs -I@ -n 1 -P 2 sh -c 'scripts/exonerate.sh @ {input.h2_one_miss} {params.proteome} {output.tmp_one_miss} {output.result_one_miss} {params.sample} {params.h2}'"
#		"{params.tidy_up}"


#		"mkdir {output.tmp}; mkdir {output.result};"
## NEXT: I'm gettin an error in the exonerate search (initial HSP -1???) doesn't happen with many, but worth looking into
## Then: generate second haplotype as well and try pipeline for more than one species
