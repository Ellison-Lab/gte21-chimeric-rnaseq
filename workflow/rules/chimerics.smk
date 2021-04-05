rule get_fusion_support:
    input:
        rules.star_align.output,
    output:
        csv = "results/breakpoints/{s}.breakpoints.csv"
    params:
        junctions = "results/star/{s}/Chimeric.out.junction"
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/get_fusion_support.R"

rule get_non_chimeric_depth:
    input:
        bam = rules.samtools_sort.output.bam,
        csv = rules.get_fusion_support.output.csv
    output:
        csv = "results/breakpoints/{s}.breakpoint-depths.csv"
    conda:
        "../envs/bioc-general.yaml"
    script:
        "../scripts/get_non_chimeric_depth.R"

rule get_total_reads_per_contig:
    input:
        bam = rules.samtools_sort.output.bam,
    output:
        csv = "results/breakpoints/{s}.total-contig-reads.csv"
    conda:
        "../envs/bioc-general.yaml"
    script:
        "../scripts/get_total_reads_per_contig.R"
