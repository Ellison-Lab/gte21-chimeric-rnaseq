rule get_fusion_support:
    input:
        rules.star_align.output,
    output:
        csv = "results/breakpoints/{s}.breakpoints.csv"
    params:
        junctions = "results/star/{s}/Chimeric.out.junction"
    conda:
        "../envs/bioc-general.yaml"
    scripts:
        "../scripts/get_fusion_support.R"
