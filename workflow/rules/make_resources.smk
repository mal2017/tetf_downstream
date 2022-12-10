rule make_gene_symbol_lookup:
    output:
        tsv = "results/resources/gene_symbol_lookup.tsv.gz"
    script:
        "../scripts/resource_making/gene_symbol_lookup.R"

rule make_pirna_gene_list:
    input:
        lkup = rules.make_gene_symbol_lookup.output.tsv,
        handler = "resources/handler2013_supp2.xlsx",
        czech = "resources/czech2013_supp2.xlsx",
    output:
        tsv = "results/resources/pirna_pathway.tsv"
    script:
        "../scripts/resource_making/get_piRNA_genes.R"


rule annotate_insertion_penetrance:
    input:
        genome_fa = config.get("GENOME_FA"),
        ins = config.get("DGRP_TIDAL"),
    output:
        all_ins = "results/resources/dgrp_tidal_insertions.bb",
        penetrance = "results/resources/dgrp_tidal_insertions.unique.bb"
    script:
        "../scripts/resource_making/annotate_insertion_penetrance.R"

rule annotate_fixed_insertions:
    """
    Insertions present in the ref (per my repeatmasker run) that are also fixed across all TIDAL strains.
    Annotated with feature overlap, gc, etc.
    """
    input:
        lms = config.get("MERGED_MODELS"),
        remap = "resources/remap2022_nr_macs2_dm6_v1_0.bed.gz",
        het = "resources/het_domains_r5todm6.bed.gz",
        s2_chromstate = "resources/chromstate_s2_9state_r5todm6.bed.gz",
        insertions = rules.annotate_insertion_penetrance.output.penetrance,
        all_ins = rules.annotate_insertion_penetrance.output.all_ins,
        gtf = config.get("GTF"),
    output:
        rds = "results/resources/annotated_fixed_insertions.gr.rds",
        remap = "results/resources/remap.gr.rds",
    script:
        "../scripts/resource_making/annotate_fixed_insertions.R"

rule ref_preproc:
    """
    make a enome fasta with stripped names
    """
    input:
        genome_fa = config.get("GENOME_FA"),
    output:
        genome_fa = "results/resources/genome.fasta",
    script:
        "../scripts/resource_making/ref_preprocessing.R"

rule make_txdb:
    """
    make a reloadable txdb a
    """
    input:
        gtf = config.get("GTF"),
    output:
        txdb = "results/resources/txdb",
    script:
        "../scripts/resource_making/make_txdb.R"