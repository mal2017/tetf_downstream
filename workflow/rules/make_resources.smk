rule make_gene_symbol_lookup:
    output:
        tsv = "results/resources/gene_symbol_lookup.tsv.gz"
    script:
        "../scripts/resource_making/gene_symbol_lookup.R"

rule make_pirna_gene_list:
    input:
        lkup = rules.make_gene_symbol_lookup.output.tsv,
        handler = "data/handler2013_supp2.xlsx",
        czech = "data/czech2013_supp2.xlsx",
    output:
        tsv = "results/resources/pirna_pathway.tsv"
    script:
        "../scripts/resource_making/get_piRNA_genes.R"


rule filter_models:
    """
    merged: all models -> mean of replicates
    filtered: merged -> all reps significant and no sign switches
    extreme: filtered -> filter at some extreme quantile of abs(coefs) (see script for exact choice)
    """
    input:
        expand("data/linear_models/{m}/{r}/lm.tidy.corrected.tsv.gz",m=["male_model_01","female_model_01"],r=[0,1,2])
    output:
        merged_tsv =  "results/resources/merged_models.tsv.gz",
        filtered_tsv =  "results/resources/filtered_models.tsv.gz",
        extreme_tsv =  "results/resources/extreme_models.tsv.gz"
    script:
        "../scripts/coexpression/consensus_lms.R"


rule annotate_fixed_insertions:
    """
    Insertions present in the ref (per my repeatmasker run) that are also fixed across all TIDAL strains.
    Annotated with feature overlap, gc, etc.
    """
    input:
        lms = rules.filter_models.output.filtered_tsv,
        remap = "data/remap2022_nr_macs2_dm6_v1_0.bed.gz",
        het = "data/het_domains_r5todm6.bed.gz",
        s2_chromstate = "data/chromstate_s2_9state_r5todm6.bed.gz",
        insertions = "data/insertions_by_strain.tsv.gz",
    output:
        rds = "results/resources/annotated_fixed_insertions.gr.rds",
    script:
        "../scripts/resource_making/annotate_fixed_insertions.R"