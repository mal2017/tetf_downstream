rule get_mash_dist:
    input:
        te_fasta = config.get("TE_FA")
    output:
        mash_dist = "results/analysis/direct_binding/te_mash.txt",
    conda:
        "../envs/mash.yaml"
    params:
        args = config.get("MASH_ARGS"),
    shell:
        """
        mash dist {params.args} -i {input.te_fasta} {input.te_fasta} > {output.mash_dist}
        """


rule each_gene_te_kmer_dist:
    input:
        te_fasta = config.get("TE_FA"),
        mash_dist = rules.get_mash_dist.output.mash_dist,
        lms = config.get("MERGED_MODELS"),
        tfs = config.get("TFS"),
    output:
        kmer_dist = "results/analysis/direct_binding/te_kmer_dist_mat.rds",
        each_gene = "results/analysis/direct_binding/each_gene_te_kmer_dist.tbl.rds",
        null_model = "results/analysis/direct_binding/null_model.rds",
    params:
        nperm = config.get("EXPECTED_MASH_DIST_NPERM")
    threads:
        8
    script:
        "../scripts/te_sequence_similarity/te_kmer_distance_each_gene_2.R"

rule plot_each_gene_te_kmer_dist:
    input:
        boots = rules.each_gene_te_kmer_dist.output.each_gene,
        kmer_dist = rules.each_gene_te_kmer_dist.output.kmer_dist,
        null_model = rules.each_gene_te_kmer_dist.output.null_model,
        tfs = config.get("TFS"),
        cofacs = config.get("COFACS"),
    output:
        rds = "results/figs/each_gene_te_kmer_dist.list.rds",
        pdf =  "results/figs/each_gene_te_kmer_dist.pdf",
    script:
        "../scripts/te_sequence_similarity/plot_intra_gene_te_kmer_distance.R"