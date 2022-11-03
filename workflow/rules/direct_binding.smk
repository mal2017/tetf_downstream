rule get_mash_dist:
    input:
        te_fasta = "data/Tidalbase_transposon_sequence.fasta",
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
        te_fasta = "data/Tidalbase_transposon_sequence.fasta",
        mash_dist = rules.get_mash_dist.output.mash_dist,
        lms = rules.filter_models.output.filtered_tsv,
    output:
        kmer_dist = "results/analysis/direct_binding/te_kmer_dist_mat.rds",
        each_gene = "results/analysis/direct_binding/each_gene_te_kmer_dist.tbl.rds",
    threads:
        8
    script:
        "../scripts/direct_binding/te_kmer_distance_each_gene.R"

rule plot_each_gene_te_kmer_dist:
    input:
        boots = rules.each_gene_te_kmer_dist.output.each_gene,
        kmer_dist = rules.each_gene_te_kmer_dist.output.kmer_dist,
        tfs = "data/Drosophila_melanogaster_TF.txt",
        cofacs = "data/Drosophila_melanogaster_TF_cofactors.txt",
    output:
        rds = "results/plots/plot_each_gene_te_kmer_dist.rds"
    script:
        "../scripts/direct_binding/plot_intra_gene_te_kmer_distance.R"


rule coregulated_te_communities:
    input:
        lkup = rules.make_gene_symbol_lookup.output.tsv,
        filtered_mods = rules.filter_models.output.extreme_tsv,
        merged_mods = rules.filter_models.output.merged_tsv,
        tfs = "data/Drosophila_melanogaster_TF.txt",
    output:
        igraph = "results/analysis/direct_binding/coregulated_te_communities.igraph_list.rds",
        comms = "results/analysis/direct_binding/coregulated_te_communities.communities_tbl.rds",
    script:
        "../scripts/direct_binding/coregulated_te_communities.R"


rule plot_coregulated_te_communities:
    input:
        igraph = rules.coregulated_te_communities.output.igraph,
        comms = rules.coregulated_te_communities.output.comms,
    output:
        rds = "results/analysis/direct_binding/coregulated_te_communities.tbl_with_graph_layouts.rds"
    script:
        "../scripts/direct_binding/plot_coregulated_te_communities.R"

rule get_intra_community_kmer_distance:
    input:
        dist = rules.each_gene_te_kmer_dist.output.kmer_dist,
        comms = rules.coregulated_te_communities.output.comms,
    output:
        by_community = "results/analysis/direct_binding/by_community_te_kmer_dist.tbl.rds",
    script:
        "../scripts/direct_binding/get_intra_community_kmer_distance.R"

rule sequence_similarity_communities:
    input:
        filtered_mods = rules.filter_models.output.filtered_tsv,
        kmer_dist = rules.get_mash_dist.output.mash_dist,
    output:
        tg = "results/analysis/direct_binding/sequence_similarity_te_communities.tbl_graph.rds",
    script:
        "../scripts/direct_binding/sequence_similarity_communities.R"

rule within_tf_nullranges:
    input:
        lms = rules.filter_models.output.extreme_tsv,
        remap = "data/remap2022_nr_macs2_dm6_v1_0.bed.gz",
        het = "data/het_domains_r5todm6.bed.gz",
        s2_chromstate = "data/chromstate_s2_9state_r5todm6.bed.gz",
        bg3_chromstate = "data/chromstate_bg3_9state_r5todm6.bed.gz",
        insertions = "data/insertions_by_strain.tsv.gz",
    output:
        rds = "results/analysis/direct_binding/within_tf_nullranges.rds",
    script:
        "../scripts/direct_binding/within-tf-grp-nullranges.R"

rule plot_within_tf_nullranges:
    input:
        each_gene = rules.each_gene_te_kmer_dist.output.each_gene,
        rds = rules.within_tf_nullranges.output.rds,
    output:
        rds = 'results/plots/plot_within_tf_grp_nullranges.rds'
    script:
        "../scripts/direct_binding/plot-within-tf-grp-nullranges.R"


rule per_pair_nullranges:
    input:
        lms = rules.filter_models.output.extreme_tsv,
        anno_ins = rules.annotate_fixed_insertions.output.rds
    output:
        rds = 'results/analysis/direct_binding/per_pair_nullranges.rds'
    script:
        "../scripts/direct_binding/per_pair_nullranges.R"