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