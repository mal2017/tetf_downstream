rule mask_cons_tes:
    input:
        tes = config.get("TE_FA"),
    output:
        masked = "results/analysis/motifs/bbmask_consensus_tes/consensus_tes.masked.fasta",
    params:
        entropy = config.get("BBMASK_ENTROPY"),
        w = config.get("BBMASK_W"),
    singularity:
        "docker://quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0"
    shell:
        """
        bbmask.sh \
            in={input.tes} \
            out={output.masked} \
            w={params.w} \
            entropy={params.entropy}
        """

checkpoint split_cons_tes_per_tf:
    input:
        tfs = config.get("TFS"),
        mods = config.get("MERGED_MODELS"),
        tes = rules.mask_cons_tes.output.masked
    output:
        odir = directory("results/analysis/motifs/consensus_tes_per_tf/")
    script:
        "../scripts/motifs/split_consensus_tes_per_tf.R"

checkpoint split_cons_tes_per_tf_unmasked:
    input:
        tfs = config.get("TFS"),
        mods = config.get("MERGED_MODELS"),
        tes = config.get("TE_FA"),
    output:
        odir = directory("results/analysis/motifs/consensus_tes_per_tf_unmasked/")
    script:
        "../scripts/motifs/split_consensus_tes_per_tf.R"

rule shuffle_cons:
    input:
        dir = rules.split_cons_tes_per_tf.output.odir,
    output:
        p = "results/analysis/motifs/shuffle_cons/{tf}/coex.{shuf_rep}.fasta",
        n = "results/analysis/motifs/shuffle_cons/{tf}/other.{shuf_rep}.fasta",
    singularity:
        "docker://memesuite/memesuite:5.5.0"
    params:
        k = 9,
    shell:
        """
        fasta-shuffle-letters -dna -kmer {params.k} --seed {wildcards.shuf_rep}{wildcards.shuf_rep} '{input.dir}/{wildcards.tf}/coex.fasta' '{output.p}'
        fasta-shuffle-letters -dna -kmer {params.k} --seed {wildcards.shuf_rep}{wildcards.shuf_rep} '{input.dir}/{wildcards.tf}/other.fasta' '{output.n}'
        """

rule mask_cons_tes_shuf:
    input:
        coex = rules.shuffle_cons.output.p,
        other = rules.shuffle_cons.output.n,
    output:
        coex = "results/analysis/motifs/shuffle_cons/{tf}/coex.{shuf_rep}.masked.fasta",
        other = "results/analysis/motifs/shuffle_cons/{tf}/other.{shuf_rep}.masked.fasta",
    params:
        entropy = config.get("BBMASK_ENTROPY"),
        w = config.get("BBMASK_W"),
    singularity:
        "docker://quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0"
    shell:
        """
        bbmask.sh \
            in={input.coex} \
            out={output.coex} \
            w={params.w} \
            entropy={params.entropy}
        
        bbmask.sh \
            in={input.other} \
            out={output.other} \
            w={params.w} \
            entropy={params.entropy}
        """

rule xstreme_per_tf_shuffled:
    input:
        p = rules.mask_cons_tes_shuf.output.coex,
        #n = rules.mask_cons_tes_shuf.output.other,
        dir = rules.split_cons_tes_per_tf.output.odir,
    output:
        odir = directory("results/analysis/motifs/xstreme_per_tf_shuffled/{tf}/{shuf_rep}")
    threads:
        1
    singularity:
        "docker://memesuite/memesuite:5.5.0"
    shell:
        """
        xstreme --oc '{output.odir}' \
            -p '{input.p}' \
            -n '{input.dir}/{wildcards.tf}/other.fasta' \
            --meme-nmotifs 0 --meme-p {threads}
        """


rule xstreme_per_tf:
    input:
        dir = rules.split_cons_tes_per_tf.output.odir,
    output:
        odir = directory("results/analysis/motifs/xstreme_per_tf/{tf}/")
    threads:
        1
    singularity:
        "docker://memesuite/memesuite:5.5.0"
    shell:
        """
        xstreme --oc '{output.odir}' \
            -p '{input.dir}/{wildcards.tf}/coex.fasta' \
            -n '{input.dir}/{wildcards.tf}/other.fasta' \
            --meme-nmotifs 0 --meme-p {threads}
        """

# def aggregate_xstreme(wildcards):
#     checkpoint_output = checkpoints.split_cons_tes_per_tf.get(**wildcards).output.odir
#     wc_path = os.path.join(checkpoint_output, "{tf}/coex.fasta")
#     tfs = glob_wildcards(wc_path).tf
#     tfs = ["pan"]
#     return expand("results/analysis/motifs/xstreme_per_tf/{tf}/", tf=tfs)

# def aggregate_xstreme_shuf(wildcards):
#     checkpoint_output = checkpoints.split_cons_tes_per_tf.get(**wildcards).output.odir
#     wc_path = os.path.join(checkpoint_output, "{tf}/coex.fasta")
#     tfs = glob_wildcards(wc_path).tf
#     tfs = ["pan"]
#     return expand("results/analysis/motifs/xstreme_per_tf_shuffled/{tf}/{shuf_rep}", tf=tfs, shuf_rep=list(range(1, 21)))

TFSOI=["pan"]

rule combine_xstreme_motifs:
    input:
        memes = expand("results/analysis/motifs/xstreme_per_tf/{tf}/", tf=TFSOI),
    output:
        meme = "results/analysis/motifs/combined_xstreme.meme"
    script:
        "../scripts/motifs/combine_xstreme_motifs.R"

rule combine_xstreme_results:
    input:
        memes = expand("results/analysis/motifs/xstreme_per_tf/{tf}/", tf=TFSOI),
        shuf = expand("results/analysis/motifs/xstreme_per_tf_shuffled/{tf}/{shuf_rep}", tf=TFSOI, shuf_rep=list(range(1, 101)))
    output:
        tsv = "results/analysis/motifs/combined_xstreme_results.tsv"
    script:
        "../scripts/motifs/combine_xstreme_results.R"

        
# checkpoint get_remap_peak_seqs:
#     input:
#         bed = rules.annotate_fixed_insertions.output.remap,
#         rpm = config.get("REPEATMASKER_BED"),
#         fa = config.get("GENOME_FA")
#     output:
#         odir = directory("results/analysis/motifs/remap_peaks/")
#     script:
#         "../scripts/motifs/get_remap_peak_seqs.R"
    
# rule get_inverted_remap_peak_seqs:
#     input:
#         bed = rules.annotate_fixed_insertions.output.remap,
#         rpm = config.get("REPEATMASKER_BED"),
#         fa = config.get("GENOME_FA")
#     params:
#         tf = "{tf}"
#     output:
#         fa = "results/analysis/motifs/inverted_remap_peaks/{tf}.fasta"
#     script:
#         "../scripts/motifs/get_inverted_remap_peak_seqs.R"


# rule sea_remap_peaks:
#     input:
#         dir = rules.get_remap_peak_seqs.output.odir,
#         nfa = rules.get_inverted_remap_peak_seqs.output.fa,
#         #xstreme = rules.combine_xstreme_motifs.output.meme # to run on all motifs, comment this out and delete '/combined.meme' from sea command
#         xstreme = rules.xstreme_per_tf.output.odir,
#     output:
#         odir = directory("results/analysis/motifs/sea_remap_peaks/{tf}")
#     singularity:
#         "docker://memesuite/memesuite:5.5.0"
#     shell:
#         """
#         sea -p '{input.dir}/{wildcards.tf}.fasta' -n {input.nfa} -m '{input.xstreme}/combined.meme' -oc '{output.odir}'
#         """

# rule xstreme_remap_peaks:
#     """
#     --n '{input.nfa}' \
#     """
#     input:
#         dir = rules.get_remap_peak_seqs.output.odir,
#         #nfa = rules.get_inverted_remap_peak_seqs.output.fa,
#     output:
#         odir = directory("results/analysis/motifs/xstreme_remap_peaks/{tf}/")
#     threads:
#         8
#     singularity:
#         "docker://memesuite/memesuite:5.5.0"
#     shell:
#         """
#         xstreme --oc '{output.odir}' \
#             -p '{input.dir}/{wildcards.tf}.fasta' \
#             --meme-p {threads}
#         """

# rule sea_inverted_remap_peaks:
#     input:
#         fa =  rules.get_inverted_remap_peak_seqs.output.fa,
#         xstreme = rules.combine_xstreme_motifs.output.meme #rules.xstreme_per_tf.output.odir,
#     output:
#         odir = directory("results/analysis/motifs/sea_inverted_remap_peaks/{tf}")
#     singularity:
#         "docker://memesuite/memesuite:5.5.0"
#     shell:
#         """
#         sea -p {input.fa} -m '{input.xstreme}' -oc '{output.odir}'
#         """

# def aggregate_sea(wildcards):
#     lms_checkpoint_output = checkpoints.split_cons_tes_per_tf.get(**wildcards).output.odir
#     remap_checkpoint_output = checkpoints.get_remap_peak_seqs.get(**wildcards).output.odir
#     #print(checkpoint_output)
#     wc_path = os.path.join(lms_checkpoint_output, "{tf}/coex.fasta")
#     tfs = glob_wildcards(wc_path).tf
#     filter_wc_path = os.path.join(remap_checkpoint_output, "{tf}.fasta")
#     filters = glob_wildcards(filter_wc_path).tf
#     filters = ["pan"]
#     return expand("results/analysis/motifs/sea_remap_peaks/{tf}", tf=[x for x in tfs if x in filters])


# rule collect_remap_peak_sea:
#     input:
#         seas = expand("results/analysis/motifs/sea_remap_peaks/{tf}", tf=TFSOI)
#     output:
#         tsv = "results/analysis/motifs/remap_peak_sea.tsv.gz"
#     script:
#         "../scripts/motifs/collect_remap_peak_sea.R"

# rule plot_remap_peak_sea:
#     input:
#         tsv = rules.collect_remap_peak_sea.output.tsv
#     output:
#         nhits_png = "results/plots/remap_peak_sea.nhits.png",
#         obs_exp_png = "results/plots/remap_peak_sea.obs_exp.png",
#         nhits_rds = "results/plots/remap_peak_sea.nhits.rds",
#         obs_exp_rds = "results/plots/remap_peak_sea.obs_exp.rds",
#     script:
#         "../scripts/motifs/plot_remap_peak_sea.R"

rule compare_hmg_motifs_from_archbold14:
    input:
        meme = rules.combine_xstreme_motifs.output.meme,
    output:
        motif_comparison = "results/analysis/motifs/archbold14_motif_comparison.rds",
        motif_similarity = "results/analysis/motifs/archbold14_motif_similarity.rds",
        archbold_motifs = "results/analysis/motifs/archbold14_motifs.tsv",
    script:
        "../scripts/motifs/compare_hmg_motifs_from_archbold14.R"

rule plot_combined_motif_analysis_table:
    input:
        memes = rules.combine_xstreme_motifs.output.meme,
        denovo = rules.combine_xstreme_results.output.tsv,
        #remap_enr = rules.collect_remap_peak_sea.output.tsv,
        archbold_compr = rules.compare_hmg_motifs_from_archbold14.output.motif_comparison,
        #shuf = expand("results/analysis/motifs/xstreme_per_tf_shuffled/{tf}/{shuf_rep}", tf=TFSOI, shuf_rep=list(range(1, 501)))
    output:
        png = "results/plots/combined_motif_analysis_table.png",
        rda = "results/plots/combined_motif_analysis_table.rda",
    script:
        "../scripts/motifs/plot_combined_motif_analysis_table.R"