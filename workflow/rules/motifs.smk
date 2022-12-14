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


rule xstreme_per_tf:
    input:
        dir = rules.split_cons_tes_per_tf.output.odir,
    output:
        odir = directory("results/analysis/motifs/xstreme_per_tf/{tf}/")
    threads:
        2
    singularity:
        "docker://memesuite/memesuite:5.5.0"
    shell:
        """
        xstreme --oc '{output.odir}' \
            -p '{input.dir}/{wildcards.tf}/coex.fasta' \
            --n '{input.dir}/{wildcards.tf}/other.fasta' \
            --meme-p {threads}
        """

def aggregate_xstreme(wildcards):
    checkpoint_output = checkpoints.split_cons_tes_per_tf.get(**wildcards).output.odir
    wc_path = os.path.join(checkpoint_output, "{tf}/coex.fasta")
    tfs = glob_wildcards(wc_path).tf
    tfs = ["pan","NfI","CG16779"]
    return expand("results/analysis/motifs/xstreme_per_tf/{tf}/", tf=tfs)


rule combine_xstreme_motifs:
    input:
        memes = aggregate_xstreme
    output:
        meme = "results/analysis/motifs/combined_xstreme.meme"
    script:
        "../scripts/motifs/combine_xstreme_motifs.R"

        
checkpoint get_remap_peak_seqs:
    input:
        bed = rules.annotate_fixed_insertions.output.remap,
        rpm = config.get("REPEATMASKER_BED"),
        fa = config.get("GENOME_FA")
    output:
        odir = directory("results/analysis/motifs/remap_peaks/")
    script:
        "../scripts/motifs/get_remap_peak_seqs.R"
    

rule sea_remap_peaks:
    input:
        dir = rules.get_remap_peak_seqs.output.odir,
        xstreme = rules.combine_xstreme_motifs.output.meme #rules.xstreme_per_tf.output.odir,
    output:
        odir = directory("results/analysis/motifs/sea_remap_peaks/{tf}")
    singularity:
        "docker://memesuite/memesuite:5.5.0"
    shell:
        """
        sea -p '{input.dir}/{wildcards.tf}.fasta' -m '{input.xstreme}' -oc '{output.odir}'
        """

def aggregate_sea(wildcards):
    lms_checkpoint_output = checkpoints.split_cons_tes_per_tf.get(**wildcards).output.odir
    remap_checkpoint_output = checkpoints.get_remap_peak_seqs.get(**wildcards).output.odir
    #print(checkpoint_output)
    wc_path = os.path.join(lms_checkpoint_output, "{tf}/coex.fasta")
    tfs = glob_wildcards(wc_path).tf
    filter_wc_path = os.path.join(remap_checkpoint_output, "{tf}.fasta")
    filters = glob_wildcards(filter_wc_path).tf
    filters = ["pan","NfI","CG16779"]
    return expand("results/analysis/motifs/sea_remap_peaks/{tf}", tf=[x for x in tfs if x in filters])


rule collect_remap_peak_sea:
    input:
        aggregate_sea
    output:
        tsv = "results/analysis/motifs/remap_peak_sea.tsv.gz"
    params:
        seas = lambda wc: [x + "/sea.tsv" for x  in aggregate_sea(wc)]
    script:
        "../scripts/motifs/collect_remap_peak_sea.R"


rule plot_remap_peak_sea:
    input:
        tsv = rules.collect_remap_peak_sea.output.tsv
    output:
        nhits_png = "results/figs/remap_peak_sea.nhits.png",
        obs_exp_png = "results/figs/remap_peak_sea.obs_exp.png",
        nhits_rds = "results/figs/remap_peak_sea.nhits.rds",
        obs_exp_rds = "results/figs/remap_peak_sea.obs_exp.rds",
    script:
        "../scripts/motifs/plot_remap_peak_sea.R"


rule compare_hmg_motifs_from_archbold14:
    input:
        meme = rules.combine_xstreme_motifs.output.meme,
    output:
        motif_comparison = "results/analysis/motifs/archbold14_motif_comparison.rds",
        motif_similarity = "results/analysis/motifs/archbold14_motif_similarity.rds",
        archbold_motifs = "results/analysis/motifs/archbold14_motifs.tsv",
    script:
        "../scripts/motifs/compare_hmg_motifs_from_archbold14.R"