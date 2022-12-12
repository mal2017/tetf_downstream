checkpoint split_cons_tes_per_tf:
    input:
        tfs = config.get("TFS"),
        mods = config.get("MERGED_MODELS"),
        tes = config.get("TE_FA"),
    output:
        odir = directory("results/analysis/motifs/consensus_tes_per_tf/")
    script:
        "../scripts/motifs/split_consensus_tes_per_tf.R"


rule xstreme_per_tf:
    input:
        dir = rules.split_cons_tes_per_tf.output.odir,
    output:
        odir = directory("results/analysis/motifs/xstreme_per_tf/{tf}")
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
        
checkpoint get_remap_peak_seqs:
    input:
        bed = rules.annotate_fixed_insertions.output.remap,
        fa = config.get("GENOME_FA")
    output:
        odir = directory("results/analysis/motifs/remap_peaks/")
    script:
        "../scripts/motifs/get_remap_peak_seqs.R"
    

rule sea_remap_peaks:
    input:
        dir = rules.get_remap_peak_seqs.output.odir,
        xstreme = rules.xstreme_per_tf.output.odir,
    output:
        odir = directory("results/analysis/motifs/sea_remap_peaks/{tf}")
    threads:
        2
    singularity:
        "docker://memesuite/memesuite:5.5.0"
    shell:
        """
        sea -p '{input.dir}/{wildcards.tf}.fasta' -m '{input.xstreme}/combined.meme' -oc '{output.odir}'
        """

def aggregate_sea(wildcards):
    lms_checkpoint_output = checkpoints.split_cons_tes_per_tf.get(**wildcards).output.odir
    remap_checkpoint_output = checkpoints.get_remap_peak_seqs.get(**wildcards).output.odir
    #print(checkpoint_output)
    wc_path = os.path.join(lms_checkpoint_output, "{tf}/coex.fasta")
    tfs = glob_wildcards(wc_path).tf
    filter_wc_path = os.path.join(remap_checkpoint_output, "{tf}.fasta")
    filters = glob_wildcards(filter_wc_path).tf
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