rule within_tf_regioner:
    input:
        remap = rules.annotate_fixed_insertions.output.remap,
        unique_ins = rules.annotate_insertion_penetrance.output.penetrance,
        genome_fa = config.get("GENOME_FA"),
        mods = config.get("MERGED_MODELS"),
    threads:
        8
    output:
        tsv = "results/analysis/direct_binding/within_tf_regioner.tsv.gz",
    script:
        "../scripts/direct_binding/regioneR.R"

rule per_pair_nullranges:
    input:
        lms = config.get("MERGED_MODELS"),
        anno_ins = rules.annotate_fixed_insertions.output.rds,
    output:
        rds = "results/analysis/direct_binding/per_pair_nullranges.rds",
    threads:
        8
    script:
        "../scripts/direct_binding/per_pair_nullranges.R"

rule plot_per_pair_nullranges:
    input:
        rds = rules.per_pair_nullranges.output.rds
    output:
        rds = 'results/plots/per_pair_nullranges.rds',
        png = 'results/plots/per_pair_nullranges.png'
    script:
        "../scripts/direct_binding/plot_per_pair_nullranges.R"

rule per_pair_bootranges:
    input:
        lms = config.get("MERGED_MODELS"),
        anno_ins = rules.annotate_fixed_insertions.output.rds,
        txdb = rules.make_txdb.output.txdb,
        genome_fa = rules.ref_preproc.output.genome_fa,
        remap = rules.annotate_fixed_insertions.output.remap,
    output:
        rds = "results/analysis/direct_binding/per_pair_bootranges.rds",
        seg_rds = "results/analysis/direct_binding/per_pair_bootranges.seg.rds",
    params:
        R = config.get("BOOTRANGES_REPS"),
        blockLength = config.get("BOOTRANGES_BLOCKLENGTH"),
        L_s = config.get("BOOTRANGES_SEGLENGTH"), 
        nseg = config.get("BOOTRANGES_NSEG"),
    threads:
        8
    script:
        "../scripts/direct_binding/per_pair_bootranges.R"


