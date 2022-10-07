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