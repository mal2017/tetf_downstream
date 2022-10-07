library(tidyverse)
library(readxl)

#lkup_path <- "results/resources/gene_symbol_lookup.tsv.gz"
lkup_path <- snakemake@input[["lkup"]]

#handler_path <- "data/handler2013_supp2.xlsx"
handler_path <- snakemake@input[["handler"]]

#czech_path <- "data/czech2013_supp2.xlsx"
czech_path <- snakemake@input[["czech"]]



lkup <- read_tsv(lkup_path)

handler <- readxl::read_xlsx(handler_path, skip = 1) %>% 
  filter(staining...6 > 1) %>% # this is a low bar
  dplyr::select(gene_symbol=Symbol,annotation_ID=CG) %>%
  distinct()

# See Note1 below
handler <- handler %>% mutate(annotation_ID = ifelse(gene_symbol == "trpm","CG44240",annotation_ID))


czech <- readxl::read_xlsx(czech_path) %>%
  mutate(across(c("HeTA","TAHRE","blood","burdock"),as.numeric)) %>%
  filter(if_all(c("HeTA","TAHRE","blood","burdock"),~{. < -2})) %>%
  dplyr::select(gene_symbol=Synonym,annotation_ID=Gene) %>%
  distinct()

czech <- czech %>% mutate(annotation_ID= ifelse(gene_symbol == "vas","CG46283",annotation_ID))

# handler2013 is old enough that some of their genes identified by
# annotation_ID (e.g. CG2183~Gasz) now have names.
handler_cg_matches <- handler %>% 
  left_join(lkup, by=c("annotation_ID"),suffix=c(".handler13","")) %>% 
  distinct() %>%
  drop_na()


czech_cg_matches <- czech %>% 
  left_join(lkup, by=c("annotation_ID"),suffix=c(".czech13","")) %>% 
  distinct() %>%
  drop_na()


cg_matches <- full_join(mutate(handler_cg_matches,in.Handler13 = T),
                        mutate(czech_cg_matches,in.Czech13=T),
                        by=c("annotation_ID", "gene_ID", "gene_symbol", "gene_type"))

# to sanity check, group by different cols and tally to make sure no doubles.
cg_matches %>%
  group_by(annotation_ID) %>%
  tally(sort = T) %>%
  pull(n) %>%
  max() %>%
  {.==1} %>%
  stopifnot()

# NOTE 1:
# that said, we lost 1 gene that wasn't findable, in each dataset, soi went back and fixed this above
stopifnot(length(handler$annotation_ID[!handler$annotation_ID %in% cg_matches$annotation_ID])==0)
stopifnot(length(czech$annotation_ID[!czech$annotation_ID %in% cg_matches$annotation_ID])==0)

# NOTE 2: 
# this leaves 30 genes with different gene symbols in handler vs flybase
# when matching by CG id.
# i resolve these manually by confirming
# that either the handler gene symbol makes sense or the flybase entry
# lists handler's cg ID as alternate
confirmed_ok <- c("CG2183","CG6768","CG9821","CG3893","CG9754","CG14749","CG3689",
                  "CG8949","CG8211","CG7504","CG11678","CG6020","CG10320","CG9140",
                  "CG5491","CG5222","CG13240","CG14213")

cg_matches %>%
  filter(gene_symbol.handler13 != gene_symbol | gene_symbol.czech13!=gene_symbol.handler13) %>%
  filter(str_remove_all(str_to_lower(gene_symbol),"[:punct:]") != str_remove_all(str_to_lower(gene_symbol.handler13),"[:punct:]")) %>%
  filter(!annotation_ID%in% confirmed_ok) %>%
  nrow() %>%
  {. == 0} %>%
  stopifnot()


cg_matches %>%
  dplyr::select(gene_ID,gene_symbol,in.Handler13,in.Czech13, gene_symbol.handler13, gene_symbol.czech13) %>%
  #filter(in.Czech13 & in.Handler13) # not very many, but confirmed this by looking at the spreadsheets
  write_tsv(snakemake@output[["tsv"]])