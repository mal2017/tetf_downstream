
# function gets mean in groups dist
# takes a kmer dist result previously generated (see above)
# takes a ig group of TEs and finds the mean distance within that group
# and excluding self-self comparisons (the diagonal)
get_dists <- function(x, d=kmer.dist) {
  as.matrix(d)[x,x] %>%
    .[upper.tri(.,diag = F)] %>%
    mean()
}


# updated to return bootstraps used (for sanity checks)
# uses the above function to asnwer the question:
# if we randomly selected n TEs, from TEs not in the group of interest,
# what is the expected in group dist? THis is measured over a number replicates
# (currently hardcoded to 10) to make this more robust than a single sample
get_boot_dists2 <- function(n,exclude) {
  #pool <- te.names[!te.names %in% exclude]
  #print(paste(n,length(pool)))
  pool <- te.names
  replicate(10,sample(pool,size = min(n,0.9*length(pool)), replace = F),simplify = F) %>%
    tibble(matched_boots = .) %>%
    mutate(matched.dist = map_dbl(matched_boots, get_dists))
}

get_boot_dists3 <- function(nTEs,n=1000) {
  # te.names is in the global env
  replicate(n,sample(te.names,size = nTEs, replace = F),simplify = F) %>%
    map_dbl(get_dists)
}
