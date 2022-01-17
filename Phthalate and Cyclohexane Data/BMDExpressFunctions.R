#### BMD Functions ####
# Clean up columns during import
cleanupcolumns <- function(x){
  x %>%
    select(Probe.Id,
           BMD,
           BMDL,
           BMDU,
           fitPValue)
}

# BMD Filtering
# BMD/BMDL > 20, BMDU/BMDL > 40, BMDU/BMD > 20, BMD > lowdose, BMD < highdose
BMDfiltering <- function(x,
                         BMD.div.BMDL = 20,
                         BMDU.div.BMDL = 40,
                         BMDU.div.BMD = 20,
                         lowdose,
                         highdose,
                         fitP = 0.1){
  x %>%
    filter(
      BMD/BMDL < BMD.div.BMDL,
      BMDU/BMDL < BMDU.div.BMDL,
      BMDU/BMD < BMDU.div.BMD,
      BMD > lowdose/10,
      BMD <= highdose,
      fitPValue > fitP
      )
}

####REACTOME Functions####
#Clean up columns during import
cleanupcolumns_reactome <- function(x){
  x %>%
    select(GO.Pathway.Gene.Set.Gene.ID,
           GO.Pathway.Gene.Set.Gene.Name,
           All.Genes..Expression.Data.,
           Genes.That.Passed.All.Filters,
           Fisher.s.Exact.Two.Tail,
           Percentage,
           BMD.Mean,
           BMD.Median
           )
}

# Reactome Filtering
# Two Tail Test < 0.05, genes that passed all filters > 3
reactome_filtering <- function(x,
                               p = 0.05,
                               min_gene = 3){
  x %>%
    filter(
      p >= Fisher.s.Exact.Two.Tail,
      min_gene <= Genes.That.Passed.All.Filters
    )
}

####GO Term Functions####
#Clean up columns
cleanupcolumns_goterm <- function(x){
  x %>%
    select(GO.Pathway.Gene.Set.Gene.ID,
           GO.Level,
           GO.Pathway.Gene.Set.Gene.Name,
           All.Genes..Expression.Data.,
           Genes.That.Passed.All.Filters,
           Fisher.s.Exact.Two.Tail,
           Percentage,
           BMD.Mean,
           BMD.Median
    )
}

# GO-term Filtering
# Two Tail Test < 0.05, genes that passed all filters > 3, GO Term level > 5
goterm_filtering <- function(x,
                               p = 0.05,
                               min_gene = 3,
                               min_level = 5){
  x %>%
    filter(
      p >= Fisher.s.Exact.Two.Tail,
      min_gene <= Genes.That.Passed.All.Filters,
      min_level <= GO.Level
    )
}
