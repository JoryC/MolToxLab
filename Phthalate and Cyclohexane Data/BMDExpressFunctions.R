
cleanupcolumns <- function(x){
  x %>%
    select(Probe.Id,
           BMD,
           BMDL,
           BMDU,
           fitPValue)
}

# BMD/BMDL > 20, BMDU/BMDL > 40, BMDU/BMD > 20, BMD > lowdose, BMD < highdose
BMDfiltering <- function(x,
                         BMD.div.BMDL,
                         BMDU.div.BMDL,
                         BMDU.div.BMD,
                         lowdose,
                         highdose){
  x %>%
    filter(
      BMD/BMDL > BMD.div.BMDL,
      BMDU/BMDL > BMDU.div.BMDL,
      BMDU/BMD > BMDU.div.BMD,
      BMD > lowdose/10,
      BMD < highdose
      )
}

