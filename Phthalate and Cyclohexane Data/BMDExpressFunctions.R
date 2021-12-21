
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

