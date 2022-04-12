# compare DMP sets from merid et al 

# load merid tables
csv.fname1 <- "dmps1-nocomplications_cord-validation_merid-et-al-2020.csv"
csv.fname2 <- "dmps2-allbirths_cord-validation_merid-et-al-2020.csv"
mt1 <- read.csv(file.path("tables", csv.fname1))
mt2 <- read.csv(file.path("tables", csv.fname2))

# load cord blood ga dmps
st <- get(load("st_gest-age-dmp_all-dmp-info.rda"))

# compare overlaps


