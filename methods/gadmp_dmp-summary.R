# load data
st.fname <- "st_gest-age-dmp_all-dmp-info.rda"
dmpdf <- get(load(file.path(st.fname)))

# evaluate cord blood DMPs
table(dmpdf$is.cord.blood.dmp)
# FALSE  TRUE 
# 16926  1097 
table(dmpdf$is.cord.blood.dmp, dmpdf$is.haftorn.2021.dmp)
58/397 # 15%
table(dmpdf$is.cord.blood.dmp, dmpdf$is.merid.2020.dmp)
1005/17686 # 6%