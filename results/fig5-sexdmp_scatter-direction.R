#!/usr/bin/env R

# Author: Sean Maden
#
# Scatter plot of DNAm means direction between sexes at DMPs. These 
# plots show the mean DNAm differences (male - female) for 
# Inoshita et al 2015 and PBMC and whole blood, highlighting where
# direction is concordant or discrepant between validations.
# 

library(ggplot2)
library(patchwork) # for inset plots

#----------
# load data
#----------
# load dmp sets info
dmp.fpath <- "st_sex-dmp_all-dmp-info.rda"
new.dmp <- get(load(file.path(dmp.fpath)))

# get plot palette
pal <- c("Other/NOS" = "#3B3B3BFF", 
         "Whole blood" = "#868686FF", 
         "Cord blood" = "#EFC000FF", 
         "PBMC" = "#CD534CFF")

#------------------
# plot whole blood
#------------------
# get plot vars
col.wb <- "#868686FF"
pt.alpha <- 0.5
filt.wb <- new.dmp$is.whole.blood.dmp == T & new.dmp$is.inoshita.2015.dmp == T
df.wb <- new.dmp[filt.wb,]

# main scatter plot
pt.wb <- ggplot() + theme_bw() +
  geom_point(data = df.wb, size = 2,
             aes(x = dnam.direction.inoshita.2015.dmp,
                 y = dnam.direction.whole.blood.dmp), 
             alpha = pt.alpha, color = col.wb) +
  theme(plot.margin = margin(0.005, 0.005, 0.08, 0.08, unit = "cm")) +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2),
                     limits = c(-0.25, 0.25)) +
  scale_y_continuous(breaks = c(-0.05, 0, 0.05), 
                     limits = c(-0.07, 0.07)) +
  xlab("Inoshita et al 2015\n(mean diff. M - F)") +
  ylab("Whole blood\n(mean diff. M - F)") +
  annotate("rect", xmin = 0, xmax = 0.25, 
           ymin = 0, ymax = 0.07, alpha = 0.2, 
           fill = "goldenrod") +
  annotate("rect", xmin = -0.25, xmax = 0, 
           ymin = -0.07, ymax = 0, alpha = 0.2, 
           fill = "goldenrod") +
  annotate("rect", xmin = -0.25, xmax = 0, 
           ymin = 0, ymax = 0.07, alpha = 0.2, 
           fill = "blue") +
  annotate("rect", xmin = 0, xmax = 0.25, 
           ymin = -0.07, ymax = 0, alpha = 0.2, 
           fill = "blue")
# make inset count matrix -- whole blood
filt.wb.pos <- df.wb$dnam.direction.whole.blood.dmp > 0
num.pos.wb <- nrow(df.wb[filt.wb.pos,])
num.neg.wb <- nrow(df.wb[!filt.wb.pos,])
dfi <- data.frame(x = c(1,1,2,2), y = c(1,2,1,2), 
                  lab = c(num.neg.wb, 0, 0, num.pos.wb),
                  value = c("a", "a","b","b"))
plot.inset <- ggplot(dfi, aes(x = x, y = y, label = lab)) + 
  geom_tile(alpha = 0.1) + geom_text(aes(label = lab)) +
annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 0.5, ymax = 1.5, 
         alpha = 0.2, fill = "goldenrod") +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = 1.5, ymax = 2.5, 
           alpha = 0.2, fill = "goldenrod") +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 1.5, ymax = 2.5, 
           alpha = 0.2, fill = "blue") +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = 0.5, ymax = 1.5, 
           alpha = 0.2, fill = "blue") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0, "cm"),
        line = element_blank(),
        axis.ticks = element_blank())

# save composite plot
plot.composite <- pt.wb + inset_element(plot.inset, 
                                        left = 0.57, 
                                        bottom = 0.08, 
                                        right = 0.92, 
                                        top = 0.46)

# save new pdf
pdf.fname <- "fig4e_sex-dmp_wb-dnam-diff-with-inset.pdf"
pdf(pdf.fname, 2.25, 2.05)
print(plot.composite)
dev.off()

#----------
# plot pbmc
#----------
# get plot vars
col.pbmc <- "#CD534CFF"
pt.alpha <- 0.5
filt.pbmc <- new.dmp$is.pbmc.dmp == T &
  new.dmp$is.inoshita.2015.dmp == T
df.pbmc <- new.dmp[filt.pbmc,]

# main scatter plot
pt.pbmc <- ggplot() + theme_bw() +
  geom_point(data = df.pbmc, size = 2,
             aes(x = dnam.direction.inoshita.2015.dmp,
                 y = dnam.direction.pbmc.dmp), 
             alpha = pt.alpha, color = col.pbmc) +
  theme(plot.margin = margin(0.005, 0.005, 0.08, 0.08, unit = "cm")) +
  scale_y_continuous(breaks = c(-0.1, 0, 0.1),
                     limits = c(-0.15, 0.15)) +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2),
                     limits = c(-0.25, 0.25)) +
  xlab("Inoshita et al 2015\n(mean diff. M - F)") +
  ylab("PBMC\n(mean diff. M - F)") +
  annotate("rect", xmin = 0, xmax = 0.25, 
           ymin = 0, ymax = 0.15, alpha = 0.2, 
           fill = "goldenrod") +
  annotate("rect", xmin = -0.25, xmax = 0, 
           ymin = -0.15, ymax = 0, alpha = 0.2, 
           fill = "goldenrod") +
  annotate("rect", xmin = -0.25, xmax = 0, 
           ymin = 0, ymax = 0.15, alpha = 0.2, 
           fill = "blue") +
  annotate("rect", xmin = 0, xmax = 0.25, 
           ymin = -0.15, ymax = 0, alpha = 0.2, 
           fill = "blue")

# make inset count matrix -- whole blood
filt.pbmc.pos <- df.pbmc$dnam.direction.pbmc.dmp > 0
num.pos.pbmc <- nrow(df.pbmc[filt.pbmc.pos,])
num.neg.pbmc <- nrow(df.pbmc[!filt.pbmc.pos,])
dfi <- data.frame(x = c(1,1,2,2), y = c(1,2,1,2), 
                  lab = c(num.neg.pbmc, 0, 0, num.pos.pbmc),
                  value = c("a", "a","b","b"))
plot.inset <- ggplot(dfi, aes(x = x, y = y, label = lab)) + #theme_bw() +
  geom_tile(alpha = 0.1) + geom_text(aes(label = lab)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 0.5, ymax = 1.5, 
           alpha = 0.2, fill = "goldenrod") +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = 1.5, ymax = 2.5, 
           alpha = 0.2, fill = "goldenrod") +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 1.5, ymax = 2.5, 
           alpha = 0.2, fill = "blue") +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = 0.5, ymax = 1.5, 
           alpha = 0.2, fill = "blue") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0, "cm"),
        line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# save composite plot
plot.composite <- pt.pbmc + inset_element(plot.inset, 
                                          left = 0.57, 
                                          bottom = 0.08, 
                                          right = 0.92, 
                                          top = 0.46)

# save new pdf
pdf.fname <- "fig4f_sex-dmp_pbmc-dnam-diff-with-inset.pdf"
pdf(pdf.fname, 2.25, 2.05)
print(plot.composite)
dev.off()