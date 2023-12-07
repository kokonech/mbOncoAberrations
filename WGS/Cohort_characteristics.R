##############################################################################################################################################
## source the settings
source("./RScripts/Settings.R")
library(wesanderson)
##############################################################################################################################################

sample.information$Gender <- toupper(sample.information$Gender)
sample.information$Group <- replace(sample.information$mnp11, is.na(sample.information$mnp11), "n.d.")
sample.information$Subgroup <- replace(sample.information$mnp12, is.na(sample.information$mnp11), "n.d.")

pdf("Plots/Cohort_overview.pdf", width=6, height=3.5)

p1 <- ggplot(sample.information, aes(x=Group, fill=Group)) + geom_bar() + scale_y_continuous(name="Number of Cases") +
  theme(axis.text.x = element_text(angle=90)) +  scale_fill_manual(values=group.colors) 

p2 <- ggplot(sample.information[sample.information$Group %in% c("MB, G3", "MB, G4"),], aes(x=Subgroup, fill=Subgroup)) + geom_bar() + scale_y_continuous(name="Number of Cases") +
  theme(axis.text.x = element_text(angle=90)) +  scale_fill_manual(values=g34.subgroup.colors) 

ggpubr::ggarrange(p1, p2, nrow=1, ncol=2)

dev.off()
