########################################################
### SHAMROCK clustering and plotting            ~~~~ ###
### VERSION: 0.1.0                              ~~~~ ###
### LAST EDIT: 10/06/25                         ~~~~ ###
### AUTHORS: Richard Edwards 2025               ~~~~ ###
### CONTACT: rich.edwards@uwa.edu.au            ~~~~ ###
### CITE: https://github.com/slimsuite/shamrock ~~~~ ###
########################################################

# This script accompanies shamrock.sh

####################################### ::: HISTORY ::: ############################################
# v0.0.0 : Initial version with basic plotting.
# v0.1.0 : Added partition=INT to set number of subgenomes [2].
version = "v0.1.0"

####################################### ::: USAGE ::: ############################################
# Example use:
# Rscript shamrock.R [basefile=FILE]
# : basefile=FILE = Prefix for main inputs and outputs []
# : k=INT = kmer length used
# : partition=INT = number of subgenomes to split into [2]
# : debug=T/F = whether to switch on additional debugging outputs [FALSE]
# : dev=T/F = whether to switch on dev mode during code updates [FALSE]
# : outlog=FILE = outlog log messages to $FILE [stdout()]

#i# Usage within R:
# Set an override vector of commandline arguments: this will replace argvec read from the commandline
# Use source() to run the script:
# source("$PATH/chromsyn.R")

####################################### ::: OUTPUTS ::: ############################################
#!# List the outputs here

####################################### ::: TODO ::: ############################################
# [ ] : Add update plan from README.md
# [ ] : Add dynamic sizing of the outputs
# [ ] : Add different colour palette options: shamrock/green/heat/inferno
# [ ] : Add more output, including output to log.
# [ ] : Update to call the shell script and perform additional calculations for assigning homeologues etc.
# [ ] : Consider wrapping up in Rmd instead, a bit like MerquryRising.

####################################### ::: SETUP ::: ############################################
### ~ Commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
defaults = list(k=31,partition=2,
                debug=FALSE,dev=FALSE,
                rdir="",runpath="",
                outlog=stdout())

settings <- defaults
argvec = commandArgs(TRUE)
if("override" %in% ls()){
  argvec = override
  settings$rscript = FALSE
}
for(cmd in argvec){
  cmdv = strsplit(cmd,'=',TRUE)[[1]]
  if(length(cmdv) > 1){
    settings[[cmdv[1]]] = cmdv[2]
  }else{
    settings[[cmdv[1]]] = ""    
  }
}
#i# integer parameters
for(cmd in c("k","partition")){
  settings[[cmd]] = as.integer(settings[[cmd]])
}
#i# other numeric parameter
# for(cmd in c("")){
#   settings[[cmd]] = as.numeric(settings[[cmd]])
# }
#i# adjust parameters where needed

#i# list parameters
# for(cmd in c("")){
#   if(sum(grep(",",settings[[cmd]],fixed=TRUE)) > 0){
#     settings[[cmd]] = strsplit(settings[[cmd]],',',TRUE)[[1]]
#   }else{
#     settings[[cmd]] <- settings[[cmd]]
#   }
# }
#i# logical parameters
for(cmd in c("debug","dev")){
  settings[[cmd]] = as.logical(settings[[cmd]])
}

#i# Set warnings based on debug
oldwarn <- getOption("warn")
if(settings$debug){
  writeLines(argvec)
}else{
  options(warn = -1)
}

### ~ logWrite function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
logWrite <- function(logstr){
  writeLines(paste0("[",date(),"] ",logstr),con=settings$outlog)
}
if(! settings$runpath == ""){
  setwd(settings$runpath)
}
logWrite(paste("#RCODE shamrock.R:",version))
logWrite(paste("#PATH Running from:",getwd()))
for(cmd in names(settings)[order(names(settings))]){
  logWrite(paste("CMD:",cmd,"=",paste0(settings[[cmd]],collapse=",")))
}
#dir.create(settings$plotdir, showWarnings = FALSE)

### ~ Load packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
if(! "tidyverse" %in% installed.packages()[,"Package"]){
  install.packages("tidyverse")
}
library(tidyverse)

####################################### ::: FUNCTIONS ::: ############################################

##### ======================== Loading data functions ======================== #####


logWrite('#RCODE Setup complete.')



####################################### ::: MAIN RUN CODE ::: ############################################

##### ======================== Load data ======================== #####


kfile <- paste0(settings$basefile,".allok",settings$k,".csv")
df <- read.csv(kfile, stringsAsFactors = FALSE)

df_self <- df %>% filter(chri == chrj) %>% select(chri,knum) %>% rename(selfk=knum)

df_plot <- df %>%
  select(chri, chrj, knum) %>% left_join(df_self) %>%
  mutate(knorm = knum/selfk)

df_max <- df_plot %>% filter(chri != chrj) %>% group_by(chri) %>% summarise(maxk=max(knorm))

df_plot <- df_plot %>% left_join(df_max) %>% mutate(kplot = knorm/maxk) %>%
  mutate(kplot=if_else(chri==chrj,1.0,kplot))

#!# Need to adapt this!
df_plot$knum[df_plot$knum == 0] <- min(df_plot$knum[df_plot$knum != 0])
df_plot$knum <- log10(df_plot$knum)

p <- ggplot(df_plot, aes(x = chrj, y = chri, fill = knum)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno") +
  theme_minimal() +
  labs(title = "Heatmap of raw knum values", x = "chrj", y = "chri") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
logWrite(paste("Saving raw heatmap to",paste0(settings$basefile,"rawk.pdf/png")))
pfile <- paste0(settings$basefile,".rawk.pdf")
ggsave(pfile,p)
pfile <- paste0(settings$basefile,".rawk.png")
ggsave(pfile,p)


# Pivot to wide format (raw matrix, chri as rows, chrj as columns)
mat_dist <- df_plot %>% mutate(kdist = 1 - kplot) %>%
  select(chri, chrj, kdist) %>%
  pivot_wider(names_from = chrj, values_from = kdist, values_fill = 0) %>%
  column_to_rownames("chri")

# Convert to matrix for heatmap
mat_dist <- as.matrix(mat_dist)

# Hierarchical clustering
hc_rows <- hclust(dist(mat_dist), method = "average")
hc_cols <- hclust(dist(t(mat_dist)), method = "average")


# Get row and column order
row_order <- hc_rows$labels[hc_rows$order]
col_order <- hc_cols$labels[hc_cols$order]

# Reorder factor levels
df_plot$chri <- factor(df_plot$chri, levels = row_order)
df_plot$chrj <- factor(df_plot$chrj, levels = col_order)

# Plot heatmap
df_plot$knorm <- df_plot$kplot
p <- ggplot(df_plot, aes(x = chrj, y = chri, fill = knorm)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno") +
  theme_minimal() +
  labs(title = "Clustered heatmap of normalised shared allokmers", x = "chrj", y = "chri") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#!# Add dynamic sizing of the outputs
logWrite(paste("Saving normalised heatmap to",paste0(settings$basefile,".shamrock.pdf/png")))
pfile <- paste0(settings$basefile,".shamrock.pdf")
ggsave(pfile,p)
pfile <- paste0(settings$basefile,".shamrock.png")
ggsave(pfile,p)

# Assign parents - group names by their assigned cluster
logWrite("Primary chromosome clustering...")
clusters <- cutree(hc_rows, k = settings$partition)
clustered_names <- split(names(clusters), clusters)
for (i in seq_along(clustered_names)) {
  logWrite(paste(sprintf("Parent %d:", i),paste(clustered_names[[i]],collapse=",")))
  writeLines(clustered_names[[i]], paste0(settings$basefile,".subgenome.",i,".txt"))
}

logWrite("Secondary chromosome clustering check...")
clusters <- cutree(hc_cols, k = settings$partition)
clustered_names <- split(names(clusters), clusters)
for (i in seq_along(clustered_names)) {
  logWrite(paste(sprintf("Parent %d:", i),paste(clustered_names[[i]],collapse=",")))
}
#!# Add checks that both sets of clusters are the same.



##### ======================== Tidy and Finish ======================== #####
if(file.exists("Rplots.pdf")){
  file.remove("Rplots.pdf")
}

### ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
options(warn = oldwarn)
logWrite("#RCODE shamrock.R finished.")
#quit("no",0)

