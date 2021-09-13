
## load libraries

library ("dplyr")
library("ggplot2")
library("eulerr")
library("ggvenn")
library("stringr")


## Read args from command line
args = commandArgs(trailingOnly=TRUE)

## Uncomment For debugging only
## Comment for production mode only
#args[1] <-"test/data/sample.targets.ref"

#args[2] <- "test/data/sample.targets.mut"

#args[3] <- "test/data/sample.changes" # output file

  ## get the mirmap tsv file from args
mirna_ref_file <- args[1]

## get targetscan tsv file from args
mirna_mut_file <- args[2]

## pass to named objects
mirna_changes <- args[3]

## Read miRNA targets
mirna_ref.df <- read.table(file= mirna_ref_file, header = T,
                           sep = "\t", stringsAsFactors = FALSE)

mirna_mut.df <- read.table(file= mirna_mut_file, header = T,
                          sep = "\t", stringsAsFactors = FALSE)


## Select mirnas targets predicted by both tools
mirna_ref_intersect.df <- mirna_ref.df %>% filter(prediction_tool ==  "both") %>% select(miRNA_ID)

mirna_mut_intersect.df <- mirna_mut.df %>% filter(prediction_tool ==  "both") %>% select(miRNA_ID)

## Select mirnas targets predicted by any tool
mirna_ref_all.df <- mirna_ref.df %>% select(miRNA_ID)

mirna_mut_all.df <- mirna_mut.df %>% select(miRNA_ID)

## Get Lost target mirna pairs
lost_targets <- mirna_ref_intersect.df %>% setdiff(mirna_mut_all.df)

## Get Gain target mirna pairs
gain_targets <- mirna_mut_intersect.df %>% setdiff(mirna_ref_all.df)

## Get remained targets r
remained_targets.df <-mirna_ref_intersect.df %>%  intersect(mirna_mut_intersect.df)


## Define if one target is lost, gained o remained
lost_targets <- lost_targets %>%  mutate(target = "lost")
gain_targets <- gain_targets  %>%  mutate(target = "gained")
remained_targets.df <- remained_targets.df  %>%  mutate(target = "remained")

## Merge the miRNA targets gained and lost into a single dataframe
target_changes.df <- full_join(x = lost_targets, y = gain_targets,
                               by = c("miRNA_ID", "target") )

## Merge all miRNA targets ids into a single dataframe
All_targets.df <- full_join(x = target_changes.df, y = remained_targets.df,
                        by = c("miRNA_ID", "target") )

## Save dataframe
write.table(All_targets.df, file = mirna_changes, sep = "\t", na = "NA", quote = F, row.names = F)


## Make a vector with the mirnas targets predicted by both tools
mirna_ref_intersect.v <- mirna_ref_intersect.df %>%  pull(1) %>%  unique()

mirna_mut_intersect.v <- mirna_mut_intersect.df %>%  pull(1) %>%  unique()

## Sort the ids list within a list for ggvenn
Venn_list <- list(
  A = mirna_ref_intersect.v,
  B = mirna_mut_intersect.v)

## Name the source of the ids
names(Venn_list) <- c("miRNAs REF","miRNAs MUT")

## Ṕlot a Venn diagram
miRNAs_Venn.p <- ggvenn(Venn_list, fill_color = c("#FF595E", "#007F5F"),
                        stroke_size = 0.5, set_name_size = 4 , text_size = 4)


## Save plot
ggsave( filename = str_interp("${mirna_changes}.png"),
        plot = miRNAs_Venn.p,
        device = "png",
        height = 7, width = 14,
        units = "in")

## Make eulerr plot
microRNAs_euler <- euler(Venn_list)

microRNAs_euler.p <- plot( x = microRNAs_euler,
                           quantities = TRUE,               
                           main = "microRNAS iDs",
                           fill = c("#FF595E", "#007F5F") )                 


# save plot
ggsave( filename = str_interp("${mirna_changes}_2.png"),        
        plot = microRNAs_euler.p,                
        device = "png",                 
        height = 7,                     
        width = 14,
        units = "in",
        dpi = 300 )                    


## Select mirnas targets predicted by TargetScan
mirna_ref_targetscan.df <- mirna_ref.df %>% filter(prediction_tool ==  "targetscan" | 
                                                     prediction_tool == "both") %>% select(miRNA_ID)

mirna_mut_targetscan.df <- mirna_mut.df %>% filter(prediction_tool ==  "targetscan" |
                                                     prediction_tool == "both") %>% select(miRNA_ID)


## Select mirnas targets predicted by mirmap
mirna_ref_mirmap.df <- mirna_ref.df %>% filter(prediction_tool ==  "mirmap" | 
                                                 prediction_tool == "both") %>% select(miRNA_ID)

mirna_mut_mirmap.df <- mirna_mut.df %>% filter(prediction_tool ==  "mirmap" |
                                                 prediction_tool == "both") %>% select(miRNA_ID)



## Make a vector with the mirnas targets predicted by both tools
mirna_ref_targetscan.v <- mirna_ref_targetscan.df %>%  pull(1) %>%  unique()

mirna_ref_mirmap.v <- mirna_ref_mirmap.df %>%  pull(1) %>%  unique()

mirna_mut_targetscan.v <- mirna_mut_targetscan.df %>%  pull(1) %>%  unique()

mirna_mut_mirmap.v <- mirna_mut_mirmap.df %>%  pull(1) %>%  unique()

## Sort the ids list within a list for ggvenn
Venn_list <- list(
  A = mirna_ref_targetscan.v,
  B = mirna_ref_mirmap.v,
  C = mirna_mut_targetscan.v,
  D = mirna_mut_mirmap.v)

## Name the source of the ids
names(Venn_list) <- c("REF_TargetScan","REF_miRmap","MUT_TargetScan","MUT_miRmap")

## Ṕlot a Venn diagram
miRNAs_Venn.p <- ggvenn(Venn_list, fill_color = c("#D9ED92", "#99D98C", "#168AAD", "#1E6091"),
                        stroke_size = 0.5, set_name_size = 4 , text_size = 4)


## Save plot
ggsave( filename = str_interp("${mirna_changes}_3.png"),
        plot = miRNAs_Venn.p,
        device = "png",
        height = 7, width = 15,
        units = "in")

