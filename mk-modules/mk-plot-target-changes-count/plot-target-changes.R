## load libraries
library ("dplyr")
library("ggplot2")
library("ggvenn")
library("stringr")
library("cowplot")

## Read args from command line
args = commandArgs(trailingOnly=TRUE)

## Uncomment For debugging only
## Comment for production mode only
#args[1] <-"test/data/sample.changes.tsv"

#args[2] <- "test/data/sample.png" # output file

## get the file with the miRNA targets and its changes
mirna_changes_file <- args[1]

## get targetscan tsv file from args
mirna_mut_file <- args[2]

## Read miRNA targets
mirna_changes.df <- read.table(file= mirna_changes_file, header = T,
                           sep = "\t", stringsAsFactors = FALSE)
#Count targets
number_of_targets <- mirna_changes.df %>%  count(miRbase_ID)

lost_targets <- mirna_changes.df %>%  filter(target_change == "lost")

number_of_lost <- lost_targets %>%  count(miRbase_ID)

number_of_lost <- number_of_lost %>% mutate(target_change = "lost")

gain_targets <- mirna_changes.df %>%  filter(target_change == "gained")

number_of_gain <- gain_targets %>%  count(miRbase_ID)

number_of_gain <- number_of_gain %>% mutate(target_change = "gained")

remained_targets <- mirna_changes.df %>%  filter(target_change == "remained")

number_of_remained <- remained_targets %>%  count(miRbase_ID)

number_of_remained <- number_of_remained %>% mutate(target_change = "remained")

count_change_target_1 <-  number_of_lost %>%  bind_rows(number_of_gain) %>% bind_rows(number_of_remained)

names(count_change_target_1)[2] <- "Number_of_targets"

count_change_target_1 <- count_change_target_1 %>%  
  mutate( miRbase_ID = (miRbase_ID %>% str_replace(":.*","")))

count_change_target_2 <- full_join(x = count_change_target_1,
                                   y = number_of_remained, by = "miRbase_ID" ,"n")

number_of_targets <- number_of_targets %>%  
  mutate(miRbase_ID = number_of_targets[,1] %>%  str_replace(":.*",""))

names(number_of_targets)[2] <- "Number_of_Targets"

# plot gain, lost and remain targets  

gain_and_lost.p <- ggplot(count_change_target_1, aes(x = miRbase_ID, 
                                                     y = Number_of_targets, 
                                                     fill = target_change)) + 
  geom_bar(position = "stack", stat = "identity") + theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 0, vjust = 0 )) + 
  scale_y_continuous(expand = c(0,0)) +
  ylab("Numero de pares miRNA/blanco") +
  xlab("ID miRbase") +
  labs(title = "Sitos blancos por miRNA y sus cambios debido a mutaciones en el miRNA") +
  labs(fill="Target change")
  

ggsave( filename = mirna_mut_file,
        plot = gain_and_lost.p,
        device = "png",
        height = 7, width = 14,
        units = "in")
