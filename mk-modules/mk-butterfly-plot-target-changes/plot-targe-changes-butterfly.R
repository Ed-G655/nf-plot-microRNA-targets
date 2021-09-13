## load libraries
library ("dplyr")
library("ggplot2")
library("ggvenn")
library("stringr")
library("cowplot")
library("tidyr")
## Read args from command line
args = commandArgs(trailingOnly=TRUE)

## Uncomment For debugging only
## Comment for production mode only
#args[1] <-"test/data/sample.changes.tsv"

#args[2] <- "test/data/sample1.png" # output file

## get the file with the miRNA targets and its changes
mirna_changes_file <- args[1]

## pass to named objects
mirna_plot <- args[2]

## Read miRNA targets
mirna_changes.df <- read.table(file= mirna_changes_file, header = T,
                           sep = "\t", stringsAsFactors = FALSE)

number_of_targets <- mirna_changes.df %>%  count(miRbase_ID)

lost_targets <- mirna_changes.df %>%  filter(target_change == "lost")

number_of_lost <- lost_targets %>%  count(miRbase_ID)

number_of_lost <- number_of_lost %>% mutate(target_change = "lost")

number_of_lost <- number_of_lost %>% mutate(n = n*-1) 

#number_of_lost_genes <- lost_targets %>%  count(GeneID)

gain_targets <- mirna_changes.df %>%  filter(target_change == "gained")

number_of_gain <- gain_targets %>%  count(miRbase_ID)

number_of_gain <- number_of_gain %>% mutate(target_change = "gained")

#number_of_gain_genes <- gain_targets %>%  count(GeneID)

remained_targets <- mirna_changes.df %>%  filter(target_change == "remained")

number_of_remained <- remained_targets %>%  count(miRbase_ID)

number_of_remained <- number_of_remained %>% mutate(target_change = "remained")

count_change_target_1 <-  number_of_lost %>%  bind_rows(number_of_gain) %>% bind_rows(number_of_remained)

#count_change_target_1 <- left_join(x = number_of_lost, y = number_of_gain, by = "miRbase_ID")

names(count_change_target_1)[2] <- "Number_of_targets"

count_change_target_1 <- count_change_target_1 %>%  mutate( miRbase_ID = (miRbase_ID %>% str_replace(":.*","")))

#names(count_change_target_1)[3] <- "Gain_targets"
number_of_targets <- number_of_targets %>%  mutate(miRbase_ID = number_of_targets[,1] %>%  str_replace(":.*",""))

names(number_of_targets)[2] <- "Number_of_Targets"

count_change_target_1 <- arrange(count_change_target_1, Number_of_targets )

paleta <- c("lost" =  "#F94144",
            "gained" = "springgreen3") 

min <- count_change_target_1$Number_of_targets %>% min()
max <- count_change_target_1$Number_of_targets %>% max()

#Plot_gain_and_lost
piramide.p <- ggplot(count_change_target_1, aes(x = miRbase_ID, y = Number_of_targets, fill = target_change )) + 
  geom_col(data = subset(count_change_target_1, target_change == "lost"), 
           width = 0.5, fill = "#F94144") + 
  geom_col(data = subset(count_change_target_1, target_change ==  "gained"), 
           width = 0.5, fill = "springgreen3") +
  coord_flip() + scale_y_continuous(
    breaks = c(seq(min, -0, by = 500), 
               seq(0, max, by = 500)),
    labels = c(seq(min, 0, by = 500)*-1, 
               seq(0, max, by = 500))
  ) + labs(y= "Numero de pares miRNA/blanco", x = "miRbase ID", color = "Legend") +
  scale_color_manual(values = paleta) +
  labs(title = "Sitos blanco por miRNA y sus cambios debido a mutaciones en el miRNA") +
  theme_minimal()

ggsave( filename = mirna_plot , 
        plot = piramide.p,
        device = "png",
        height = 7, width = 14,
        units = "in")

