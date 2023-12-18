###########Yaun a la Wu

#Step 1: translate

#Step2: assign mutation and mutation types

#Step3: filter input (mutant plasmid) count >= 10

#Step4: filter input frequency >= 6*WT frequency.

# ***here we could add other error control: passaged frequency >= 6*WT_passaged frequency

#Step5: calculate enrichment for each mutation (and log10)

#Step6: subtract mean enrichment of silent mutations

library(tidyverse)

############################## Step 1-2 ###############################
codons <- tibble(
  codon = c("TTT","TTC",
            "TTA","TTG","CTT","CTC","CTA","CTG",
            "ATT","ATC","ATA",
            "ATG",
            "GTT","GTC","GTA","GTG",
            "TCT","TCC","TCA","TCG","AGT","AGC",
            "CCT","CCC","CCA","CCG",
            "ACT","ACC","ACA","ACG",
            "GCT","GCC","GCA","GCG",
            "TAT","TAC",
            "TAA","TAG","TGA",
            "CAT","CAC",
            "CAA","CAG",
            "AAT","AAC",
            "AAA","AAG",
            "GAT","GAC",
            "GAA","GAG",
            "TGT","TGC",
            "TGG",
            "CGT","CGC","CGA","CGG","AGA","AGG",
            "GGT","GGC","GGA","GGG"),
  aa = c("Phe","Phe",
         "Leu","Leu","Leu","Leu","Leu","Leu",
         "Ile","Ile","Ile",
         "Met",
         "Val","Val","Val","Val",
         "Ser","Ser","Ser","Ser","Ser","Ser",
         "Pro","Pro","Pro","Pro",
         "Thr","Thr","Thr","Thr",
         "Ala","Ala","Ala","Ala",
         "Tyr","Tyr",
         "***","***","***",
         "His","His",
         "Gln","Gln",
         "Asn","Asn",
         "Lys","Lys",
         "Asp","Asp",
         "Glu","Glu",
         "Cys","Cys",
         "Trp",
         "Arg","Arg","Arg","Arg","Arg","Arg",
         "Gly","Gly","Gly","Gly")
)

read_counts <- function(filename){
  counts <- read_csv(filename) %>%
    #translate by joining to codons dataframe
    pivot_longer(cols = -c(site, wildtype),
                 names_to = "codon", values_to = "count") %>%
    left_join(codons, by = "codon") %>%
    left_join(codons, by = c("wildtype" = "codon")) %>%
    #assign mutation types
    mutate(mutation = paste0(aa.y, site, aa.x, sep = ""),
           mutation_type = case_when(
             aa.x == aa.y & codon == wildtype ~ "wildtype",
             aa.x == aa.y ~ "silent",
             aa.x == "***" ~ "nonsense",
             TRUE ~ "missense"
           )) %>%
    #collapse mutations
    group_by(mutation, mutation_type) %>%
    summarize(count = sum(count))
  
  return(counts)
}

Mut_plasmid <- read_counts("C:/Users/dogcatfrogspider/Desktop/Desktop/py/20231122_AH/Mplasmid_codoncounts.csv")
WT_plasmid <- read_counts("C:/Users/dogcatfrogspider/Desktop/Desktop/py/20231122_AH/WTplasmid_codoncounts.csv")
passaged <- read_counts("C:/Users/dogcatfrogspider/Desktop/Desktop/py/20231122_AH/P1C_codoncounts.csv")

############################## Step 3-4 ###############################

#filter for counts >10 in plasmid libary
Input <- Mut_plasmid %>%
  filter(count >= 10)

input_subamplicons <- Input %>%
  mutate(site = as.numeric(str_extract(mutation, "[0-9]+"))) %>%
  mutate(subamplicon = case_when(
    site >= 1 & site <= 90 ~ "1",
    site >= 92 & site <= 187 ~ "2",
    site >= 188 & site <= 284 ~ "3",
    site >= 285 & site <= 382 ~ "4",
    site >= 383 & site <= 479 ~ "5",
    site >= 481 & site <= 578 ~ "6",
    site >= 579 & site <= 675 ~ "7",
    site >= 676 & site <= 758 ~ "8",
    TRUE ~ "invalid"
  ))

#combine all 3 counts datasets into a single dataframe
combined <- input_subamplicons %>%
  left_join(WT_plasmid, by = c("mutation", "mutation_type"),
            suffix = c(".input", ".wt")) %>%
  left_join(passaged, by = c("mutation", "mutation_type")) %>%
  rename("count.passaged" = count) %>%
  ungroup() %>%
  group_by(subamplicon) %>%
  #calculate frequencies (with a pseudocount added)
  mutate(freq.input = (count.input+1)/(sum(count.input) + length(count.input)),
         freq.wt = (count.wt+1)/(sum(count.wt) + length(count.wt)),
         freq.passaged = (count.passaged+1)/(sum(count.passaged) + length(count.passaged)),
         length_test = length(count.passaged)) %>%
  #filter for frequency in input >= 6 * frequency in wt plasmid
  filter(freq.input >= 6*freq.wt |
           mutation_type == "wildtype")


############################## Step 5 ###############################

#calculate enrichment from **frequencies** 
#note that this is different than wu paper, where they used counts
#but I think frequencies makes more sense since it would account for depth
enrichment <- combined %>%
  mutate(enrichment = (freq.passaged)/(freq.input))

############################## Step 6 ###############################

#calculate mean enrichment of silent mutations
silent_mean <- enrichment %>%
  filter(mutation_type == "silent") %>%
  pull(enrichment) %>%
  log10() %>%
  mean()

#Calculate fitness as log10 of enrichment - mean of silent mutations
fitness <- enrichment %>%
  mutate(fitness = ifelse(mutation_type == "wildtype", log10(enrichment),
                        log10(enrichment) - silent_mean))


############################## Plot ###############################

ggplot(fitness, aes(x = fitness, group = mutation_type,
                    color = mutation_type,
                    fill = mutation_type)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  facet_wrap(~mutation_type, scales = "free_y",
             ncol = 1)

fitness_plotting <- fitness %>%
  mutate(site = as.numeric(str_extract(mutation, "[0-9]+")))

ggplot(fitness_plotting, aes(x = site, y = fitness, color = mutation_type)) +
  geom_point(alpha = 0.7)+
  facet_wrap(~mutation_type,
             ncol = 1) +
  labs(title = "P1C Mutation class fitness by site")

####################### Normalize by subamplicon #####################


subamplicon_enrichment <- enrichment %>%
  mutate(site = as.numeric(str_extract(mutation, "[0-9]+"))) %>%
  mutate(subamplicon = case_when(
    site >= 1 & site <= 90 ~ "1",
    site >= 92 & site <= 187 ~ "2",
    site >= 188 & site <= 284 ~ "3",
    site >= 285 & site <= 382 ~ "4",
    site >= 383 & site <= 479 ~ "5",
    site >= 481 & site <= 578 ~ "6",
    site >= 579 & site <= 675 ~ "7",
    site >= 676 & site <= 758 ~ "8",
    TRUE ~ "invalid"
    ))

silent_fitness <- enrichment %>%
  filter(mutation_type == "silent") %>%
  mutate(fitness = log10(enrichment)) %>%
  group_by(subamplicon) %>%
  summarize(mean_silent_fitness = mean(fitness))

subamplicon_fitness <- enrichment %>%
  left_join(silent_fitness, by = "subamplicon") %>%
  #mutate(normalized_fitness = log10(enrichment) - mean_silent_fitness) %>%
  mutate(normalized_fitness = ifelse(mutation_type == "wildtype", log10(enrichment),
                                     log10(enrichment) - mean_silent_fitness))

palette <- c("#5C4B51", "#1A3A48", "#8CBEB2", "#048789", "#F2EBBF", #using Yuan's paper's palette
             "#F3B562", "#F06060", "#96505B", "#7E8AA2")

ggplot(subamplicon_fitness, aes(x = normalized_fitness, y = after_stat(scaled),
                                group = mutation_type,
                                color = mutation_type,
                                fill = mutation_type)) +
  geom_density(alpha = 0.5, linewidth = 1.2) +
  facet_wrap(~mutation_type,
             ncol = 1) +
  scale_color_manual(name = "Mutation Types",
                     values=c(palette[8], palette[6], palette[4], palette[9]),
                     breaks=c("missense", "nonsense", "silent", "wildtype"),
                     labels=c("Missense", "Nonsense", "Silent", "Wildtype")) +
  scale_fill_manual(name = "Mutation Types",
                    values=c(palette[8], palette[6], palette[4], palette[9]),
                    breaks=c("missense", "nonsense", "silent", "wildtype"),
                    labels=c("Missense", "Nonsense", "Silent", "Wildtype")) +
  labs(x = "Fitness", y = "Scaled density", title = "P1C Mutation class density and fitness ") +
  scale_x_continuous(limits = c(-4, 2), breaks = -4:2) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(subamplicon_fitness, aes(x = site, y = normalized_fitness, color = mutation_type)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~mutation_type, scales = "free_y",
             ncol = 1)

###Troubleshooting

ggplot(subamplicon_fitness, aes(x = site, y = normalized_fitness, color = mutation_type)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~mutation_type,
             ncol = 1)

ggplot(fitness_plotting, aes(x = site, y = fitness, color = mutation_type)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~mutation_type,
             ncol = 1)

#Compare normalizing to subamplicon vs whole data

normalization_comparison <- fitness %>%
  left_join(subamplicon_fitness, by = c("mutation", "mutation_type"))

ggplot(normalization_comparison, aes(x = fitness, y = normalized_fitness)) +
  geom_point(alpha = 0.2) +
  labs(x = "nonnormalized fitness", y ="fitness", title = "P1C Fitness comparison to nonnormalized")

####Hotspot investigation####

hotspot <- combined %>%
  mutate(site = as.numeric(str_extract(mutation, "[0-9]+"))) %>%
  filter(mutation_type != "wildtype") %>%
  group_by(site) %>%
  summarize(total_mutfreq = sum(freq.passaged))

fitness_200zoom <- fitness_plotting %>%
  filter(site %in% c(190:210))

ggplot(fitness_200zoom, aes(x = site, y = fitness, color = mutation_type)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~mutation_type,
             ncol = 1)

passaged_ntchecks <- read_csv("C:/Users/dogcatfrogspider/Desktop/Desktop/py/20231122_AH/P1C_codoncounts.csv") %>%
  filter(site == 206) %>%
  pivot_longer(cols = -c("site","wildtype"),
               names_to = "codon", values_to = "count")%>%
  arrange(desc(count))
  

  
