
# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(scales)
library(gridExtra)
library(grid)
library(cowplot)

# Load the frames done in advance that are in a good structure for plotting
steps_per_frag <- read_tsv("walks_tsvs/walk_summary_ns.tsv")
permutated_steps_per_frag <- read_tsv("walks_tsvs/permutated_walk_summary_ns.tsv")

# Create object just with walks that are of class 2 and have one or more inter-chromosomal interaction
walks_with_one_or_more_interchr_steps <- steps_per_frag %>%
  ungroup() %>%
  filter(class == 2, inter_chr >= 1)

# Get the number of intra-chromosomal walks and fraction of cis steps for the permutations
permutations_info_per_step_by_state <- permutated_steps_per_frag %>%
  separate(walk_id, c("pacbio_frag_ID", "iteration_clone"), sep = " ") %>%
  filter(class == 2, pacbio_frag_ID %in% walks_with_one_or_more_interchr_steps$pacbio_frag_ID) %>%
  group_by(iteration, state, steps) %>%
  summarize(intra_chromosome = sum(num_chr == 1)/length(num_chr),
            fraction_cis_hops = 1 - (sum(inter_chr)/sum(steps)))

# Get the number of intra-chromosomal walks and fraction of cis steps for the real walks
info_per_step_by_state <- steps_per_frag %>%
  filter(class == 2, pacbio_frag_ID %in% walks_with_one_or_more_interchr_steps$pacbio_frag_ID) %>%
  group_by(state, steps) %>%
  summarize(intra_chromosome = sum(num_chr == 1)/length(num_chr),
            fraction_cis_hops = 1 - (sum(inter_chr)/sum(steps)))

# Do the first plot, of the number of inter-chrosomomal interactions per walk
p1 <- steps_per_frag %>%
  ungroup() %>%
  filter(class == 2, state == "ns", steps <= 10) %>%
  ggplot(aes(x = steps, fill = as.factor(inter_chr))) +
  geom_bar() +
  scale_fill_viridis_d("Number of\ntrans interactions") +
  scale_x_continuous(breaks = 1:10) +
  theme_classic()
print(p1)

# Do the second plot, of the fraction of cis steps in the real and permutated walks
# Subset for walks with three or more steps
p2 <- ggplot(filter(permutations_info_per_step_by_state, steps >= 3, steps <= 10, state == "ns"),
             aes(x=steps, y=fraction_cis_hops)) +
  geom_jitter(size = 1, alpha = I(1/4) ) +
  labs(y = "Fraction of cis steps") +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous() +
  geom_line(data = filter(info_per_step_by_state, steps >= 3, steps <= 10, state == "ns"),
            size = 1) +
  theme_classic() +
  labs(title = "Walks with one or more\ninterchromosomal steps")
print(p2)

# Print the two plots together
cowplot::plot_grid(p1, ggplot() + theme_classic(), p2,
                   align = "h", axis = "l",
                   rel_widths = c(0.45, 0.1, 0.45),
                   ncol = 3)
