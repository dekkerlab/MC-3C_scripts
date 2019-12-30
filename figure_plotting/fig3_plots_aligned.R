
library(ggplot2)
library(dplyr)
library(tidyr)

library(ggpubr)

# Set global theme to classic.
theme_set(theme_classic())

# Fig 3A

load("frames_for_plotting_stricter_keep_nas.RData")

# Subset the frames so that we remove the cluster samples
direct_ints <- filter(direct_ints, state != "cluster")
steps_per_frag <- filter(steps_per_frag, state != "cluster")

# Plot the step distance for NS samples
# Separate between same domain and different domain steps
# Only intra-chromosomal walks
p <- ggplot(filter(left_join(direct_ints,
                             select(ungroup(steps_per_frag), pacbio_frag_ID, compartment)),
                   state == "ns",
                   is.na(first_frag_compartment_index) == FALSE, is.na(second_frag_compartment_index) == FALSE,
                   class == 1),
            aes(x = abs(dist), y = ..density..)) +
  geom_density(mapping = aes(y = ..density.. * 2), colour = "grey") +
  geom_density(mapping = aes(colour = first_frag_compartment_index == second_frag_compartment_index)) +
  scale_x_log10("Step size (bp)") +
  scale_colour_manual("Interactions between\nfragments in:",
                      labels = c("different compartment domains", "same compartment domain"), values = c("red", "blue")) +
  theme_classic()
print(p)

fig3a <- p


# Fig 3B

load("frames_for_plotting_stricter.RData")
load("frames_for_plotting_simulations.RData")

# Plot number of full intra-compartment walks in real walks and in Hi-C based simulations
p <- ggplot(mapping = aes(x=steps, linetype = type)) +
  geom_line(data = filter(steps_per_frag,
                          state == "ns",
                          number_of_compartment_types < 3,
                          steps <= 7,
                          steps > 1,
                          class == 1,
                          is.na(compartment) == FALSE) %>%
              ungroup() %>%
              count(state, number_of_compartment_types, steps) %>%
              group_by(state, steps)  %>%
              mutate(sum = sum(n), type = "Real walks") %>%
              filter(number_of_compartment_types == 1), mapping = aes( y = n / sum)) +
  geom_line(data = filter(sim_steps_per_frag,
                          state == "ns",
                          number_of_compartment_types < 3,
                          chr == "chr4",
                          steps <= 7,
                          steps > 1,
                          class == 1,
                          is.na(compartment)== FALSE) %>%
              count(state, number_of_compartment_types, steps) %>%
              group_by(state, steps) %>%
              mutate(sum = sum(n), type = "Simulated walks") %>%
              filter(number_of_compartment_types == 1),
            mapping = aes( y = n / sum), alpha = I(1/2)) +
  scale_y_continuous("Proportion of walks\ninside one compartment type") +
  scale_colour_discrete("", labels = c("One compartment type", "Both compartment types")) +
  theme_classic()
print(p)

fig3b <- p

# New figure 3B

#Also, I think we need to rethink figure 3B:
#  that figure needs to show that the walks are enriched in A-A and B-B interactions.
# THe simulation you show now does not really show that in the right way. 
#
# Can we show: for all walks that touch 2 or more domains:
#   more A-A and B-B steps than expected
#  (based on number of fragments that are in A and B for the set of walks),
#   fewer A-B steps?

load("frames_for_plotting_stricter_keep_nas.RData")

unique_frags_per_comp_type <- unique(rbind(direct_ints %>%
                                             ungroup() %>%
                                             filter(#pacbio_frag_ID %in% walks_two_or_more_distinct_compartments$pacbio_frag_ID,
                                               state == "ns",
                                               class == 1) %>%
                                             select(V1, V2, V3, first_frag_compartment_type, first_frag_compartment_index),
                                           direct_ints %>%
                                             ungroup() %>%
                                             filter(#pacbio_frag_ID %in% walks_two_or_more_distinct_compartments$pacbio_frag_ID,
                                               state == "ns",
                                               class == 1) %>%
                                             select(V1 = V4, V2 = V5, V3 = V6,
                                                    first_frag_compartment_type = second_frag_compartment_type,
                                                    first_frag_compartment_index = second_frag_compartment_index))) %>%
  count(first_frag_compartment_type) %>%
  drop_na()

a_frags <- filter(unique_frags_per_comp_type, first_frag_compartment_type == "plus")$n
b_frags <- filter(unique_frags_per_comp_type, first_frag_compartment_type == "minus")$n
total_frags <- a_frags + b_frags

exp_prob_aa <- (a_frags / total_frags) * (a_frags / total_frags)
exp_prob_bb <- (b_frags / total_frags) * (b_frags / total_frags)
exp_prob_ab <- (a_frags / total_frags) * (b_frags / total_frags)
exp_prob_ba <- (b_frags / total_frags) * (a_frags / total_frags)

exp_probs <- tibble(first_frag_compartment_type = c("minus", "minus", "plus", "plus"),
                    second_frag_compartment_type = c("minus", "plus", "minus", "plus"),
                    exp_prob = c(exp_prob_bb, exp_prob_ba, exp_prob_ab, exp_prob_aa)) %>%
  mutate(interaction_class = ifelse(first_frag_compartment_type == second_frag_compartment_type,
                                    first_frag_compartment_type,
                                    "AB"),
         exp_prob = ifelse(interaction_class == "AB",
                           exp_prob,
                           exp_prob)) %>%
  select(interaction_class, exp_prob) %>%
  unique()

observed_probs <- direct_ints %>%
  ungroup() %>%
  filter(#pacbio_frag_ID %in% walks_two_or_more_distinct_compartments$pacbio_frag_ID,
    state == "ns",
    class == 1) %>%
  count(first_frag_compartment_type, second_frag_compartment_type) %>%
  drop_na() %>%
  mutate(observed_prob = n / sum(n)) %>%
  select(first_frag_compartment_type,
         second_frag_compartment_type,
         observed_prob) %>%
  mutate(interaction_class = ifelse(first_frag_compartment_type == second_frag_compartment_type,
                                    first_frag_compartment_type,
                                    "AB"),
         observed_prob = ifelse(interaction_class == "AB",
                                observed_prob,
                                observed_prob)) %>%
  select(interaction_class, observed_prob) %>%
  unique()

new_fig3b <- left_join(exp_probs, observed_probs) %>%
  group_by(interaction_class) %>%
  summarise(exp_prob = sum(exp_prob),
            observed_prob = sum(observed_prob)) %>%
  ggplot(aes(x = interaction_class, y = log2(observed_prob / exp_prob))) +
  geom_col() +
  scale_x_discrete("Interaction", labels = c("AB", "BB", "AA")) +
  scale_y_continuous(expression(log2("Observed probability" / "Expected probability")))

new_fig3b



# Fig 3C

# Load the frames done in advance that are in a good structure for plotting
load("frames_for_plotting_stricter.RData")
load("frames_for_plotting_permutated.RData")

# Subset the frames so that we remove the cluster samples
direct_ints <- filter(direct_ints, state != "cluster")
steps_per_frag <- filter(steps_per_frag, state != "cluster")

supp.labs <- c("A-B", "A-A", "B-B")
names(supp.labs) <- c("non_unique_compartment", "plus", "minus")

# Number of compartment domains per walk
# For walks that have two or more domains (akin to interchromosomal walks in figure 2)
# Walks must have at least one inter-compartment domain interaction
# Separate A-B, A and B walks
p <- steps_per_frag %>%
  ungroup() %>%
  filter(state == "ns",
         steps <= 10,
         class == 1,
         is.na(compartment) == FALSE,
         number_of_compartments >= 2,
         different_compartment_index_ints_number > 0) %>%
  mutate(number_of_compartments_for_label = ifelse(number_of_compartments >= 5,
                                                   "5+", as.character(number_of_compartments))) %>%
  ggplot(aes(x = steps,
             fill = factor(number_of_compartments_for_label,
                              levels = c("2", "3", "4", "5+"),
                              ordered = TRUE))) +
  geom_bar() +
  scale_fill_viridis_d(name = "Number of\ncompartment\ndomains") +
  facet_wrap(~ compartment, labeller = labeller(compartment = supp.labs)) +
  scale_x_continuous(breaks = 1:10, labels = 1:10) +
  theme_classic()

print(p)

fig3c <- p

# Fig 3D

# Number of inter-compartment domain interactions per walk
# Separete for walks that visit only two domains and that visit two or more domains
# Walks must have at least one inter-compartment domain interaction
p <- rbind(steps_per_frag %>%
             ungroup() %>%
             filter(state == "ns",
                    steps <= 10,
                    class == 1,
                    is.na(compartment) == FALSE,
                    number_of_compartments == 2,
                    different_compartment_index_ints_number > 0) %>%
             mutate(set = factor("Only two domains",
                                 levels = c("Two or more domains", "Only two domains"),
                                 ordered = TRUE)),
           steps_per_frag %>%
             ungroup() %>%
             filter(state == "ns",
                    steps <= 10,
                    class == 1,
                    is.na(compartment) == FALSE,
                    number_of_compartments >= 2,
                    different_compartment_index_ints_number > 0) %>%
             mutate(set = factor("Two or more domains",
                                 levels = c("Two or more domains", "Only two domains"),
                                 ordered = TRUE))) %>%
  mutate(different_compartment_index_ints_number_for_label = ifelse(different_compartment_index_ints_number >= 5,
                                                                    "5+", as.character(different_compartment_index_ints_number))) %>%
  ggplot(aes(x = steps, fill = factor(different_compartment_index_ints_number_for_label,
                                      levels = c("1", "2", "3", "4", "5+"),
                                      ordered = TRUE))) +
  geom_bar() +
  scale_fill_viridis_d(name = "Number of\ninteractions\nbetween\ncompartment\ndomains") +
  scale_x_continuous(breaks = 1:10, labels = 1:10) +
  facet_wrap(~ set, scales = "free_y") +
  theme_classic()

print(p)

fig3d <- p


# Fig 3E

# Select walks that visit two or more distinct compartment domains
walks_two_or_more_distinct_compartments <- steps_per_frag %>%
  ungroup() %>%
  filter(state == "ns",
         steps <= 10,
         steps >= 3,
         class == 1,
         number_of_compartments >= 2,
         different_compartment_index_ints_number > 0)

# Plot fraction of steps within the same compartment domain
# For real walks and permutations
p <- ggplot(mapping = aes(x = steps, y = prop_same_comp, alpha = set)) +
  geom_line(data = walks_two_or_more_distinct_compartments %>%
              filter(steps <= 6) %>%
              group_by(state, steps) %>%
              summarise(prop_same_comp = 1 - sum(different_compartment_index_ints_number)/sum(steps)) %>%
              mutate(set = "Real walk")
  ) +
  geom_jitter(data = permutated_steps_per_frag %>%
                ungroup() %>%
                mutate(compartment = factor(compartment, levels = c("non_unique_compartment",
                                                                    "plus", "minus"),
                                            ordered = TRUE)) %>%
                separate(walk_id, into = c("walk_id", "iteration_clone"), sep = " ") %>%
                filter(walk_id %in% walks_two_or_more_distinct_compartments$pacbio_frag_ID) %>%
                filter(steps <= 6) %>%
                group_by(iteration, state, steps) %>%
                summarise(prop_same_comp = 1 - sum(different_compartment_index_ints_number)/sum(steps)) %>%
                mutate(set = "Permutated walk")) +
  theme_classic() +
  facet_wrap(~ "All such walks") +
  ggtitle("") +
  scale_x_continuous(breaks = 3:10, labels = 3:10) +
  scale_y_continuous("Fraction of intra-domain steps", limits = c(0.2, 0.7)) +
  scale_alpha_discrete("Set of walks")

print(p)

fig3ea <- p

# Same as above, but separate A-B, A-A and A-B walks
p <- ggplot(mapping = aes(x = steps, y = prop_same_comp, alpha = set)) +
  geom_line(data = walks_two_or_more_distinct_compartments %>%
              filter(steps <= 6) %>%
              group_by(state, compartment, steps) %>%
              summarise(prop_same_comp = 1- sum(different_compartment_index_ints_number)/sum(steps)) %>%
              mutate(set = "Real walk")
            ) +
  geom_jitter(data = permutated_steps_per_frag %>%
                ungroup() %>%
                mutate(compartment = factor(compartment, levels = c("non_unique_compartment",
                                                                    "plus", "minus"),
                                            ordered = TRUE)) %>%
                separate(walk_id, into = c("walk_id", "iteration_clone"), sep = " ") %>%
                filter(walk_id %in% walks_two_or_more_distinct_compartments$pacbio_frag_ID) %>%
                filter(steps <= 6) %>%
                group_by(iteration, state, compartment, steps) %>%
                summarise(prop_same_comp = 1 - sum(different_compartment_index_ints_number)/sum(steps)) %>%
                mutate(set = "Permutated walk")) +
  theme_classic() +
  facet_wrap(~ compartment, labeller = labeller(compartment = supp.labs)) +
  ggtitle("Walks that touch two or more distinct compartments") +
  scale_x_continuous(breaks = 3:10, labels = 3:10) +
  scale_y_continuous("", limits = c(0.2, 0.7)) +
  scale_alpha_discrete("Set of walks")

print(p)

fig3e <- p


# Fig 3F

# Select walks that visit only two distinct compartment domains
walks_just_two_distinct_compartments <- steps_per_frag %>%
  ungroup() %>%
  filter(state == "ns",
         steps <= 10,
         steps >= 3,
         class == 1,
         number_of_compartments == 2,
         different_compartment_index_ints_number > 0)

# Plot fraction of steps within the same compartment domain
# For real walks and permutations
p <- ggplot(mapping = aes(x = steps, y = prop_same_comp, alpha = set)) +
  geom_line(data = walks_just_two_distinct_compartments %>%
              filter(steps <= 6) %>%
              group_by(state, steps) %>%
              summarise(prop_same_comp = 1 - sum(different_compartment_index_ints_number)/sum(steps)) %>%
              mutate(set = "Real walk")) +
  geom_jitter(data = permutated_steps_per_frag %>%
                ungroup() %>%
                mutate(compartment = factor(compartment, levels = c("non_unique_compartment",
                                                                    "plus", "minus"),
                                            ordered = TRUE)) %>%
                separate(walk_id, into = c("walk_id", "iteration_clone"), sep = " ") %>%
                filter(walk_id %in% walks_just_two_distinct_compartments$pacbio_frag_ID) %>%
                filter(steps <= 6) %>%
                group_by(iteration, state, steps) %>%
                summarise(prop_same_comp = 1 - sum(different_compartment_index_ints_number)/sum(steps)) %>%
                mutate(set = "Permutated walk")) +
  theme_classic() +
  ggtitle("") +
  facet_wrap(~ "All such walks") +
  scale_x_continuous(breaks = 3:10, labels = 3:10) +
  scale_y_continuous("Fraction of intra-domain steps", limits = c(0.2, 0.7)) +
  scale_alpha_discrete("Set of walks")

print(p)

fig3fa <- p

# Same as above, but separate A-B, A-A and A-B walks
p <- ggplot(mapping = aes(x = steps, y = prop_same_comp, alpha = set)) +
  geom_line(data = walks_just_two_distinct_compartments %>%
              filter(steps <= 6) %>%
              group_by(state, compartment, steps) %>%
              summarise(prop_same_comp = 1 - sum(different_compartment_index_ints_number)/sum(steps)) %>%
              mutate(set = "Real walk")) +
  geom_jitter(data = permutated_steps_per_frag %>%
                ungroup() %>%
                mutate(compartment = factor(compartment, levels = c("non_unique_compartment",
                                                                    "plus", "minus"),
                                            ordered = TRUE)) %>%
                separate(walk_id, into = c("walk_id", "iteration_clone"), sep = " ") %>%
                filter(walk_id %in% walks_just_two_distinct_compartments$pacbio_frag_ID) %>%
                filter(steps <= 6) %>%
                group_by(iteration, state, compartment, steps) %>%
                summarise(prop_same_comp = 1 - sum(different_compartment_index_ints_number)/sum(steps)) %>%
                mutate(set = "Permutated walk")) +
  theme_classic() +
  facet_wrap(~ compartment, labeller = labeller(compartment = supp.labs)) +
  ggtitle("Walks that touch only two distinct compartments") +
  scale_x_continuous(breaks = 3:10, labels = 3:10) +
  scale_y_continuous("", limits = c(0.2, 0.7)) +
  scale_alpha_discrete("Set of walks")

print(p)

fig3f <- p





cowplot::plot_grid(
cowplot::plot_grid(
  cowplot::plot_grid(
    fig3a + theme(legend.position = "none")
    , ggplot()
    , fig3c + theme(legend.position = "none")
    , ggplot()
    , cowplot::plot_grid(fig3ea + theme(legend.position = "none"),
                         fig3e + theme(legend.position = "none"),
                         ncol = 2, align = "hv", axis = "lb",
                         rel_widths = c(0.3, 0.7))
    , ncol = 1
    , align = "h", axis = "l",
    rel_heights = c(0.3, 0.03, 0.3, 0.03, 0.3),
    labels = c("A", "", "C", "", "E"))
  , cowplot::plot_grid(
    get_legend(fig3a)
    , ggplot()
    , get_legend(fig3c)
    , ggplot()
    , get_legend(fig3e)
    , ncol =1,
    rel_heights = c(0.3, 0.03, 0.3, 0.03, 0.3))
  , rel_widths = c(7,3)
),
ggplot(),
cowplot::plot_grid(
  cowplot::plot_grid(
    cowplot::plot_grid(new_fig3b + theme(legend.position = "none"),
                       ggplot(),
                       fig3b + theme(legend.position = "none"),
                       ncol = 3, align = "hv", axis = "b",
                       rel_widths = c(0.3, 0.1, 0.6))
    , ggplot()
    , fig3d + theme(legend.position = "none")
    , ggplot()
    , cowplot::plot_grid(fig3fa + theme(legend.position = "none"),
                         fig3f + theme(legend.position = "none"),
                         ncol = 2, align = "hv", axis = "b",
                         rel_widths = c(0.3, 0.7))
    , ncol = 1
    , align = "h", axis = "l",
    rel_heights = c(0.3, 0.03, 0.3, 0.03, 0.3),
    labels = c("B", "", "D", "", "F"))
  , cowplot::plot_grid(
    get_legend(fig3b)
    , ggplot()
    , get_legend(fig3d)
    , ggplot()
    , get_legend(fig3f)
    , ncol =1,
    rel_heights = c(0.3, 0.03, 0.3, 0.03, 0.3))
  , rel_widths = c(7,3)
),
ncol = 3,
rel_widths = c(4.75, 0.5, 4.75)
)


