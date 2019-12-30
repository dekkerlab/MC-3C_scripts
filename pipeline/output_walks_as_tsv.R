# Script to output walks as tables

# Load the libraries
library(dplyr)
library(tidyr)
library(readr)

# Start with the walks that still have the NAs.
load("frames_for_plotting_stricter_keep_nas.RData")

# Filter out the walks from the "cluster" samples and those with class == 3
direct_ints <- ungroup(direct_ints) %>%
               dplyr::filter(state != "cluster",
                             class != 3)

steps_per_frag <- ungroup(steps_per_frag) %>%
                  dplyr::filter(state != "cluster",
                                class != 3)

# Output walks as tsv files
write_tsv(filter(direct_ints, state == "ns"),
          path = "walks_tsvs/walks_ns_with_nas.tsv")

write_tsv(filter(direct_ints, state == "mitotic"),
          path = "walks_tsvs/walks_mitotic_with_nas.tsv")

# Output walk summary data as tsv
write_tsv(filter(steps_per_frag, state == "ns"),
          path = "walks_tsvs/walk_summary_ns_with_nas.tsv")

write_tsv(filter(steps_per_frag, state == "mitotic"),
          path = "walks_tsvs/walks_summary_mitotic_with_nas.tsv")

# Remove existing objects
rm(list=ls())


# Load walks without NAs
load("frames_for_plotting_stricter.RData")

# Filter out the walks from the "cluster" samples and those with class == 3
direct_ints <- ungroup(direct_ints) %>%
  dplyr::filter(state != "cluster",
                class != 3)

steps_per_frag <- ungroup(steps_per_frag) %>%
  dplyr::filter(state != "cluster",
                class != 3)

# Output walks as tsv files
write_tsv(filter(direct_ints, state == "ns"),
          path = "walks_tsvs/walks_ns.tsv")

write_tsv(filter(direct_ints, state == "mitotic"),
          path = "walks_tsvs/walks_mitotic.tsv")

# Output walk summary data as tsv
write_tsv(filter(steps_per_frag, state == "ns"),
          path = "walks_tsvs/walk_summary_ns.tsv")

write_tsv(filter(steps_per_frag, state == "mitotic"),
          path = "walks_tsvs/walks_summary_mitotic.tsv")

# remove objects
rm(list=ls())

# Load permutated walks
load("frames_for_plotting_permutated.RData")


# Filter out the walks from the "cluster" samples and those with class == 3
permutated_frame <- ungroup(permutated_frame) %>%
  dplyr::filter(state != "cluster",
                class != 3)

permutated_steps_per_frag <- ungroup(permutated_steps_per_frag) %>%
  dplyr::filter(state != "cluster",
                class != 3)

# Output walks as tsv files
write_tsv(filter(permutated_frame, state == "ns"),
          path = "walks_tsvs/permutated_walks_ns.tsv")

write_tsv(filter(permutated_frame, state == "mitotic"),
          path = "walks_tsvs/permutated_walks_mitotic.tsv")

# Output walk summary data as tsv
write_tsv(filter(permutated_steps_per_frag, state == "ns"),
          path = "walks_tsvs/permutated_walk_summary_ns.tsv")

write_tsv(filter(permutated_steps_per_frag, state == "mitotic"),
          path = "walks_tsvs/permutated_walks_summary_mitotic.tsv")

