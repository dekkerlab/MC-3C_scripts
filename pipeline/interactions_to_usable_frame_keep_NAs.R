
# Script to load and process the interactions, then save the dataframes for later plotting

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(rtracklayer)

# Load the interactions from the different datasets and concatenate them into a data frame
direct_ints <- data.frame()

long_mitotic_07JUL16 <- read.table("../bwa_outputs/07JUL16_Mitotic_HeLa_3CPacbio_long_digested.sorted_stricter.interactions.txt",
                                   sep = "\t", stringsAsFactors = FALSE)

direct_ints <- rbind(direct_ints, data.frame(long_mitotic_07JUL16, size = "long", state = "mitotic", dataset = "07JUL16"))

short_mitotic_07JUL16 <- read.table("../bwa_outputs/07JUL16_Mitotic_HeLa_3CPacbio_short_digested.sorted_stricter.interactions.txt",
                                    sep = "\t", stringsAsFactors = FALSE)

direct_ints <- rbind(direct_ints, data.frame(short_mitotic_07JUL16, size = "short", state = "mitotic", dataset = "07JUL16"))

long_ns_15AUG16 <- read.table("../bwa_outputs/15AUG16_PacBio_HelaCluster1_long_digested.sorted_stricter.interactions.txt",
                              sep = "\t", stringsAsFactors = FALSE)

direct_ints <- rbind(direct_ints, data.frame(long_ns_15AUG16, size = "long", state = "cluster", dataset = "15AUG16CL"))

long_ns_15AUG16_resequence <- read.table("../bwa_outputs/15AUG16_PacBio_HelaCluster1_long_resequence_digested.sorted_stricter.interactions.txt",
                                         sep = "\t", stringsAsFactors = FALSE)

direct_ints <- rbind(direct_ints, data.frame(long_ns_15AUG16_resequence, size = "long", state = "cluster", dataset = "15AUG16CL_resequence"))

short_ns_15AUG16 <- read.table("../bwa_outputs/15AUG16_PacBio_HelaCluster1_short_digested.sorted_stricter.interactions.txt",
                               sep = "\t", stringsAsFactors = FALSE)

direct_ints <- rbind(direct_ints, data.frame(short_ns_15AUG16, size = "short", state = "cluster", dataset = "15AUG16CL"))

short_ns_15AUG16_resequence <- read.table("../bwa_outputs/15AUG16_PacBio_HelaCluster1_short_resequence_digested.sorted_stricter.interactions.txt",
                                          sep = "\t", stringsAsFactors = FALSE)

direct_ints <- rbind(direct_ints, data.frame(short_ns_15AUG16_resequence, size = "short", state = "cluster", dataset = "15AUG16CL_resequence"))

long_mitotic_15AUG16 <- read.table("../bwa_outputs/15AUG16_PacBio_MitoticHela_long_digested.sorted_stricter.interactions.txt",
                                   sep = "\t", stringsAsFactors = FALSE)

direct_ints <- rbind(direct_ints, data.frame(long_mitotic_15AUG16, size = "long", state = "mitotic", dataset = "15AUG16MT"))

short_mitotic_15AUG16 <- read.table("../bwa_outputs/15AUG16_PacBio_MitoticHela_short_digested.sorted_stricter.interactions.txt",
                                    sep = "\t", stringsAsFactors = FALSE)

direct_ints <- rbind(direct_ints, data.frame(short_mitotic_15AUG16, size = "short", state = "mitotic", dataset = "15AUG16MT"))

long_ns_26APR16 <- read.table("../bwa_outputs/26APR16_PacBio_HelaS3_unsyn_long_digested.sorted_stricter.interactions.txt",
                              sep = "\t", stringsAsFactors = FALSE)

direct_ints <- rbind(direct_ints, data.frame(long_ns_26APR16, size = "long", state = "ns", dataset = "26APR16"))

long_ns_26APR16_resequence <- read.table("../bwa_outputs/26APR16_PacBio_HelaS3_unsyn_long_resequence_digested.sorted_stricter.interactions.txt",
                                         sep = "\t", stringsAsFactors = FALSE)

direct_ints <- rbind(direct_ints, data.frame(long_ns_26APR16_resequence, size = "long", state = "ns", dataset = "26APR16_resequence"))

short_ns_26APR16 <- read.table("../bwa_outputs/26APR16_PacBio_HelaS3_unsyn_short_digested.sorted_stricter.interactions.txt",
                               sep = "\t", stringsAsFactors = FALSE)

direct_ints <- rbind(direct_ints, data.frame(short_ns_26APR16, size = "short", state = "ns", dataset = "26APR16"))

short_ns_26APR16_resequence <- read.table("../bwa_outputs/26APR16_PacBio_HelaS3_unsyn_short_resequence_digested.sorted_stricter.interactions.txt",
                                          sep = "\t", stringsAsFactors = FALSE)

direct_ints <- rbind(direct_ints, data.frame(short_ns_26APR16_resequence, size = "short", state = "ns", dataset = "26APR16_resequence"))

short_ns_21SEP <- read.table("../bwa_outputs/21SEP_PacBio_short_digested.sorted_stricter.interactions.txt",
                               sep = "\t", stringsAsFactors = FALSE)

direct_ints <- rbind(direct_ints, data.frame(short_ns_21SEP, size = "short", state = "ns", dataset = "21SEP"))

long_ns_21SEP <- read.table("../bwa_outputs/21SEP_PacBio_long_digested.sorted_stricter.interactions.txt",
                               sep = "\t", stringsAsFactors = FALSE)

direct_ints <- rbind(direct_ints, data.frame(long_ns_21SEP, size = "long", state = "ns", dataset = "21SEP"))



#Set an extra column with unique pacbio frag IDs
direct_ints$pacbio_frag_ID <- with(direct_ints, paste(state, dataset, size, V7, sep="_"))

direct_ints[direct_ints$V1 == "*", "V1"] <- NA
direct_ints[direct_ints$V2 == "*", "V2"] <- NA
direct_ints[direct_ints$V3 == "*", "V3"] <- NA
direct_ints[direct_ints$V4 == "*", "V4"] <- NA
direct_ints[direct_ints$V5 == "*", "V5"] <- NA
direct_ints[direct_ints$V6 == "*", "V6"] <- NA

# Remove any self-interactions
# Explicity keeping the NAs
direct_ints <- subset(direct_ints, !(V1 == V4 & (V2==V5 | V3==V6)) | is.na(V1) | is.na(V4))

# Create filters for cases where the first or last fragment of a walk are NA
direct_ints <- direct_ints %>%
  group_by(pacbio_frag_ID) %>%
  mutate(first_or_NA = is.na(V1[[1]]),
         last_is_NA = is.na(V4[[n()]]),
         position = seq(from = 1, to = n()),
         first = position == 1,
         last = position == n(),
         is_first_and_NA = is.na(V1[[1]]) & first,
         is_last_and_NA = is.na(V4[[n()]]) & last,
         first_last_na = is_first_and_NA | is_last_and_NA)

# Use that filter to remove cases where the first or last fragment of a walk are NA
direct_ints <- filter(direct_ints, first_last_na == FALSE)

# Add a column for interaction distance, taking the end point of the restriction fragments as reference
direct_ints$dist <- with(direct_ints, as.integer(V6) - as.integer(V3))

# Set the distance to NA for trans interactions
direct_ints[which(direct_ints$V1 != direct_ints$V4), "dist"] <- NA

# Add a column to discriminate trans interactions
direct_ints$inter_chr <- with(direct_ints, V1 != V4)

# Add a column for the directionality of the interactions
direct_ints$direction <- with(direct_ints, ifelse(V6 > V3, "downstream", "upstream"))

# Set trans interaction directionality as NA
direct_ints$direction[direct_ints$inter_chr] <- NA

# Add a column for a unique ID for pacbio fragments
direct_ints$pacbio_frag_ID <- with(direct_ints, paste(state, dataset, size, V7, sep="_"))

## Load frame that has the number of unique ranges per walk
#load("ranges_lengths_frame.RData")
#
## Convert the pacbio_frag_ID to character
#ranges_lengths_frame$pacbio_frag_ID <- as.character(ranges_lengths_frame$pacbio_frag_ID)
#
## Merge direct_ints and ranges_lengths_frame
#direct_ints <- merge(direct_ints, ranges_lengths_frame, by="pacbio_frag_ID")
#
## Subset direct_ints by removing those where the number of reduced ranges is not the number of steps plus one
#direct_ints <- subset(direct_ints, unique_reduced_ranges == (walk_length + 1))

# Split direct_ints by the pacbio_frag_ID
direct_ints_split <- split(direct_ints, direct_ints$pacbio_frag_ID)

# Get the fragment information, avoiding repeating frags that are only visited once
fragment_rangeslist <- lapply(direct_ints_split,
                                function(x) {

                                  n_rows <- dim(x)[[1]]
                                  
                                  ranges <- append(with(subset(x, is.na(V1) == FALSE),
                                                        GRanges(seqnames = V1,
                                                                ranges = IRanges(start = as.numeric(V2),
                                                                                 end = as.numeric(V3)),
                                                                pacbio_frag_ID = pacbio_frag_ID)),
                                                   with(subset(x[n_rows,], is.na(V4) == FALSE),
                                                        GRanges(seqnames = V4,
                                                                ranges = IRanges(start = as.numeric(V5),
                                                                                 end = as.numeric(V6)),
                                                                pacbio_frag_ID = pacbio_frag_ID)))
                                  
                                  ranges$steps_with_NA <- n_rows
                                  ranges$number_of_non_unique_ranges <- length(ranges)
                                  
                                  unique(ranges)
                                  
                                })


unique_reduced_ranges <- sapply(fragment_rangeslist, function(x) length(reduce(x)))

steps_with_NA <- sapply(fragment_rangeslist, function(x) x$steps_with_NA[[1]])

number_of_non_unique_ranges <- sapply(fragment_rangeslist, function(x) x$number_of_non_unique_ranges[[1]])

number_of_unique_ranges <- sapply(fragment_rangeslist, function(x) length(x))

ranges_lengths_frame <- data.frame(pacbio_frag_ID = names(unique_reduced_ranges),
                                   unique_reduced_ranges = unique_reduced_ranges,
                                   number_of_non_unique_ranges = number_of_non_unique_ranges,
                                   number_of_unique_ranges = number_of_unique_ranges)

# Merge direct_ints and ranges_lengths_frame
direct_ints <- merge(direct_ints, ranges_lengths_frame, by="pacbio_frag_ID")

# Subset direct_ints by removing those where the number of reduced ranges is not the number of total ranges
direct_ints <- subset(direct_ints, unique_reduced_ranges == number_of_non_unique_ranges)


# Start processing now for compartment analysis

# Initialise variable with chromosome names for later looping
chrs <- paste0("chr", c(1:22, "X"))

# Initialise the variable that will hold the processed ranges
hela_compartments <- GRanges()

# Start the loop per chromosome
for (chr in chrs) {
  
  # Import the table with the compartment eigenvectors
  hela_chr_compartments <- read.table(paste0("../HeLaS3_DMSO_compartments/TB-HiC-Dpn-DMSOasy2-ACAGTG__hg19__genome__C-250000-iced__",
                                             chr, ".zScore.compartments"),
                                      comment.char = "", header = TRUE)
  
  # Get the ranges for the positive eigenvalue bins and reduce them
  hela_chr_positive_bins <- with(subset(hela_chr_compartments,
                                        eigen1 > 0),
                                 reduce(GRanges(seqnames = X.chr,
                                                ranges = IRanges(start = start,
                                                                 end = end),
                                                strand = "*",
                                                eigen1 = eigen1,
                                                compartment = "plus")))
  
  # Add column to reduced object with the sign
  hela_chr_positive_bins$compartment <- "plus"
  
  # Get the ranges for the negative eigenvalue bins and reduce them
  hela_chr_negative_bins <- with(subset(hela_chr_compartments,
                                        eigen1 < 0),
                                 reduce(GRanges(seqnames = X.chr,
                                                ranges = IRanges(start = start,
                                                                 end = end),
                                                strand = "*",
                                                eigen1 = eigen1,
                                                compartment = "minus")))
  
  # Add column to reduced object with the sign
  hela_chr_negative_bins$compartment <- "minus"
  
  
  # Append this chromosome's compartment ranges to the overall compartment object
  hela_compartments <- suppressWarnings(append(hela_compartments, hela_chr_positive_bins))
  hela_compartments <- suppressWarnings(append(hela_compartments, hela_chr_negative_bins))
  
}


# Add a column with the original direct_ints index number
direct_ints$index <- 1:dim(direct_ints)[[1]]

# Create a GRanges object with the coordinates for the first fragment of each interaction
direct_ints_first_frag_ranges <- with(subset(direct_ints, is.na(V1) == FALSE),
                                      GRanges(seqnames = V1,
                                              ranges = IRanges(start = as.integer(V2),
                                                               end = as.integer(V3)),
                                              direct_ints_index = index))

# Find the overlaps between the first fragment of an interaction and the compartments
first_frag_compartment_overlaps <- as.data.frame(findOverlaps(direct_ints_first_frag_ranges,
                                                              hela_compartments))

# Add the compartment type (plus or minus sign on the eigenvector) to the overlap frame
first_frag_compartment_overlaps$compartment_type <- hela_compartments[first_frag_compartment_overlaps$subjectHits]$compartment

# Get only only the fragments that overlap with a single compartment
first_frag_compartment_overlaps <- first_frag_compartment_overlaps %>%
  group_by(queryHits) %>%
  summarise(number_hits = length(compartment_type),
            compartment_type = unique(compartment_type)[[1]],
            compartment_index = unique(subjectHits)[[1]]) %>%
  subset(number_hits == 1)

# Connect the queryHits to the original direct_ints index
first_frag_compartment_overlaps$direct_ints_index <- subset(direct_ints, is.na(V1) == FALSE)[first_frag_compartment_overlaps$queryHits, "index"]

# Add the compartment type of the first interaction fragment to the direct_ints frame
direct_ints[first_frag_compartment_overlaps$direct_ints_index,
            "first_frag_compartment_type"] <- first_frag_compartment_overlaps$compartment_type

# Add the compartment index of the first interaction fragment to the direct_ints frame
direct_ints[first_frag_compartment_overlaps$direct_ints_index,
            "first_frag_compartment_index"] <- first_frag_compartment_overlaps$compartment_index



# Create a GRanges object with the coordinates of the second fragment of each interaction
direct_ints_second_frag_ranges <- with(subset(direct_ints, is.na(V4) == FALSE),
                                       GRanges(seqnames = V4,
                                               ranges = IRanges(start = as.integer(V5),
                                                                end = as.integer(V6)),
                                               direct_ints_index = index))

# Find the overlaps between the second fragment of an interaction and the compartments
second_frag_compartment_overlaps <- as.data.frame(findOverlaps(direct_ints_second_frag_ranges,
                                                               hela_compartments))

# Add the compartment type (plus or minus sign on the eigenvector) to the overlap frame
second_frag_compartment_overlaps$compartment_type <- hela_compartments[second_frag_compartment_overlaps$subjectHits]$compartment

# Get only only the fragments that overlap with a single compartment
second_frag_compartment_overlaps <- second_frag_compartment_overlaps %>%
  group_by(queryHits) %>%
  summarise(number_hits = length(compartment_type),
            compartment_type = unique(compartment_type)[[1]],
            compartment_index = unique(subjectHits)[[1]]) %>%
  subset(number_hits == 1)

# Connect the queryHits to the original direct_ints index
second_frag_compartment_overlaps$direct_ints_index <- subset(direct_ints, is.na(V4) == FALSE)[second_frag_compartment_overlaps$queryHits, "index"]

# Add the compartment type of the second interaction fragment to the direct_ints frame
direct_ints[second_frag_compartment_overlaps$direct_ints_index,
            "second_frag_compartment_type"] <- second_frag_compartment_overlaps$compartment_type

# Add the compartment index of the second interaction fragment to the direct_ints frame
direct_ints[second_frag_compartment_overlaps$direct_ints_index,
            "second_frag_compartment_index"] <- second_frag_compartment_overlaps$compartment_index



# Add information relative to ChromHMM states

chromhmm <- import.bed("../wgEncodeAwgSegmentationChromhmmHelas3.bed.gz")

state_to_simplified <- c("Tss" = "Active", "TssF" = "Active",
                         "PromF" = "Active", "PromP" = "Active",
                         "Enh" = "Active", "EnhF" = "Active",
                         "EnhWF" = "Active", "EnhW" = "Active",
                         "DnaseU" = "Active", "DnaseD" = "Active",
                         "FaireW" = "Active",
                         "CtcfO" = "Active", "Ctcf" = "Active",
                         "Gen5'" = "Active", "Gen3'" = "Active",
                         "Elon" = "Active", "ElonW" = "Active",
                         "Pol2" = "Active", "H4K20" = "Active",
                         "Low" = "Active",
                         "ReprD" = "Repressed", "Repr" = "Repressed", "ReprW" = "Repressed",
                         "Quies" = "Inactive", "Art" = "Inactive")

chromhmm$simplified <- state_to_simplified[chromhmm$name]



# Create a GRanges object with the coordinates for the first fragment of each interaction
direct_ints_first_frag_ranges <- with(subset(direct_ints, is.na(V1) == FALSE),
                                      GRanges(seqnames = V1,
                                              ranges = IRanges(start = as.integer(V2),
                                                               end = as.integer(V3)),
                                              direct_ints_index = index))

# Find the overlaps between the first fragment of an interaction and the compartments
first_frag_state_overlaps <- as.data.frame(findOverlaps(direct_ints_first_frag_ranges,
                                                        chromhmm))

first_frag_state_overlaps$simplified <- chromhmm[first_frag_state_overlaps$subjectHits]$simplified

first_frag_state_overlaps <- first_frag_state_overlaps %>%
  select(queryHits, simplified) %>%
  unique() %>% group_by(queryHits) %>%
  mutate(n = n()) %>% filter(n == 1)

direct_ints[first_frag_state_overlaps$queryHits, "first_frag_simplified_state"] <- first_frag_state_overlaps$simplified



direct_ints_second_frag_ranges <- with(subset(direct_ints, is.na(V4) == FALSE),
                                       GRanges(seqnames = V4,
                                               ranges = IRanges(start = as.integer(V5),
                                                                end = as.integer(V6)),
                                               direct_ints_index = index))

# Find the overlaps between the second fragment of an interaction and the compartments
second_frag_state_overlaps <- as.data.frame(findOverlaps(direct_ints_second_frag_ranges,
                                                         chromhmm))

second_frag_state_overlaps$simplified <- chromhmm[second_frag_state_overlaps$subjectHits]$simplified

second_frag_state_overlaps <- second_frag_state_overlaps %>%
  select(queryHits, simplified) %>%
  unique() %>% group_by(queryHits) %>%
  mutate(n = n()) %>% filter(n == 1)

direct_ints[second_frag_state_overlaps$queryHits, "second_frag_simplified_state"] <- second_frag_state_overlaps$simplified


# Create a new dataframe that summarises some information per pacbio frag
steps_per_frag <- direct_ints %>%
  group_by(size, state, dataset, V7, pacbio_frag_ID) %>%
  summarise(steps = length(V1),
            inter_chr = sum(inter_chr, na.rm = TRUE),
            inter_char_proportion = sum(inter_chr, na.rm = TRUE)  / length(V1),
            num_chr = length(unique(c(V1, V4)[!is.na(c(V1, V4))])),
            max_span = ifelse(num_chr == 1, abs(max(as.integer(c(V3, V6)), na.rm = TRUE) - min(as.integer(c(V3, V6)), na.rm = TRUE)), NA),
            min_coordinate = min(as.integer(c(V2, V3, V5, V6)), na.rm = TRUE),
            max_coordinate = max(as.integer(c(V2, V3, V5, V6)), na.rm = TRUE),
            max_int_dist = suppressWarnings(max(abs(dist), na.rm=TRUE)),
            min_int_dist = suppressWarnings(min(abs(dist), na.rm=TRUE)),
            sum_dists = sum(abs(dist), na.rm = TRUE),
            first_to_last_frag_dist = ifelse(num_chr == 1, abs(as.integer(last(V6)) - as.integer(first(V3))), NA),
            chr = ifelse( length(unique(c(V1, V4)[!is.na(c(V1, V4))]))==1, as.character(V1[!is.na(V1)]), "*"),
            number_of_compartment_types = length(unique(c(first_frag_compartment_type,
                                                          second_frag_compartment_type)[!is.na(c(first_frag_compartment_type,
                                                                                                 second_frag_compartment_type))])),
            compartment = ifelse(length(unique(c(first_frag_compartment_type,
                                                 second_frag_compartment_type)[!is.na(c(first_frag_compartment_type,
                                                                                        second_frag_compartment_type))])) == 1,
                                 unique(c(first_frag_compartment_type,
                                          second_frag_compartment_type)[!is.na(c(first_frag_compartment_type,
                                                                                 second_frag_compartment_type))]),
                                 "non_unique_compartment"),
            number_of_compartments = length(unique(c(first_frag_compartment_index,
                                                     second_frag_compartment_index)[!is.na(c(first_frag_compartment_type,
                                                                                             second_frag_compartment_type))])),
            first_compartment = first_frag_compartment_type[[1]],
            percentage_downstream = sum(direction == "downstream") / length(direction),
            number_of_simplified_chrom_states = length(unique(c(first_frag_simplified_state, second_frag_simplified_state)[!is.na(c(first_frag_simplified_state, second_frag_simplified_state))])),
            simplified_chrom_state = ifelse(length(unique(c(first_frag_simplified_state, second_frag_simplified_state)[!is.na(c(first_frag_simplified_state, second_frag_simplified_state))])) == 1,
                                            unique(c(first_frag_simplified_state, second_frag_simplified_state, na.rm = TRUE)),
                                            "non_unique_chrom_state"))

# Add a class column to the steps_per_frag data frame
# Class 1 has frags that are only in one chromosome
# Class 2 has frags that are only in two chromosomes
# Class 3 has frags in three or more chromosomes
steps_per_frag$class <- steps_per_frag$num_chr
steps_per_frag$class[steps_per_frag$class > 3] <- 3

# Set max span to NA if equals 0
steps_per_frag <- mutate(steps_per_frag, max_span = ifelse(max_span == 0, NA, max_span))

# Order the factors in the compartment column
steps_per_frag$compartment <- factor(steps_per_frag$compartment,
                                     levels = c("non_unique_compartment", "plus", "minus"),
                                     ordered=TRUE)

# Summarise the steps_per_frag data frame to get a frame with information per number of steps in a pacbio frag
info_per_step <- steps_per_frag %>%
  subset(class != 3) %>%
  group_by(size, state, dataset, steps) %>%
  summarize(intra_chromosome = sum(num_chr == 1)/length(num_chr),
            fraction_cis_hops = 1 - (sum(inter_chr)/sum(steps)),
            fraction_one_compartment_type = sum(number_of_compartment_types == 1) / length(number_of_compartment_types))

# Summarise the steps_per_frag data frame to get a frame with information per number of steps in a pacbio frag
# Do not separate datasets from the same condition
info_per_step_by_state <- steps_per_frag %>%
  subset(class != 3) %>%
  group_by(state, steps) %>%
  summarize(intra_chromosome = sum(num_chr == 1)/length(num_chr),
            fraction_cis_hops = 1 - (sum(inter_chr)/sum(steps)),
            fraction_one_compartment_type = sum(number_of_compartment_types == 1) / length(number_of_compartment_types))


# Merge the information of the info_per_step data frame to the direct_ints frame
direct_ints <- merge(direct_ints, steps_per_frag[, c("pacbio_frag_ID", "class", "steps", "max_span", "sum_dists")])

# Add columns to the direct_ints frame to have the ratios of certain distances
direct_ints$ratio_80kb <- with(direct_ints, abs(dist) / 80e3)
direct_ints$ratio_100kb <- with(direct_ints, abs(dist) / 100e3)
direct_ints$ratio_120kb <- with(direct_ints, abs(dist) / 120e3)


# Summarise the steps_per_frag data frame to get a frame with information per number of steps and per class in a pacbio frag
info_per_step_per_class <- steps_per_frag %>%
  subset(class != 3) %>%
  group_by(size, state, dataset, steps, class) %>%
  summarize(intra_chromosome = sum(num_chr == 1)/length(num_chr),
            fraction_cis_hops = 1 - (sum(inter_chr)/sum(steps)),
            fraction_one_compartment = sum(number_of_compartment_types == 1) / length(number_of_compartment_types))


# Get the fraction of fragments restricted to a certain size
size_restriction <- steps_per_frag %>%
  subset(class == 1) %>%
  mutate(within_1Mb = (max_span <= 1e6),
         within_2Mb = (max_span <= 2e6),
         within_3Mb = (max_span <= 3e6)) %>%
  group_by(size, state, dataset, steps) %>%
  summarise(fraction_within_1Mb = sum(within_1Mb) / length(within_1Mb),
            fraction_within_2Mb = sum(within_2Mb) / length(within_2Mb),
            fraction_within_3Mb = sum(within_3Mb) / length(within_3Mb)) %>%
  gather(size_selection, fraction, fraction_within_1Mb:fraction_within_3Mb)

size_restriction_by_state <- steps_per_frag %>%
  subset(class == 1) %>%
  mutate(within_1Mb = (max_span <= 1e6),
         within_2Mb = (max_span <= 2e6),
         within_3Mb = (max_span <= 3e6)) %>%
  group_by(state, steps) %>%
  summarise(fraction_within_1Mb = sum(within_1Mb) / length(within_1Mb),
            fraction_within_2Mb = sum(within_2Mb) / length(within_2Mb),
            fraction_within_3Mb = sum(within_3Mb) / length(within_3Mb)) %>%
  gather(size_selection, fraction, fraction_within_1Mb:fraction_within_3Mb)


# Add a column with a unique dataset ID to each data frame to help with plot colouring and faceting
direct_ints$full_dataset_ID <- with(direct_ints, paste(state, size, dataset))
steps_per_frag$full_dataset_ID <- with(steps_per_frag, paste(state, size, dataset))
info_per_step$full_dataset_ID <- with(info_per_step, paste(state, size, dataset))
info_per_step_per_class$full_dataset_ID <- with(info_per_step_per_class, paste(state, size, dataset))
size_restriction$full_dataset_ID <- with(size_restriction, paste(state, size, dataset))

# Start processing now for TAD analysis

# Initialise variable with chromosome names for later looping
chrs <- paste0("chr", c(1:22, "X"))

# Initialise the variable that will hold the processed ranges
hela_tads <- GRanges()

# Start the loop per chromosome
for (chr in chrs) {
  
  # Import the bedgraph with the insulation values
  hela_chr_insulation <- import(paste0("../HeLaS3_DMSO_insulation/TB-HiC-Hela-DMSOasy2__hg19__genome__C-40000-iced__",
                                       chr, "--is520000--nt0--ids320000--ss0--immean.insulation.bedGraph"))
  
  # Import the bed file with the TAD borders
  hela_chr_border <- import(paste0("../HeLaS3_DMSO_insulation/TB-HiC-Hela-DMSOasy2__hg19__genome__C-40000-iced__",
                                   chr, "--is520000--nt0--ids320000--ss0--immean.insulation.boundaries.bed"))
  
  # Get the bins that are inside a TAD, ie. that are not borders
  hela_chr_tads <- hela_chr_insulation[which(countOverlaps(hela_chr_insulation, hela_chr_border)==FALSE)]
  
  # Slightly resize the bins so that contiguous bins have a one bp overlap and then reduce them to merge them into a single TAD range
  hela_chr_tads <- reduce(resize(hela_chr_tads, 40001))
  
  # Append this chromosome's TAD ranges to the overall TADs object
  hela_tads <- suppressWarnings(append(hela_tads, hela_chr_tads))
  
}

# Convert class 1 fragments to GRanges
class_1_ranges <- with(subset(steps_per_frag, class == 1 & !is.na(max_span)),
                       GRanges(seqnames = chr,
                               ranges = IRanges(start = min_coordinate,
                                                end = max_coordinate),
                               strand = "*",
                               steps = steps,
                               size = size,
                               state = state,
                               dataset = dataset,
                               full_dataset_ID = full_dataset_ID,
                               pacbio_frag_ID = pacbio_frag_ID))

# Count how many TADs each class 1 pacbio frag overlaps with
class_1_ranges$within_one_TAD <- countOverlaps(class_1_ranges, hela_tads) == 1

# Create a GRanges object with the coordinates for the first fragment of each interaction
direct_ints_first_frag_ranges <- with(subset(direct_ints, is.na(V1) == FALSE),
                                      GRanges(seqnames = V1,
                                              ranges = IRanges(start = as.integer(V2),
                                                               end = as.integer(V3)),
                                              direct_ints_index = index))

# Find the overlaps between the first fragment of an interaction and the TADs
first_frag_tad_overlaps <- as.data.frame(findOverlaps(direct_ints_first_frag_ranges,
                                                        hela_tads))

# Get only only the fragments that overlap with a single tad
first_frag_tad_overlaps <- first_frag_tad_overlaps %>%
  group_by(queryHits) %>%
  summarise(number_hits = n(),
            subjectHits = unique(subjectHits)[[1]]) %>%
  subset(number_hits == 1)

# Connect the queryHits to the original direct_ints index
first_frag_tad_overlaps$direct_ints_index <- subset(direct_ints, is.na(V1) == FALSE)[first_frag_tad_overlaps$queryHits, "index"]

# Add the tad index of the first interaction fragment to the direct_ints frame
direct_ints[first_frag_tad_overlaps$direct_ints_index,
            "first_frag_tad_index"] <- first_frag_tad_overlaps$subjectHits

# Create a GRanges object with the coordinates of the second fragment of each interaction
direct_ints_second_frag_ranges <- with(subset(direct_ints, is.na(V4) == FALSE),
                                       GRanges(seqnames = V4,
                                               ranges = IRanges(start = as.integer(V5),
                                                                end = as.integer(V6)),
                                               direct_ints_index = index))

# Find the overlaps between the second fragment of an interaction and the TADs
second_frag_tad_overlaps <- as.data.frame(findOverlaps(direct_ints_second_frag_ranges,
                                                      hela_tads))

# Get only only the fragments that overlap with a single tad
second_frag_tad_overlaps <- second_frag_tad_overlaps %>%
  group_by(queryHits) %>%
  summarise(number_hits = n(),
            subjectHits = unique(subjectHits)[[1]]) %>%
  subset(number_hits == 1)

# Connect the queryHits to the original direct_ints index
second_frag_tad_overlaps$direct_ints_index <- subset(direct_ints, is.na(V4) == FALSE)[second_frag_tad_overlaps$queryHits, "index"]

# Add the tad index of the second interaction fragment to the direct_ints frame
direct_ints[second_frag_tad_overlaps$direct_ints_index,
            "second_frag_tad_index"] <- second_frag_tad_overlaps$subjectHits


# Let me try something with the minimum distance between consecutive cis fragments
# Split the class 1 interactions so we can do a lapply on it
class_1_ints_split <- split(subset(direct_ints, class == "1"), subset(direct_ints, class == "1")$pacbio_frag_ID)

spacings_frame <- lapply(class_1_ints_split, function(x){
  
  frag_id <- x$pacbio_frag_ID[[1]]
  
  state <- x$state[[1]]
  
  steps <- x$steps[[1]]
  
  max_span <- x$max_span[[1]]
  
  full_dataset_ID <- x$full_dataset_ID[[1]]
  
  x <- rbind(x[, c("V1", "V3")], data.frame(V1 = x[dim(x)[1], "V4"], V3 = x[dim(x)[1], "V6"]))
  
  x$V3 <- as.integer(x$V3)
  
  x <- unique(x[order(x$V3, x$V1), c("V1", "V3")])
  
  data.frame(frag = frag_id, spacings = diff(x$V3),
             state = state, full_dataset_ID = full_dataset_ID,
             steps = steps, max_span = max_span)
  
})

spacings_frame <- do.call(rbind, spacings_frame)

spacings_frame <- spacings_frame[order(spacings_frame$full_dataset_ID, spacings_frame$max_span),]

spacings_frame$frag_factor <- factor(spacings_frame$frag, levels = unique(spacings_frame$frag), ordered = TRUE)

spacings_frame <- spacings_frame %>%
  group_by(state, full_dataset_ID, steps, frag) %>%
  do(data.frame(., ranks = rank(-.$spacings)))

# Get the max_spans for each group of cis fragments in a walk
direct_ints_split <- split(direct_ints, direct_ints$pacbio_frag_ID)

max_span_by_chr <- lapply(direct_ints_split, function(x) {
  last_frag <- x[dim(x)[1], c("V4", "V5", "V6")]
  colnames(last_frag) <- c("V1", "V2", "V3")
  fragments <- rbind(x[, c("V1", "V2", "V3")], last_frag)
  
  fragments$index <- 1:dim(fragments)[1]
  
  fragment_ranges <- with(subset(fragments, is.na(V1) == FALSE),
                          GRanges(seqnames = V1,
                                  ranges = IRanges(start = as.integer(V2),
                                                   end = as.integer(V3)),
                                  index = index))
  
  fragments_compartment_overlaps <- as.data.frame(findOverlaps(fragment_ranges,
                                                               hela_compartments))
  
  fragments_compartment_overlaps$index <- subset(fragments, is.na(V1) == FALSE)[fragments_compartment_overlaps$queryHits, "index"]
  

  fragments[fragments_compartment_overlaps$index, "compartment"] <- hela_compartments[fragments_compartment_overlaps$subjectHits]$compartment

  fragments %>%
    filter(is.na(V1) == FALSE) %>%
    group_by(V1) %>%
    summarise(max_span_per_chr = max(as.integer(V3), na.rm = TRUE) - min(as.integer(V3), na.rm = TRUE),
              class = unique(x$class),
              pacbio_frag_ID = unique(x$pacbio_frag_ID),
              state = unique(x$state),
              steps = unique(x$steps),
              number_compartments = length(compartment),
              number_a_compartments = sum(compartment == "plus"),
              number_b_compartments = sum(compartment == "minus"),
              compartment = ifelse(length(unique(compartment)) == 1,
                                   unique(compartment),
                                   "non_unique_compartment"))
    
})

max_span_by_chr <- do.call(rbind, max_span_by_chr)

# Create frames with the data for the full pairwise (indirect) interactions between fragments of the same walk
pairwise_cis_distances <- lapply(direct_ints_split, function(x) {
  
  non_direct_ints <- combn(1:(length(x["V1"][[1]])+1), 2, function(x) x[[1]] != (x[[2]] - 1))
  
  pairwise_distances <- combn(as.integer(c(x["V3"][[1]], x["V6"][dim(x)[1], "V6"])), 2, function(x) abs(diff(x)))
  pairwise_intrachr <- combn(c(x["V1"][[1]], x["V4"][dim(x)[1], "V4"]), 2, function(x) x[[1]] == x[[2]])
  
  return(data.frame(dist = pairwise_distances[which(pairwise_intrachr & non_direct_ints)],
                    state = rep(x["state"][[1]][[1]], length(which(pairwise_intrachr & non_direct_ints))),
                    steps = rep(x["steps"][[1]][[1]], length(which(pairwise_intrachr & non_direct_ints)))))
  
})

pairwise_cis_distances <- do.call(rbind, pairwise_cis_distances)

pairwise_intrachr <- lapply(direct_ints_split, function(x) {

  non_direct_ints <- combn(1:(length(x["V1"][[1]])+1), 2, function(x) x[[1]] != (x[[2]] - 1))
  
  pairwise_intrachr <- combn(c(x["V1"][[1]], x["V4"][dim(x)[1], "V4"]), 2, function(x) x[[1]] == x[[2]])
  
  return(data.frame(intrachr = pairwise_intrachr[non_direct_ints],
                    state = rep(x["state"][[1]][[1]], length(pairwise_intrachr[non_direct_ints])),
                    steps = rep(x["steps"][[1]][[1]], length(pairwise_intrachr[non_direct_ints])),
                    pacbio_frag_ID = rep(x["pacbio_frag_ID"][[1]][[1]], length(pairwise_intrachr[non_direct_ints]))))
  
})

pairwise_intrachr <- do.call(rbind, pairwise_intrachr)



# Get the max_spans for each group of cis fragments in a walk
direct_ints_split <- split(direct_ints, direct_ints$pacbio_frag_ID)

max_span_by_compartment <- lapply(direct_ints_split, function(x) {
  last_frag <- x[dim(x)[1], c("V4", "V5", "V6", "second_frag_compartment_index", "second_frag_compartment_type")]
  colnames(last_frag) <- c("V1", "V2", "V3", "first_frag_compartment_index", "first_frag_compartment_type")
  fragments <- rbind(x[, c("V1", "V2", "V3", "first_frag_compartment_index", "first_frag_compartment_type")], last_frag)
  
  fragments %>%
    filter(is.na(V1) == FALSE) %>%
    group_by(first_frag_compartment_index) %>%
    summarise(max_span_per_compartment = max(as.integer(V3), na.rm = TRUE) - min(as.integer(V3), na.rm = TRUE),
              class = unique(x$class),
              pacbio_frag_ID = unique(x$pacbio_frag_ID),
              state = unique(x$state),
              steps = unique(x$steps),
              compartment_type = unique(first_frag_compartment_type))
  
})

max_span_by_compartment <- do.call(rbind, max_span_by_compartment)

# Save the objects for later use
save(direct_ints, steps_per_frag,
     info_per_step, info_per_step_per_class,
     info_per_step_by_state,
     size_restriction, class_1_ranges,
     size_restriction_by_state,
     hela_tads, hela_compartments,
     spacings_frame,
     max_span_by_chr,
     pairwise_cis_distances,
     pairwise_intrachr,
     max_span_by_compartment,
     file = "frames_for_plotting_stricter_keep_nas.RData")
