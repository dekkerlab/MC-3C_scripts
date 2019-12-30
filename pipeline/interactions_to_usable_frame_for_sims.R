library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

seed <- commandArgs(trailingOnly = TRUE)[[1]]
#mitotic_sim_files <- list.files()[grep("mitotic_class1_sim", list.files(), perl = TRUE)]



load(paste0("sims/mitotic_class1_sim_", seed, ".RData"))
  
mitotic_class1_sim$sim_rep <- as.integer(seed) - 17122016 + 1
  



load(paste0("sims/ns_class1_sim_", seed, ".RData"))
     
ns_class1_sim$sim_rep <- as.integer(seed) - 17122016 + 1
  

mitotic_class1_sim$state <- "mitotic"
ns_class1_sim$state <- "ns"

class1_sim <- rbind(mitotic_class1_sim, ns_class1_sim)

class1_sim <- class1_sim %>%
  separate(start_bin,
           c("frag1_bin_id", "frag1_genome", "frag1_chr", "frag1_start", "frag1_end"),
           convert = TRUE) %>%
  separate(end_bin,
           c("frag2_bin_id", "frag2_genome", "frag2_chr", "frag2_start", "frag2_end"),
           convert = TRUE)

class1_sim <- class1_sim %>%
  mutate(dist = frag1_start - frag2_start,
         walk_id = paste(walk_id, state, sim_rep))




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


# Create a GRanges object with the coordinates for the first fragment of each interaction
class1_sim_first_frag_ranges <- with(class1_sim,
                                     GRanges(seqnames = frag1_chr,
                                             ranges = IRanges(start = as.integer(frag1_start),
                                                              end = as.integer(frag1_end)),
                                             strand = "*"))

# Find the overlaps between the first fragment of an interaction and the compartments
first_frag_compartment_overlaps <- as.data.frame(findOverlaps(class1_sim_first_frag_ranges,
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

# Add the compartment type of the first interaction fragment to the class1_sim frame
class1_sim[first_frag_compartment_overlaps$queryHits,
           "first_frag_compartment_type"] <- first_frag_compartment_overlaps$compartment_type

# Add the compartment index of the first interaction fragment to the class1_sim frame
class1_sim[first_frag_compartment_overlaps$queryHits,
           "first_frag_compartment_index"] <- first_frag_compartment_overlaps$compartment_index



# Create a GRanges object with the coordinates of the second fragment of each interaction
class1_sim_second_frag_ranges <- with(class1_sim,
                                      GRanges(seqnames = frag2_chr,
                                              ranges = IRanges(start = as.integer(frag2_start),
                                                               end = as.integer(frag2_end)),
                                              strand = "*"))

# Find the overlaps between the second fragment of an interaction and the compartments
second_frag_compartment_overlaps <- as.data.frame(findOverlaps(class1_sim_second_frag_ranges,
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

# Add the compartment type of the second interaction fragment to the class1_sim frame
class1_sim[second_frag_compartment_overlaps$queryHits,
           "second_frag_compartment_type"] <- second_frag_compartment_overlaps$compartment_type

# Add the compartment index of the second interaction fragment to the class1_sim frame
class1_sim[second_frag_compartment_overlaps$queryHits,
           "second_frag_compartment_index"] <- second_frag_compartment_overlaps$compartment_index



class1_sim_steps_per_frag <- class1_sim %>%
  group_by(sim_rep, state, walk_id) %>%
  summarise(steps = length(frag1_chr),
            max_span = abs(max(as.integer(c(frag1_start, frag2_start)) - min(as.integer(c(frag1_start, frag2_start))))),
            min_coordinate = min(as.integer(c(frag1_start, frag1_end, frag2_start, frag2_end))),
            max_coordinate = max(as.integer(c(frag1_start, frag1_end, frag2_start, frag2_end))),
            max_int_dist = suppressWarnings(max(abs(dist), na.rm=TRUE)),
            sum_dists = sum(abs(dist)),
            chr = unique(frag1_chr),
            number_of_compartment_types = length(unique(c(first_frag_compartment_type,
                                                          second_frag_compartment_type))),
            compartment = ifelse(length(unique(c(first_frag_compartment_type,
                                                 second_frag_compartment_type))) == 1,
                                 unique(c(first_frag_compartment_type,
                                          second_frag_compartment_type)),
                                 "non_unique_compartment"),
            number_of_compartments = length(unique(c(first_frag_compartment_index,
                                                     second_frag_compartment_index))),
            first_compartment = first_frag_compartment_type[[1]])


class1_sim <- merge(class1_sim,
                    class1_sim_steps_per_frag[, c("walk_id", "steps", "max_span", "sum_dists")])



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
sim_class_1_ranges <- with(class1_sim_steps_per_frag,
                           GRanges(seqnames = chr,
                                   ranges = IRanges(start = min_coordinate,
                                                    end = max_coordinate),
                                   strand = "*",
                                   steps = steps,
                                   state = state,
                                   sim_rep = sim_rep))

# Count how many TADs each class 1 pacbio frag overlaps with
sim_class_1_ranges$within_one_TAD <- countOverlaps(sim_class_1_ranges, hela_tads) == 1



class_1_sim_split <- split(class1_sim, class1_sim$walk_id)

sim_spacings_frame <- lapply(class_1_sim_split, function(x){
  
  x <- as.data.frame(x)
  
  frag_id <- x$walk_id[[1]]
  
  state <- x$state[[1]]
  
  steps <- length(x$walk_id)
  
  sim_rep <- x$sim_rep[[1]]
  
  max_span <- x$max_span[[1]]
  
  x <- rbind(x[, c("frag1_chr", "frag1_end")],
             data.frame(frag1_chr = x[dim(x)[1], "frag2_chr"],
                        frag1_end = x[dim(x)[1], "frag2_end"]))
  
  x$frag1_end <- as.integer(x$frag1_end)
  
  x <- unique(x[order(x$frag1_end, x$frag1_chr), c("frag1_chr", "frag1_end")])
  
  data.frame(frag = frag_id, spacings = diff(x$frag1_end),
             state = state, steps = steps,
             max_span = max_span, sim_rep = sim_rep)
  
})

sim_spacings_frame <- do.call(rbind, sim_spacings_frame)

sim_spacings_frame <- sim_spacings_frame[order(sim_spacings_frame$state, sim_spacings_frame$max_span),]

save(sim_spacings_frame, class1_sim,
     class1_sim_steps_per_frag, sim_class_1_ranges,
     file = paste0("sim_frames/sim_frames_", seed, ".RData"))

