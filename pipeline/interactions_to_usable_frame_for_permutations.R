
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(scales)


seed <- commandArgs(trailingOnly = TRUE)[[1]]
#mitotic_permutated_files <- list.files()[grep("mitotic_permutated_direct_ints", list.files(), perl = TRUE)]

load(paste0("../permutated_direct_ints/permutated_direct_ints", seed, ".RData"))

colnames(permutated_direct_ints)[1:6] <- paste0("V", 1:6)
colnames(permutated_direct_ints)[11] <- "iteration"


permutated_direct_ints <- permutated_direct_ints %>%
  mutate(dist = ifelse(V1 == V4, V2 - V5, NA),
         walk_id = paste(pacbio_frag_ID, iteration))

permutated_direct_ints <- permutated_direct_ints %>% separate(full_dataset_ID, c("state", "size", "dataset"), remove = FALSE)


permutated_direct_ints$inter_chr <- with(permutated_direct_ints, V1 != V4)


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
permutated_direct_ints_first_frag_ranges <- with(permutated_direct_ints,
                                     GRanges(seqnames = V1,
                                             ranges = IRanges(start = as.integer(V2),
                                                              end = as.integer(V3)),
                                             strand = "*"))

# Find the overlaps between the first fragment of an interaction and the compartments
first_frag_compartment_overlaps <- as.data.frame(findOverlaps(permutated_direct_ints_first_frag_ranges,
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

# Add the compartment type of the first interaction fragment to the permutated_direct_ints frame
permutated_direct_ints[first_frag_compartment_overlaps$queryHits,
           "first_frag_compartment_type"] <- first_frag_compartment_overlaps$compartment_type

# Add the compartment index of the first interaction fragment to the permutated_direct_ints frame
permutated_direct_ints[first_frag_compartment_overlaps$queryHits,
           "first_frag_compartment_index"] <- first_frag_compartment_overlaps$compartment_index



# Create a GRanges object with the coordinates of the second fragment of each interaction
permutated_direct_ints_second_frag_ranges <- with(permutated_direct_ints,
                                      GRanges(seqnames = V4,
                                              ranges = IRanges(start = as.integer(V5),
                                                               end = as.integer(V6)),
                                              strand = "*"))

# Find the overlaps between the second fragment of an interaction and the compartments
second_frag_compartment_overlaps <- as.data.frame(findOverlaps(permutated_direct_ints_second_frag_ranges,
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

# Add the compartment type of the second interaction fragment to the permutated_direct_ints frame
permutated_direct_ints[second_frag_compartment_overlaps$queryHits,
           "second_frag_compartment_type"] <- second_frag_compartment_overlaps$compartment_type

# Add the compartment index of the second interaction fragment to the permutated_direct_ints frame
permutated_direct_ints[second_frag_compartment_overlaps$queryHits,
           "second_frag_compartment_index"] <- second_frag_compartment_overlaps$compartment_index



permutated_direct_ints_steps_per_frag <- permutated_direct_ints %>%
  group_by(iteration, state, walk_id) %>%
  summarise(steps = length(V1),
            inter_chr = sum(inter_chr, na.rm = TRUE),
            inter_char_proportion = sum(inter_chr, na.rm = TRUE)  / length(V1),
            max_span = abs(max(as.integer(c(V2, V5)) - min(as.integer(c(V2, V5))))),
            min_coordinate = min(as.integer(c(V2, V3, V5, V6))),
            max_coordinate = max(as.integer(c(V2, V3, V5, V6))),
            max_int_dist = suppressWarnings(max(abs(dist), na.rm=TRUE)),
            num_chr = length(unique(c(V1, V4))),
            sum_dists = sum(abs(dist)),
            chr = ifelse( length(unique(c(V1, V4)))==1, as.character(V1), "*"),
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

# Add a class column to the steps_per_frag data frame
# Class 1 has frags that are only in one chromosome
# Class 2 has frags that are only in two chromosomes
# Class 3 has frags in three or more chromosomes
permutated_direct_ints_steps_per_frag$class <- permutated_direct_ints_steps_per_frag$num_chr
permutated_direct_ints_steps_per_frag$class[permutated_direct_ints_steps_per_frag$class > 3] <- 3

permutated_direct_ints <- merge(permutated_direct_ints,
                    permutated_direct_ints_steps_per_frag[, c("walk_id", "steps", "max_span", "sum_dists", "class")])



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
permutated_class_1_ranges <- with(permutated_direct_ints_steps_per_frag,
                           GRanges(seqnames = chr,
                                   ranges = IRanges(start = min_coordinate,
                                                    end = max_coordinate),
                                   strand = "*",
                                   steps = steps,
                                   state = state,
                                   iteration = iteration))

# Count how many TADs each class 1 pacbio frag overlaps with
permutated_class_1_ranges$within_one_TAD <- countOverlaps(permutated_class_1_ranges, hela_tads) == 1

# Create a GRanges object with the coordinates for the first fragment of each interaction
direct_ints_first_frag_ranges <- with(subset(permutated_direct_ints, is.na(V1) == FALSE),
                                      GRanges(seqnames = V1,
                                              ranges = IRanges(start = as.integer(V2),
                                                               end = as.integer(V3))))

# Find the overlaps between the first fragment of an interaction and the TADs
first_frag_tad_overlaps <- as.data.frame(findOverlaps(direct_ints_first_frag_ranges,
                                                      hela_tads))

# Get only only the fragments that overlap with a single tad
first_frag_tad_overlaps <- first_frag_tad_overlaps %>%
  group_by(queryHits) %>%
  summarise(number_hits = n(),
            subjectHits = unique(subjectHits)[[1]]) %>%
  subset(number_hits == 1)


# Add the tad index of the first interaction fragment to the direct_ints frame
permutated_direct_ints[first_frag_tad_overlaps$queryHits,
            "first_frag_tad_index"] <- first_frag_tad_overlaps$subjectHits

# Create a GRanges object with the coordinates of the second fragment of each interaction
direct_ints_second_frag_ranges <- with(subset(permutated_direct_ints, is.na(V4) == FALSE),
                                       GRanges(seqnames = V4,
                                               ranges = IRanges(start = as.integer(V5),
                                                                end = as.integer(V6))))

# Find the overlaps between the second fragment of an interaction and the TADs
second_frag_tad_overlaps <- as.data.frame(findOverlaps(direct_ints_second_frag_ranges,
                                                       hela_tads))

# Get only only the fragments that overlap with a single tad
second_frag_tad_overlaps <- second_frag_tad_overlaps %>%
  group_by(queryHits) %>%
  summarise(number_hits = n(),
            subjectHits = unique(subjectHits)[[1]]) %>%
  subset(number_hits == 1)


# Add the tad index of the second interaction fragment to the direct_ints frame
permutated_direct_ints[second_frag_tad_overlaps$queryHits,
            "second_frag_tad_index"] <- second_frag_tad_overlaps$subjectHits



class_1_permutated_split <- split(permutated_direct_ints, permutated_direct_ints$walk_id)

permutated_spacings_frame <- lapply(class_1_permutated_split, function(x){
  
  x <- as.data.frame(x)
  
  frag_id <- x$walk_id[[1]]
  
  state <- x$state[[1]]
  
  steps <- length(x$walk_id)
  
  iteration <- x$iteration[[1]]
  
  max_span <- x$max_span[[1]]
  
  x <- rbind(x[, c("V1", "V3")],
             data.frame(V1 = x[dim(x)[1], "V4"],
                        V3 = x[dim(x)[1], "V6"]))
  
  x$V3 <- as.integer(x$V3)
  
  x <- unique(x[order(x$V3, x$V1), c("V1", "V3")])
  
  data.frame(frag = frag_id, spacings = diff(x$V3),
             state = state, steps = steps,
             max_span = max_span, iteration = iteration)
  
})

permutated_spacings_frame <- do.call(rbind, permutated_spacings_frame)

permutated_spacings_frame <- permutated_spacings_frame[order(permutated_spacings_frame$state, permutated_spacings_frame$max_span),]

# Get the max_spans for each group of cis fragments in a walk
permutated_direct_ints_split <- split(permutated_direct_ints, permutated_direct_ints$walk_id)

permutated_direct_ints_max_span_by_chr <- lapply(permutated_direct_ints_split, function(x) {
  last_frag <- x[dim(x)[1], c("V4", "V5", "V6")]
  colnames(last_frag) <- c("V1", "V2", "V3")
  fragments <- rbind(x[, c("V1", "V2", "V3")], last_frag)
  
  fragments %>%
    filter(is.na(V1) == FALSE) %>%
    group_by(V1) %>%
    summarise(max_span_per_chr = max(as.integer(V3), na.rm = TRUE) - min(as.integer(V3), na.rm = TRUE),
              class = unique(x$class),
              walk_id = unique(x$walk_id),
              state = unique(x$state),
              steps = unique(x$steps),
	      iteration = unique(x$iteration))
    
})

permutated_direct_ints_max_span_by_chr <- do.call(rbind, permutated_direct_ints_max_span_by_chr)


save(permutated_spacings_frame, permutated_direct_ints,
     permutated_direct_ints_steps_per_frag, permutated_class_1_ranges,
     permutated_direct_ints_max_span_by_chr,
     file = paste0("permutated_frames/permutated_frames_", seed, ".RData"))


