#!/usr/bin/env Rscript
# Calls the pinned upstream SeedMatchR sequence functions. No substitute matcher.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop("Expected repo, candidates, reference, sites, provenance, and seed classes")
if (getRversion() < "4.3.0") stop("SeedMatchR requires R >= 4.3.0")
if (!requireNamespace("Biostrings", quietly = TRUE)) stop("Biostrings is missing; run scripts/install_r.sh")

# The iranges branch requires Biostrings only. Source upstream functions without
# installing optional annotation retrieval and RNA-seq report dependencies.
source(file.path(args[[1]], "R", "get_seed.R"))
source(file.path(args[[1]], "R", "SeedMatchR.R"))
candidates <- read.delim(args[[2]], stringsAsFactors = FALSE, check.names = FALSE)
reference <- Biostrings::readDNAStringSet(args[[3]])
names(reference) <- sub("[[:space:]].*$", "", names(reference))
if (anyDuplicated(names(reference))) stop("Transcript FASTA identifiers must be unique")
seed.classes <- strsplit(args[[6]], ",", fixed = TRUE)[[1]]
output <- data.frame(candidate_id = character(), transcript_id = character(),
                     start = integer(), end = integer(), width = integer(),
                     g2_position = integer(), seed_class = character(),
                     target_sequence = character(), stringsAsFactors = FALSE)

# Reuse results when guides share a seed. Calling SeedMatchR's iranges branch
# avoids its unguarded zero-row data.frame conversion and preserves hit names.
cached <- new.env(parent = emptyenv())
chunks <- list()
for (i in seq_len(nrow(candidates))) {
  guide <- gsub("T", "U", candidates$guide[[i]], fixed = TRUE)
  for (seed in seed.classes) {
    pattern <- as.character(get_seed(guide.seq = guide, seed.name = seed)$Target.Seq)
    if (!exists(pattern, envir = cached, inherits = FALSE)) {
      ranges <- SeedMatchR(seqs = reference, sequence = guide, seed.name = seed,
                           res.format = "iranges", max.mismatch = 0L,
                           with.indels = FALSE, fixed = TRUE, allow_wobbles = FALSE)
      if (length(ranges) == 0L) {
        matches <- data.frame(transcript_id = character(), start = integer(),
                              end = integer(), width = integer())
      } else {
        matches <- data.frame(transcript_id = names(ranges),
                              start = BiocGenerics::start(ranges),
                              end = BiocGenerics::end(ranges),
                              width = BiocGenerics::width(ranges),
                              stringsAsFactors = FALSE)
      }
      assign(pattern, matches, envir = cached)
    }
    matches <- get(pattern, envir = cached, inherits = FALSE)
    if (nrow(matches) == 0L) next
    matches$candidate_id <- candidates$candidate_id[[i]]
    matches$seed_class <- seed
    matches$target_sequence <- pattern
    # A1 classes include a terminal A opposite g1, one base after the g2 anchor.
    matches$g2_position <- matches$end - as.integer(seed %in% c("mer8", "mer7A1"))
    chunks[[length(chunks) + 1L]] <- matches[names(output)]
  }
}
if (length(chunks)) output <- do.call(rbind, chunks)
write.table(output, args[[4]], sep = "\t", quote = TRUE, row.names = FALSE)
info <- data.frame(key = c("r_version", "biostrings_version", "reference_transcripts"),
                   value = c(as.character(getRversion()),
                             as.character(utils::packageVersion("Biostrings")),
                             as.character(length(reference))))
write.table(info, args[[5]], sep = "\t", quote = TRUE, row.names = FALSE)
