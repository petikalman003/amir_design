#!/usr/bin/env Rscript
# Call the pinned upstream SeedMatchR sequence API on lossless prefiltered sets.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11L) {
  stop("Expected repo, candidates, reference, sites, provenance, seeds, mapping, target, summary, progress, and detail mode")
}
if (getRversion() < "4.3.0") stop("SeedMatchR requires R >= 4.3.0")
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Biostrings is missing; run scripts/install_r.sh")
}
source(file.path(args[[1]], "R", "get_seed.R"))
source(file.path(args[[1]], "R", "SeedMatchR.R"))
notify <- function(message) cat(message, "\n", file = args[[10]], append = TRUE, sep = "")
notify("Loading the full transcript reference")
candidates <- read.delim(args[[2]], colClasses = "character", na.strings = character(),
                         stringsAsFactors = FALSE, check.names = FALSE)
candidates <- candidates[order(candidates$candidate_id, method = "radix"), , drop = FALSE]
reference <- Biostrings::readDNAStringSet(args[[3]])
names(reference) <- sub("[[:space:]].*$", "", names(reference))
if (anyDuplicated(names(reference))) stop("Transcript FASTA identifiers must be unique")
reference.count <- length(reference)
mapping <- read.delim(args[[7]], colClasses = "character", na.strings = character(),
                         stringsAsFactors = FALSE, check.names = FALSE)
gene.index <- match(names(reference), mapping$transcript_id)
if (anyNA(gene.index)) stop("A reference transcript has no gene mapping")
genes <- mapping$gene_id[gene.index]
# The intended gene contributes no off-targets: exclude all its transcripts first.
reference <- reference[genes != args[[8]]]
genes <- genes[genes != args[[8]]]
seed.classes <- strsplit(args[[6]], ",", fixed = TRUE)[[1]]
priority <- c(mer8 = 0L, mer7m8 = 1L, mer7A1 = 2L, mer6 = 3L)
if (!length(seed.classes) || any(!seed.classes %in% names(priority))) stop("Unsupported seed classes")
output <- data.frame(candidate_id = character(), transcript_id = character(),
                     start = integer(), end = integer(), width = integer(),
                     g2_position = integer(), seed_class = character(),
                     target_sequence = character(), gene_id = character(),
                     stringsAsFactors = FALSE)
include.sites <- identical(args[[11]], "sites")
write.table(output, args[[4]], sep = "\t", quote = TRUE, row.names = FALSE)
summary.names <- c("candidate_id", "off_target_genes", "off_target_transcripts", "off_target_sites",
                   paste0(seed.classes, "_sites"))
summary <- as.data.frame(setNames(rep(list(integer()), length(summary.names)), summary.names))
summary$candidate_id <- character()
write.table(summary, args[[9]], sep = "\t", quote = TRUE, row.names = FALSE)

# Index every distinct requested target pattern together, grouped by width for
# PDict's exact dictionary API. The index only filters transcripts with no match;
# all final coordinates still come from the unmodified SeedMatchR function.
guides <- gsub("T", "U", candidates$guide, fixed = TRUE)
patterns <- lapply(guides, function(guide) {
  vapply(seed.classes, function(seed) {
    as.character(get_seed(guide.seq = guide, seed.name = seed)$Target.Seq)
  }, character(1))
})
unique.patterns <- unique(unlist(patterns, use.names = FALSE))
pattern.index <- lapply(patterns, match, table = unique.patterns)
selected.transcripts <- vector("list", length(unique.patterns))
if (length(reference)) {
  notify(sprintf("Indexing %d unique exact seeds across %d non-target transcripts",
                 length(unique.patterns), length(reference)))
  for (width in unique(nchar(unique.patterns))) {
    group <- which(nchar(unique.patterns) == width)
    dictionary <- Biostrings::PDict(Biostrings::DNAStringSet(unique.patterns[group]))
    # At most four million integer cells (16 MB), even for large guide sets.
    chunk.size <- min(25000L, max(1L, as.integer(4000000 / length(group))))
    for (start in seq.int(1L, length(reference), by = chunk.size)) {
      indices <- seq.int(start, min(start + chunk.size - 1L, length(reference)))
      counts <- Biostrings::vcountPDict(dictionary, reference[indices], fixed = TRUE)
      for (j in seq_along(group)) {
        found <- indices[which(counts[j, ] > 0L)]
        if (length(found)) {
          selected.transcripts[[group[[j]]]] <- c(selected.transcripts[[group[[j]]]], found)
        }
      }
    }
  }
}
notify("Exact seed prefilter ready; screening guides")
for (i in seq_len(nrow(candidates))) {
  chunks <- vector("list", length(seed.classes))
  for (j in seq_along(seed.classes)) {
    seed <- seed.classes[[j]]
    pattern <- patterns[[i]][[j]]
    subset.indices <- selected.transcripts[[pattern.index[[i]][[j]]]]
    if (!length(subset.indices)) next
    subset <- reference[subset.indices]
    ranges <- SeedMatchR(seqs = subset, sequence = guides[[i]], seed.name = seed,
                         res.format = "iranges", max.mismatch = 0L,
                         with.indels = FALSE, fixed = TRUE, allow_wobbles = FALSE)
    if (!length(ranges)) next
    matches <- data.frame(candidate_id = candidates$candidate_id[[i]],
                          transcript_id = names(ranges),
                          start = BiocGenerics::start(ranges),
                          end = BiocGenerics::end(ranges),
                          width = BiocGenerics::width(ranges),
                          seed_class = seed, target_sequence = pattern,
                          stringsAsFactors = FALSE)
    matches$g2_position <- matches$end - as.integer(seed %in% c("mer8", "mer7A1"))
    matches$gene_id <- genes[subset.indices[match(matches$transcript_id, names(subset))]]
    chunks[[j]] <- matches[names(output)]
  }
  # Retain only the strongest class at each g2 location, then release this guide's
  # rows after writing. Do not accumulate millions of sites across all candidates.
  hits <- if (all(vapply(chunks, is.null, logical(1)))) output else do.call(rbind, chunks)
  if (nrow(hits)) {
    hits <- hits[order(priority[hits$seed_class]), ]
    key <- paste(hits$transcript_id, hits$g2_position, sep = ":")
    hits <- hits[!duplicated(key), ]
    if (include.sites) {
      hits <- hits[order(hits$transcript_id, hits$g2_position, method = "radix"), ]
      write.table(hits, args[[4]], sep = "\t", quote = TRUE, row.names = FALSE,
                  col.names = FALSE, append = TRUE)
    }
  }
  counts <- c(length(unique(hits$gene_id)), length(unique(hits$transcript_id)), nrow(hits),
              vapply(seed.classes, function(seed) sum(hits$seed_class == seed), integer(1)))
  row <- as.data.frame(as.list(c(candidate_id = candidates$candidate_id[[i]],
                                  setNames(counts, summary.names[-1]))), stringsAsFactors = FALSE)
  write.table(row, args[[9]], sep = "\t", quote = TRUE, row.names = FALSE,
              col.names = FALSE, append = TRUE)
  notify(sprintf("Screened guide %d/%d (%d off-target sites)", i, nrow(candidates), nrow(hits)))
}
info <- data.frame(key = c("r_version", "biostrings_version", "reference_transcripts", "search_strategy"),
                   value = c(as.character(getRversion()),
                             as.character(utils::packageVersion("Biostrings")),
                             as.character(reference.count),
                             "lossless batched exact seed prefilter; upstream SeedMatchR per canonical class"))
write.table(info, args[[5]], sep = "\t", quote = TRUE, row.names = FALSE)
