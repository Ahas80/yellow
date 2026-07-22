#!/usr/bin/env Rscript

# Prepare within-host variants for scientifically valid reference-aware plots.
#
# The raw show-snps coordinate is local to a reference contig.  This script
# verifies the complete provenance chain (show-snps -> reference FASTA -> GFF),
# resolves Prokka-renamed contigs only by exact sequence identity, and derives
# a cumulative coordinate separately within each reference assembly.  It never
# uses a same-length contig as sufficient evidence for a mapping.

suppressPackageStartupMessages({
  library(Biostrings)
  library(digest)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

normalise_sequence_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^>", "", x)
  sub("[[:space:]].*$", "", x)
}

sha256_file <- function(path) {
  unname(digest::digest(path, algo = "sha256", file = TRUE))
}

sha256_sequence <- function(sequence) {
  digest::digest(toupper(as.character(sequence)), algo = "sha256", serialize = FALSE)
}

collapse_unique <- function(x) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (!length(x)) NA_character_ else paste(sort(x), collapse = "; ")
}

contig_evidence_from_sequences <- function(headers, sequences, source) {
  if (!length(headers) || length(headers) != length(sequences)) {
    stop(source, " contains no usable FASTA records or mismatched headers/sequences.", call. = FALSE)
  }
  ids <- normalise_sequence_id(headers)
  if (any(!nzchar(ids)) || anyDuplicated(ids)) {
    duplicates <- unique(ids[duplicated(ids) | duplicated(ids, fromLast = TRUE)])
    stop(source, " has empty or duplicate normalized contig IDs: ",
         paste(duplicates, collapse = ", "), call. = FALSE)
  }
  sequences <- toupper(gsub("[[:space:]]", "", as.character(sequences)))
  invalid <- !grepl("^[ACGTMRWSYKVHDBN.-]+$", sequences)
  if (any(invalid)) {
    stop(source, " contains non-IUPAC sequence symbols in contig(s): ",
         paste(ids[invalid], collapse = ", "), call. = FALSE)
  }
  lengths <- nchar(sequences, type = "chars")
  if (any(lengths < 1L)) stop(source, " contains an empty contig.", call. = FALSE)

  tibble(
    contig_header = as.character(headers),
    contig_id = ids,
    contig_index = seq_along(ids),
    contig_length = as.numeric(lengths),
    contig_sha256 = vapply(sequences, sha256_sequence, character(1)),
    contig_offset = c(0, head(cumsum(as.numeric(lengths)), -1L)),
    sequence = sequences
  )
}

read_reference_fasta <- function(path) {
  dna <- Biostrings::readDNAStringSet(path, format = "fasta", use.names = TRUE)
  contig_evidence_from_sequences(names(dna), as.character(dna), basename(path))
}

parse_embedded_fasta <- function(lines, source) {
  header_idx <- which(startsWith(lines, ">"))
  if (!length(header_idx)) {
    stop(source, " has no embedded FASTA records after ##FASTA.", call. = FALSE)
  }
  next_header <- c(header_idx[-1L], length(lines) + 1L)
  headers <- sub("^>", "", lines[header_idx])
  sequences <- vapply(seq_along(header_idx), function(i) {
    first <- header_idx[i] + 1L
    last <- next_header[i] - 1L
    if (first > last) return("")
    paste(lines[first:last], collapse = "")
  }, character(1))
  contig_evidence_from_sequences(headers, sequences, source)
}

extract_gff_attribute <- function(attributes, key) {
  pattern <- paste0("(?:^|;)", key, "=([^;]*)")
  hit <- stringr::str_match(attributes, pattern)[, 2]
  hit <- ifelse(is.na(hit) | !nzchar(hit), NA_character_, vapply(hit, utils::URLdecode, character(1)))
  as.character(hit)
}

read_gff_bundle <- function(path) {
  lines <- readLines(path, warn = FALSE)
  fasta_marker <- which(lines == "##FASTA")
  if (length(fasta_marker) != 1L) {
    stop(basename(path), " must contain exactly one ##FASTA marker.", call. = FALSE)
  }
  marker <- fasta_marker[[1L]]
  gff_lines <- if (marker > 1L) lines[seq_len(marker - 1L)] else character()
  gff_lines <- gff_lines[nzchar(gff_lines) & !startsWith(gff_lines, "#")]
  fields <- strsplit(gff_lines, "\t", fixed = TRUE)
  field_count <- lengths(fields)
  if (any(field_count != 9L)) {
    stop(basename(path), " contains malformed non-comment GFF rows.", call. = FALSE)
  }
  feature_matrix <- if (length(fields)) do.call(rbind, fields) else matrix(character(), nrow = 0L, ncol = 9L)
  features <- as_tibble(feature_matrix, .name_repair = "minimal")
  names(features) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  features <- features %>%
    transmute(
      seqid = normalise_sequence_id(.data$seqid),
      source = .data$source,
      type = .data$type,
      start = suppressWarnings(as.numeric(.data$start)),
      end = suppressWarnings(as.numeric(.data$end)),
      strand = .data$strand,
      attributes = .data$attributes,
      gene = extract_gff_attribute(.data$attributes, "gene"),
      product = extract_gff_attribute(.data$attributes, "product"),
      locus_tag = extract_gff_attribute(.data$attributes, "locus_tag")
    )
  if (any(!is.finite(features$start) | !is.finite(features$end) |
          features$start < 1 | features$end < features$start)) {
    stop(basename(path), " contains invalid feature coordinates.", call. = FALSE)
  }

  fasta_lines <- if (marker < length(lines)) lines[(marker + 1L):length(lines)] else character()
  contigs <- parse_embedded_fasta(fasta_lines, paste0(basename(path), " embedded FASTA"))
  contig_lengths <- setNames(contigs$contig_length, contigs$contig_id)
  feature_limits <- unname(contig_lengths[features$seqid])
  if (any(is.na(feature_limits)) || any(features$end > feature_limits)) {
    bad <- unique(features$seqid[is.na(feature_limits) | features$end > feature_limits])
    stop(basename(path), " has features outside its embedded FASTA: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  list(contigs = contigs, features = features)
}

resolve_reference_to_gff_contigs <- function(reference_contigs, gff_contigs) {
  rows <- lapply(seq_len(nrow(reference_contigs)), function(i) {
    ref <- reference_contigs[i, ]
    exact_idx <- which(gff_contigs$contig_id == ref$contig_id)
    hash_idx <- which(gff_contigs$contig_sha256 == ref$contig_sha256 &
                        gff_contigs$contig_length == ref$contig_length)
    length_idx <- which(gff_contigs$contig_length == ref$contig_length)

    status <- "unmapped_no_sequence_identity"
    chosen <- NA_integer_
    if (length(exact_idx)) {
      exact_valid <- exact_idx[
        gff_contigs$contig_sha256[exact_idx] == ref$contig_sha256 &
          gff_contigs$contig_length[exact_idx] == ref$contig_length
      ]
      if (length(exact_valid) == 1L) {
        status <- "exact_normalized_id_and_sequence_sha256"
        chosen <- exact_valid
      } else if (length(exact_valid) > 1L) {
        status <- "ambiguous_exact_id"
      } else {
        status <- "exact_id_sequence_mismatch"
      }
    } else if (length(hash_idx) == 1L) {
      status <- "unique_sequence_sha256"
      chosen <- hash_idx
    } else if (length(hash_idx) > 1L) {
      status <- "ambiguous_sequence_sha256"
    }

    tibble(
      Reference_Contig_ID = ref$contig_id,
      Reference_Contig_Index = ref$contig_index,
      Reference_Contig_Length = ref$contig_length,
      Reference_Contig_SHA256 = ref$contig_sha256,
      Reference_Contig_Offset = ref$contig_offset,
      GFF_Contig_ID = if (is.na(chosen)) NA_character_ else gff_contigs$contig_id[chosen],
      GFF_Contig_Length = if (is.na(chosen)) NA_real_ else gff_contigs$contig_length[chosen],
      GFF_Contig_SHA256 = if (is.na(chosen)) NA_character_ else gff_contigs$contig_sha256[chosen],
      Contig_Mapping_Method = status,
      Contig_Mapping_Valid = status %in% c(
        "exact_normalized_id_and_sequence_sha256", "unique_sequence_sha256"
      ),
      Length_Only_Candidate_Count = length(length_idx)
    )
  })
  mapping <- bind_rows(rows)
  duplicated_target <- !is.na(mapping$GFF_Contig_ID) &
    (duplicated(mapping$GFF_Contig_ID) | duplicated(mapping$GFF_Contig_ID, fromLast = TRUE))
  if (any(duplicated_target)) {
    mapping$Contig_Mapping_Method[duplicated_target] <- "ambiguous_non_unique_gff_target"
    mapping$Contig_Mapping_Valid[duplicated_target] <- FALSE
  }
  mapping
}

read_show_snps <- function(path) {
  raw <- readr::read_tsv(
    path, col_names = FALSE, col_types = cols(.default = col_character()),
    show_col_types = FALSE, progress = FALSE
  )
  if (ncol(raw) < 12L) stop(basename(path), " has fewer than 12 show-snps columns.", call. = FALSE)
  tibble(
    Pos_Ref = suppressWarnings(as.numeric(raw[[1L]])),
    Ref_Base = as.character(raw[[2L]]),
    Qry_Base = as.character(raw[[3L]]),
    Pos_Qry = suppressWarnings(as.numeric(raw[[4L]])),
    Show_Snps_Reference_Contig_Length = suppressWarnings(as.numeric(raw[[7L]])),
    Show_Snps_Query_Contig_Length = suppressWarnings(as.numeric(raw[[8L]])),
    Ref_Seqid = normalise_sequence_id(raw[[11L]]),
    Qry_Seqid = normalise_sequence_id(raw[[12L]]),
    Show_Snps_Source_Line = seq_len(nrow(raw))
  )
}

variant_key <- function(pos_ref, ref_base, qry_base, pos_qry, ref_seqid, qry_seqid) {
  paste(
    format(as.numeric(pos_ref), scientific = FALSE, trim = TRUE),
    as.character(ref_base), as.character(qry_base),
    format(as.numeric(pos_qry), scientific = FALSE, trim = TRUE),
    normalise_sequence_id(ref_seqid), normalise_sequence_id(qry_seqid),
    sep = "\037"
  )
}

attach_show_snps_evidence <- function(variants) {
  variants$Show_Snps_Match_Count <- NA_integer_
  variants$Show_Snps_Source_Line <- NA_integer_
  variants$Show_Snps_Reference_Contig_Length <- NA_real_
  variants$Show_Snps_Query_Contig_Length <- NA_real_

  path_hashes <- variants %>% distinct(.data$SNP_Path, .data$SNP_SHA256)
  if (anyDuplicated(path_hashes$SNP_Path)) {
    stop("A show-snps path is associated with more than one recorded SHA-256.", call. = FALSE)
  }
  for (i in seq_len(nrow(path_hashes))) {
    path <- normalizePath(path_hashes$SNP_Path[i], winslash = "/", mustWork = TRUE)
    expected_hash <- path_hashes$SNP_SHA256[i]
    observed_hash <- sha256_file(path)
    if (!identical(observed_hash, expected_hash)) {
      stop("show-snps SHA-256 mismatch: ", path, call. = FALSE)
    }
    raw <- read_show_snps(path)
    raw_key <- variant_key(raw$Pos_Ref, raw$Ref_Base, raw$Qry_Base, raw$Pos_Qry,
                           raw$Ref_Seqid, raw$Qry_Seqid)
    idx <- which(variants$SNP_Path == path_hashes$SNP_Path[i] |
                   normalizePath(variants$SNP_Path, winslash = "/", mustWork = FALSE) == path)
    input_key <- variant_key(
      variants$Pos_Ref[idx], variants$Ref_Base[idx], variants$Qry_Base[idx],
      variants$Pos_Qry[idx], variants$Ref_Seqid[idx], variants$Qry_Seqid[idx]
    )
    if (anyDuplicated(input_key)) {
      stop("Input variants duplicate an exact show-snps row in ", basename(path), ".", call. = FALSE)
    }
    matches <- lapply(input_key, function(key) which(raw_key == key))
    variants$Show_Snps_Match_Count[idx] <- lengths(matches)
    unique_match <- lengths(matches) == 1L
    if (any(unique_match)) {
      raw_idx <- vapply(matches[unique_match], `[[`, integer(1), 1L)
      target <- idx[unique_match]
      variants$Show_Snps_Source_Line[target] <- raw$Show_Snps_Source_Line[raw_idx]
      variants$Show_Snps_Reference_Contig_Length[target] <- raw$Show_Snps_Reference_Contig_Length[raw_idx]
      variants$Show_Snps_Query_Contig_Length[target] <- raw$Show_Snps_Query_Contig_Length[raw_idx]
    }
  }
  variants
}

annotate_one_position <- function(features, gff_seqid, position) {
  hits <- features %>%
    filter(.data$seqid == gff_seqid, .data$start <= position, .data$end >= position)
  preferred <- if (any(hits$type == "CDS")) {
    hits %>% filter(.data$type == "CDS")
  } else if (any(hits$type == "gene")) {
    hits %>% filter(.data$type == "gene")
  } else {
    hits[0, ]
  }
  if (!nrow(preferred)) {
    return(tibble(
      Region = "Intergenic", Gene = NA_character_, Product = NA_character_,
      Locus = NA_character_, Overlapping_Feature_Count = 0L,
      Annotation_Quality_Flag = "validated_intergenic"
    ))
  }
  tibble(
    Region = if (any(preferred$type == "CDS")) "CDS" else "gene",
    Gene = collapse_unique(preferred$gene),
    Product = collapse_unique(preferred$product),
    Locus = collapse_unique(preferred$locus_tag),
    Overlapping_Feature_Count = nrow(preferred),
    Annotation_Quality_Flag = if (nrow(preferred) == 1L) {
      if (all(is.na(preferred$gene)) && all(is.na(preferred$product))) {
        "validated_feature_missing_gene_and_product"
      } else {
        "validated_unique_feature"
      }
    } else {
      "validated_multiple_overlapping_features"
    }
  )
}

prepare_reference_aware_variants <- function(
    input_file = file.path("results", "longitudinal", "variant_annotation_detailed.csv"),
    output_dir = file.path("results", "figure_audit"),
    strict = TRUE) {
  required_packages <- c("Biostrings", "digest", "dplyr", "readr", "stringr", "tibble")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
  if (!file.exists(input_file)) stop("Missing detailed variant input: ", input_file, call. = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  output_csv <- file.path(output_dir, "reference_aware_variants.csv")
  validation_csv <- file.path(output_dir, "reference_aware_variants_validation.csv")
  validation_txt <- file.path(output_dir, "reference_aware_variants_validation.txt")

  variants <- readr::read_csv(
    input_file, show_col_types = FALSE,
    col_types = cols(Participant_id = col_character(), .default = col_guess())
  )
  required <- c(
    "Pos_Ref", "Ref_Base", "Qry_Base", "Pos_Qry", "Ref_Seqid", "Qry_Seqid",
    "Participant_id", "From_Time", "To_Time", "SNP_Path", "SNP_SHA256",
    "Reference_FASTA_Path", "Reference_FASTA_SHA256", "Type",
    "GFF_Path", "GFF_SHA256", "Region", "Gene", "Product", "Locus"
  )
  missing <- setdiff(required, names(variants))
  if (length(missing)) stop("Detailed variant input lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!nrow(variants)) stop("Detailed variant input is empty; no coordinate map can be validated.", call. = FALSE)

  variants <- variants %>%
    mutate(
      Variant_Row_ID = row_number(),
      Participant_id = as.character(.data$Participant_id),
      Ref_Seqid_Normalized = normalise_sequence_id(.data$Ref_Seqid),
      Qry_Seqid_Normalized = normalise_sequence_id(.data$Qry_Seqid),
      Original_Region = as.character(.data$Region),
      Original_Gene = as.character(.data$Gene),
      Original_Product = as.character(.data$Product),
      Original_Locus = as.character(.data$Locus)
    )
  variants <- attach_show_snps_evidence(variants)
  variants$Reference_ID <- NA_character_
  variants$Reference_Total_Length <- NA_real_
  variants$Reference_Contig_ID <- NA_character_
  variants$Reference_Contig_Index <- NA_integer_
  variants$Reference_Contig_Length <- NA_real_
  variants$Reference_Contig_SHA256 <- NA_character_
  variants$Reference_Contig_Offset <- NA_real_
  variants$Reference_Cumulative_Position <- NA_real_
  variants$GFF_Contig_ID <- NA_character_
  variants$Contig_Mapping_Method <- NA_character_
  variants$Contig_Mapping_Valid <- FALSE
  variants$Position_Within_Contig_Valid <- FALSE
  variants$Reference_Base_From_FASTA <- NA_character_
  variants$Reference_Base_Matches_FASTA <- NA
  variants$Show_Snps_Length_Matches_FASTA <- FALSE
  variants$Overlapping_Feature_Count <- NA_integer_
  variants$Annotation_Quality_Flag <- "unvalidated"

  original_names <- c("Region", "Gene", "Product", "Locus")
  variants[original_names] <- lapply(variants[original_names], function(x) NA_character_)
  validation_rows <- list()
  reference_sets <- variants %>%
    distinct(.data$Reference_FASTA_Path, .data$Reference_FASTA_SHA256,
             .data$GFF_Path, .data$GFF_SHA256)
  if (anyDuplicated(reference_sets$Reference_FASTA_Path)) {
    stop("A reference FASTA maps to multiple path/hash/GFF provenance records.", call. = FALSE)
  }

  for (i in seq_len(nrow(reference_sets))) {
    ref_path <- normalizePath(reference_sets$Reference_FASTA_Path[i], winslash = "/", mustWork = TRUE)
    gff_path <- normalizePath(reference_sets$GFF_Path[i], winslash = "/", mustWork = TRUE)
    ref_hash_ok <- identical(sha256_file(ref_path), reference_sets$Reference_FASTA_SHA256[i])
    gff_hash_ok <- identical(sha256_file(gff_path), reference_sets$GFF_SHA256[i])
    if (!ref_hash_ok) stop("Reference FASTA SHA-256 mismatch: ", ref_path, call. = FALSE)
    if (!gff_hash_ok) stop("GFF SHA-256 mismatch: ", gff_path, call. = FALSE)

    ref_contigs <- read_reference_fasta(ref_path)
    gff_bundle <- read_gff_bundle(gff_path)
    mapping <- resolve_reference_to_gff_contigs(ref_contigs, gff_bundle$contigs)
    reference_total <- sum(ref_contigs$contig_length)
    reference_id <- paste0(tools::file_path_sans_ext(basename(ref_path)), "__", substr(reference_sets$Reference_FASTA_SHA256[i], 1L, 12L))
    idx <- which(variants$Reference_FASTA_Path == reference_sets$Reference_FASTA_Path[i])

    for (j in idx) {
      ref_match <- which(ref_contigs$contig_id == variants$Ref_Seqid_Normalized[j])
      if (length(ref_match) != 1L) next
      rc <- ref_contigs[ref_match, ]
      map_row <- mapping %>% filter(.data$Reference_Contig_ID == rc$contig_id)
      pos <- suppressWarnings(as.numeric(variants$Pos_Ref[j]))
      position_valid <- length(pos) == 1L && is.finite(pos) && pos == floor(pos) && pos >= 1 && pos <= rc$contig_length

      variants$Reference_ID[j] <- reference_id
      variants$Reference_Total_Length[j] <- reference_total
      variants$Reference_Contig_ID[j] <- rc$contig_id
      variants$Reference_Contig_Index[j] <- rc$contig_index
      variants$Reference_Contig_Length[j] <- rc$contig_length
      variants$Reference_Contig_SHA256[j] <- rc$contig_sha256
      variants$Reference_Contig_Offset[j] <- rc$contig_offset
      variants$Reference_Cumulative_Position[j] <- if (position_valid) rc$contig_offset + pos else NA_real_
      variants$GFF_Contig_ID[j] <- map_row$GFF_Contig_ID
      variants$Contig_Mapping_Method[j] <- map_row$Contig_Mapping_Method
      variants$Contig_Mapping_Valid[j] <- map_row$Contig_Mapping_Valid
      variants$Position_Within_Contig_Valid[j] <- position_valid
      variants$Show_Snps_Length_Matches_FASTA[j] <-
        isTRUE(variants$Show_Snps_Match_Count[j] == 1L) &&
        isTRUE(variants$Show_Snps_Reference_Contig_Length[j] == rc$contig_length)

      if (position_valid) {
        observed_base <- substr(rc$sequence, pos, pos)
        variants$Reference_Base_From_FASTA[j] <- observed_base
        if (identical(as.character(variants$Ref_Base[j]), ".")) {
          variants$Reference_Base_Matches_FASTA[j] <- NA
        } else {
          variants$Reference_Base_Matches_FASTA[j] <-
            identical(toupper(as.character(variants$Ref_Base[j])), observed_base)
        }
      }

      if (position_valid && isTRUE(map_row$Contig_Mapping_Valid)) {
        ann <- annotate_one_position(gff_bundle$features, map_row$GFF_Contig_ID, pos)
        variants$Region[j] <- ann$Region
        variants$Gene[j] <- ann$Gene
        variants$Product[j] <- ann$Product
        variants$Locus[j] <- ann$Locus
        variants$Overlapping_Feature_Count[j] <- ann$Overlapping_Feature_Count
        variants$Annotation_Quality_Flag[j] <- ann$Annotation_Quality_Flag
      } else if (!position_valid) {
        variants$Annotation_Quality_Flag[j] <- "invalid_reference_coordinate"
      } else {
        variants$Annotation_Quality_Flag[j] <- "unvalidated_gff_contig_mapping"
      }
    }

    per_contig <- mapping %>%
      mutate(
        Reference_ID = reference_id,
        Reference_FASTA_Path = ref_path,
        Reference_FASTA_SHA256 = reference_sets$Reference_FASTA_SHA256[i],
        Reference_FASTA_Hash_Valid = ref_hash_ok,
        GFF_Path = gff_path,
        GFF_SHA256 = reference_sets$GFF_SHA256[i],
        GFF_Hash_Valid = gff_hash_ok,
        Reference_Total_Length = reference_total,
        Variant_Count = 0L,
        All_Positions_Valid = TRUE,
        All_Reference_Bases_Valid = TRUE,
        All_Show_Snps_Rows_Validated = TRUE,
        All_Show_Snps_Lengths_Match = TRUE
      )
    for (k in seq_len(nrow(per_contig))) {
      used <- idx[variants$Ref_Seqid_Normalized[idx] == per_contig$Reference_Contig_ID[k]]
      per_contig$Variant_Count[k] <- length(used)
      per_contig$All_Positions_Valid[k] <- if (!length(used)) TRUE else all(variants$Position_Within_Contig_Valid[used])
      per_contig$All_Reference_Bases_Valid[k] <- if (!length(used)) TRUE else all(
        is.na(variants$Reference_Base_Matches_FASTA[used]) |
          variants$Reference_Base_Matches_FASTA[used]
      )
      per_contig$All_Show_Snps_Rows_Validated[k] <- if (!length(used)) TRUE else all(variants$Show_Snps_Match_Count[used] == 1L)
      per_contig$All_Show_Snps_Lengths_Match[k] <- if (!length(used)) TRUE else all(variants$Show_Snps_Length_Matches_FASTA[used])
    }
    validation_rows[[length(validation_rows) + 1L]] <- per_contig
  }

  variants <- variants %>%
    mutate(
      Variant_Validation_Status = case_when(
        .data$Show_Snps_Match_Count != 1L ~ "fail_show_snps_row_not_unique",
        is.na(.data$Reference_Contig_Index) ~ "fail_reference_contig_not_exactly_resolved",
        !.data$Position_Within_Contig_Valid ~ "fail_position_outside_reference_contig",
        !.data$Show_Snps_Length_Matches_FASTA ~ "fail_show_snps_length_mismatch",
        !is.na(.data$Reference_Base_Matches_FASTA) & !.data$Reference_Base_Matches_FASTA ~ "fail_reference_base_mismatch",
        !.data$Contig_Mapping_Valid ~ "fail_gff_contig_mapping_unvalidated",
        TRUE ~ "pass"
      ),
      Figure_Eligible = .data$Variant_Validation_Status == "pass"
    ) %>%
    relocate(
      Variant_Row_ID, Reference_ID, Participant_id,
      From_Time, To_Time, Reference_Cumulative_Position,
      Reference_Total_Length, Reference_Contig_ID,
      Reference_Contig_Index, Reference_Contig_Length,
      Reference_Contig_Offset, Pos_Ref
    )

  validation <- bind_rows(validation_rows) %>%
    mutate(
      Figure_Eligible = .data$Reference_FASTA_Hash_Valid & .data$GFF_Hash_Valid &
        .data$Contig_Mapping_Valid & .data$All_Positions_Valid &
        .data$All_Reference_Bases_Valid & .data$All_Show_Snps_Rows_Validated &
        .data$All_Show_Snps_Lengths_Match
    ) %>%
    select(all_of(c(
      "Reference_ID", "Reference_FASTA_Path", "Reference_FASTA_SHA256",
      "Reference_FASTA_Hash_Valid", "Reference_Total_Length",
      "Reference_Contig_ID", "Reference_Contig_Index",
      "Reference_Contig_Length", "Reference_Contig_Offset",
      "Reference_Contig_SHA256", "GFF_Path", "GFF_SHA256",
      "GFF_Hash_Valid", "GFF_Contig_ID", "GFF_Contig_Length",
      "GFF_Contig_SHA256", "Contig_Mapping_Method",
      "Contig_Mapping_Valid", "Length_Only_Candidate_Count",
      "Variant_Count", "All_Positions_Valid",
      "All_Reference_Bases_Valid", "All_Show_Snps_Rows_Validated",
      "All_Show_Snps_Lengths_Match", "Figure_Eligible"
    )))

  all_valid <- all(variants$Figure_Eligible) && all(validation$Figure_Eligible)
  input_hash <- sha256_file(input_file)
  readr::write_csv(variants, output_csv, na = "")
  readr::write_csv(validation, validation_csv, na = "")

  status <- if (all_valid) "RETAIN_WITH_REFERENCE_AWARE_REPLOT" else "DO_NOT_RETAIN"
  annotation_counts <- variants %>% count(.data$Annotation_Quality_Flag, name = "n")
  mapping_counts <- validation %>% count(.data$Contig_Mapping_Method, name = "n")
  report <- c(
    "Reference-aware within-host variant validation",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste0("Input: ", normalizePath(input_file, winslash = "/", mustWork = TRUE)),
    paste0("Input SHA-256: ", input_hash),
    paste0("Variant rows: ", nrow(variants)),
    paste0("Reference assemblies: ", n_distinct(variants$Reference_ID)),
    paste0("Reference contigs: ", nrow(validation)),
    paste0("Validated variant rows: ", sum(variants$Figure_Eligible), "/", nrow(variants)),
    paste0("Unique show-snps row matches: ", sum(variants$Show_Snps_Match_Count == 1L), "/", nrow(variants)),
    paste0("Coordinates within actual reference contigs: ", sum(variants$Position_Within_Contig_Valid), "/", nrow(variants)),
    paste0("Reference bases matching FASTA (where applicable): ",
           sum(variants$Reference_Base_Matches_FASTA %in% TRUE, na.rm = TRUE), "/",
           sum(!is.na(variants$Reference_Base_Matches_FASTA))),
    paste0("GFF contig mappings validated by ID+sequence or unique sequence SHA-256: ",
           sum(validation$Contig_Mapping_Valid), "/", nrow(validation)),
    "",
    "Contig mapping methods:",
    paste0("- ", mapping_counts$Contig_Mapping_Method, ": ", mapping_counts$n),
    "",
    "Annotation quality flags:",
    paste0("- ", annotation_counts$Annotation_Quality_Flag, ": ", annotation_counts$n),
    "",
    paste0("FIG08_STATUS=", status),
    paste0("FIG08_VALIDATABLE=", ifelse(all_valid, "TRUE", "FALSE")),
    "FIG08_COORDINATE=Reference_Cumulative_Position",
    "FIG08_REQUIREMENT=Facet by reference/comparison; coordinates from different references are not homologous and must not share an implied common genome axis.",
    "FIG08_REQUIREMENT=Use Reference_Total_Length and per-reference contig boundaries; do not use raw Pos_Ref as a genome-wide coordinate.",
    "FIG08_REQUIREMENT=Only Figure_Eligible rows may be plotted."
  )
  writeLines(report, validation_txt, useBytes = TRUE)

  if (isTRUE(strict) && !all_valid) {
    stop("Reference-aware variant validation failed. See ", validation_txt, call. = FALSE)
  }
  invisible(list(
    variants = variants, validation = validation, all_valid = all_valid,
    output_csv = output_csv, validation_csv = validation_csv,
    validation_txt = validation_txt
  ))
}

parse_cli <- function(args) {
  values <- list(
    input_file = file.path("results", "longitudinal", "variant_annotation_detailed.csv"),
    output_dir = file.path("results", "figure_audit")
  )
  for (arg in args) {
    if (startsWith(arg, "--input=")) values$input_file <- sub("^--input=", "", arg)
    else if (startsWith(arg, "--output-dir=")) values$output_dir <- sub("^--output-dir=", "", arg)
    else stop("Unknown argument: ", arg, call. = FALSE)
  }
  values
}

if (sys.nframe() == 0L) {
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))
  result <- prepare_reference_aware_variants(cli$input_file, cli$output_dir, strict = TRUE)
  message("Reference-aware variants written to ", result$output_csv)
  message("Validation status: ", ifelse(result$all_valid, "PASS", "FAIL"))
}
