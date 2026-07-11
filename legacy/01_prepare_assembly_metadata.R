#!/usr/bin/env Rscript
if (!identical(Sys.getenv("RUTIS_ALLOW_OBSOLETE_LEGACY_STATUS"), "1")) {
  stop(
    "legacy/01_prepare_assembly_metadata.R uses obsolete ASB/UTI/Negative ",
    "strict-CFU status logic. Use 00a_load_clean_clinical.R + ",
    "00b_classify_episodes.R for the catheter-aware UTI vs Not_UTI definition, ",
    "or set RUTIS_ALLOW_OBSOLETE_LEGACY_STATUS=1 to run this legacy script intentionally.",
    call. = FALSE
  )
}
# =============================================================
# 01_prepare_assembly_metadata.R — Population-first S&S + strict CFU + QA
# -------------------------------------------------------------
# Culture positive:
#   - CFU text explicitly "> / ≥ / 'meer dan' 100,000" => POSITIVE
#   - If CFU missing, Beoord "+++" => POSITIVE; "+" or "++" => NEGATIVE
#   - Any ranges (e.g., "10k-100k", "100000–200000") => NOT positive
#   - Single numerics like "200000" are NOT positive by default (toggle below)
#
# Signs & Symptoms (S&S) resolution order (ALL batches):
#   1) Population (authoritative if determinable)
#        - "pt zonder UWI"                    -> FALSE
#        - "pt met UWI voor/na meting"        -> FALSE at THIS timepoint
#        - "pt met UWI" (no "voor/na meting") -> TRUE
#      otherwise…
#   2) Explicit symptom columns (any present => TRUE)
#   3) SnS_status fallback (B1–B2 only): 2=TRUE, 0=FALSE, else=UNKNOWN
#
# Final status:
#   culture_pos + S&S TRUE   -> "UTI"
#   culture_pos + S&S FALSE  -> "ASB"
#   culture_pos + S&S NA     -> "Culture-positive, S&S unknown"
#   culture negative         -> "Negative"
#   culture indeterminate    -> "None"
# =============================================================

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tidyr); library(readr)
  library(Biostrings); library(glue); library(forcats); library(ggplot2); library(purrr)
})

# ---------- helpers ----------
nz_chr <- function(x) ifelse(is.na(x), "", x)

# Canonicalize timepoints to "T<number>" or "Uricult"
canon_tp <- function(x) {
  x   <- stringr::str_trim(as.character(x))
  low <- stringr::str_to_lower(x)
  d   <- suppressWarnings(readr::parse_number(low))
  d   <- ifelse(is.na(d), NA_integer_, as.integer(d))
  dplyr::case_when(
    grepl("^t?\\d+$", low) ~ paste0("T", d),
    grepl("uricult", low)  ~ "Uricult",
    TRUE                   ~ x
  )
}

# ---- CFU parsing controls (strict) ----
CFU_THRESHOLD <- 100000L
CFU_ACCEPT_SINGLE_NUMERIC_GT <- FALSE  # set TRUE only if plain "200000" should count

th_gt <- function(v) v > CFU_THRESHOLD

# Robust CFU bucketer (strict, vectorized safely)
# - EXPLICIT "meer dan 100000" OR symbols ">" / "≥" 100000 => ">=1e5"
# - ANY RANGE "x-y" => "<1e5_or_other"
# - SINGLE NUMERIC >100000 => ">=1e5" only if toggle TRUE; else "<1e5_or_other"
# - non-empty everything else => "<1e5_or_other"
cfu_bucket <- function(x) {
  s <- tolower(stringr::str_trim(as.character(x)))
  s_dotless <- gsub("\\.", "", s)  # strip thousand separators
  out <- rep(NA_character_, length(s))
  
  # explicit textual/symbolic "more than 100000"
  explicit_meer   <- grepl("meer\\s*dan\\s*100\\s*000", s_dotless)
  explicit_symbol <- grepl("(>|≥|>=)\\s*100\\s*000",   s_dotless)
  explicit_pos    <- explicit_meer | explicit_symbol
  out[explicit_pos] <- ">=1e5"
  
  # ranges → NOT positive
  idx_range <- which(grepl("[-–]", s_dotless))
  if (length(idx_range)) {
    to_fill <- idx_range[is.na(out[idx_range])]
    out[to_fill] <- "<1e5_or_other"
  }
  
  # single numeric (no hyphen)
  idx_need <- which(is.na(out) & grepl("\\d", s_dotless) & !grepl("[-–]", s_dotless))
  if (length(idx_need)) {
    nums <- regmatches(s_dotless[idx_need], gregexpr("\\d+", s_dotless[idx_need]))
    last_num <- vapply(nums, function(v) if (length(v)) tail(v,1) else NA_character_, character(1))
    val <- suppressWarnings(as.integer(last_num))
    ok  <- !is.na(val)
    if (any(ok)) {
      pos <- idx_need[ ok & CFU_ACCEPT_SINGLE_NUMERIC_GT & th_gt(val[ok]) ]
      neg <- idx_need[ ok & (!CFU_ACCEPT_SINGLE_NUMERIC_GT | !th_gt(val[ok])) ]
      if (length(pos)) out[pos] <- ">=1e5"
      if (length(neg)) out[neg] <- "<1e5_or_other"
    }
  }
  
  # everything else that’s non-empty → other
  still_na <- is.na(out) & nz_chr(s_dotless) != ""
  out[still_na] <- "<1e5_or_other"
  out
}

# Batch-3 Population mapping (also used for all batches if present)
population_to_sns <- function(pop) {
  p <- stringr::str_to_lower(nz_chr(pop))
  has_zonder  <- stringr::str_detect(p, "\\bzonder\\s*uwi\\b")
  has_met     <- stringr::str_detect(p, "\\bmet\\s*uwi\\b")
  has_voorna  <- stringr::str_detect(p, "(voor\\s*/?\\s*na|voor-?na)")
  has_meting  <- stringr::str_detect(p, "\\bmeting\\b")
  voorna_meting <- has_voorna & has_meting
  dplyr::case_when(
    has_zonder               ~ FALSE,  # no S&S
    has_met & voorna_meting  ~ FALSE,  # no S&S at THIS timepoint
    has_met                  ~ TRUE,   # S&S present
    TRUE                     ~ NA
  )
}

# Always-return schema for SnS
empty_sns_frame <- function() {
  tibble(Participant_id=character(), Timepoint=character(),
         SnS_status=integer(), has_SnS=logical(), Batch=integer())
}

# Tri-state episode collapse: all TRUE -> TRUE; all FALSE -> FALSE; mix -> NA
collapse_tristate <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA)
  ux <- unique(x)
  if (length(ux) == 1) ux else NA
}

# ---- OUTPUT CONTROL & RELIABLE WRITES ----
base_dir <- normalizePath(getwd(), mustWork = TRUE)
out_dir  <- file.path(base_dir, "results")
plot_dir <- file.path(out_dir, "plots")
dbg_dir  <- file.path(out_dir, "debug")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dbg_dir,  recursive = TRUE, showWarnings = FALSE)

msg <- function(fmt, ...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), sprintf(fmt, ...), "\n")
  flush.console()
}
write_csv_confirm <- function(x, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp"); readr::write_csv(x, tmp, ...)
  if (file.exists(path)) unlink(path)
  ok <- file.rename(tmp, path)
  fi <- tryCatch(file.info(path), error = function(e) NULL)
  msg("wrote: %s  (size=%.1f KB, mtime=%s)",
      normalizePath(path), if (!is.null(fi)) fi$size/1024 else NA,
      if (!is.null(fi)) as.character(fi$mtime) else "NA")
  invisible(ok)
}
write_csv <- write_csv_confirm
safe_save_plot <- function(filename, plot, width, height, dpi = 300, units = "in") {
  tryCatch({
    if (grepl("\\.png$", filename, ignore.case = TRUE) && requireNamespace("ragg", quietly = TRUE)) {
      ggsave(filename, plot, device = ragg::agg_png, width = width, height = height, dpi = dpi, units = units)
    } else {
      ggsave(filename, plot, width = width, height = height, dpi = dpi, units = units)
    }
  }, error = function(e) {
    warning("Fallback saver for ", filename, ": ", conditionMessage(e))
    ggsave(filename, plot, width = width, height = height, dpi = dpi, units = units)
  })
}
safe_pdf_begin <- function(file, width, height) grDevices::pdf(file, width = width, height = height)

# ---------- 1) Load metadata (batch1 + batch2 + batch3) ----------
batch1 <- read_csv("batch1.csv", show_col_types = FALSE)
batch2 <- read_csv("batch2.csv", show_col_types = FALSE)
batch3 <- read_csv("batch3.csv", show_col_types = FALSE)

# Drop any in-line S&S columns from the "metadata" frame (we’ll pull them separately and standardise)
drop_SnS_cols <- function(df) {
  df %>%
    select(
      -matches(
        paste(
          "S&S|sns.*status|signs.*symptoms.*status|sympt|dysur|urg|urge|fever|koorts|pijn|pain|burn",
          "|(^|[_ -])sx[_ -]?present($|[_ -])",
          sep = ""
        ),
        ignore.case = TRUE
      )
    )
}
b1_clean <- drop_SnS_cols(batch1)
b2_clean <- drop_SnS_cols(batch2)
b3_clean <- drop_SnS_cols(batch3)

common_cols <- Reduce(intersect, list(names(b1_clean), names(b2_clean), names(b3_clean)))
metadata <- bind_rows(
  mutate(b1_clean, across(all_of(common_cols), as.character), Batch = 1L),
  mutate(b2_clean, across(all_of(common_cols), as.character), Batch = 2L),
  mutate(b3_clean, across(all_of(common_cols), as.character), Batch = 3L)
)

# Standardize ID & Timepoint
standardize_id_tp <- function(df) {
  nm <- names(df)
  id_col <- nm[grepl("^participant[_ ]?id$", nm, ignore.case = TRUE)][1]
  tp_col <- nm[grepl("^time\\s*point$|^timepoint$", nm, ignore.case = TRUE)][1]
  if (is.na(id_col) || is.na(tp_col)) stop("Could not find Participant_id / Timepoint columns in metadata.")
  df %>%
    dplyr::rename(Participant_id = all_of(id_col), Timepoint = all_of(tp_col)) %>%
    dplyr::mutate(
      Participant_id = as.character(Participant_id),
      Timepoint      = canon_tp(Timepoint)
    )
}
metadata <- metadata %>% standardize_id_tp()
msg("✓ metadata loaded  (n = %d)", nrow(metadata))

# ---------- 2) FASTA discovery & metrics ----------
fasta_dir  <- "ont-yellow-routine-fastas"   # adjust if needed
all_fasta  <- list.files(fasta_dir, "\\.fasta$", full.names = TRUE)
if (length(all_fasta) == 0) stop("No FASTA files found in ", fasta_dir)

assembly_tbl <- tibble(full_path = all_fasta) %>%
  mutate(
    file_name  = basename(full_path),
    Isolate_ID = str_extract(str_split_fixed(file_name, "_", 4)[, 3], "[0-9A-Za-z]+-[0-9]+"),
    assembler  = case_when(
      str_detect(file_name, "flye")       ~ "flye",
      str_detect(file_name, "longcycler") ~ "longcycler",
      TRUE                                ~ "unknown"
    )
  )
summarise_fasta <- function(fp) {
  x  <- readDNAStringSet(fp)
  af <- colSums(alphabetFrequency(x, baseOnly = TRUE))
  tibble(
    num_contigs = length(x),
    total_bases = sum(width(x)),
    gc_content  = round((af["G"] + af["C"]) / sum(af) * 100, 2)
  )
}
metrics_tbl <- purrr::map_dfr(assembly_tbl$full_path, summarise_fasta)
assembly_df <- bind_cols(assembly_tbl %>% select(-full_path), metrics_tbl) %>%
  left_join(metadata, by = c("Isolate_ID" = "isolate_ID"), relationship = "many-to-many")

# ---------- 3) Recover per-symptom columns (B1–B3) + SnS_status (B1–B2 fallback only) ----------
# IMPORTANT: explicitly exclude any "*status*" field from the symptom matrix
symptom_cols2 <- function(df) {
  nm <- names(df)
  id_col <- nm[grepl("^participant[_ ]?id$", nm, ignore.case = TRUE)][1]
  tp_col <- nm[grepl("^time\\s*point$|^timepoint$", nm, ignore.case = TRUE)][1]
  if (is.na(id_col) || is.na(tp_col)) return(NULL)
  df %>%
    rename(Participant_id = all_of(id_col), Timepoint = all_of(tp_col)) %>%
    mutate(Participant_id = as.character(Participant_id),
           Timepoint      = canon_tp(Timepoint)) %>%
    select(
      Participant_id, Timepoint,
      matches("sympt|dysur|urg|urge|fever|koorts|pijn|pain|burn", ignore.case = TRUE),
      -matches("status|sns.*status|signs.*symptoms.*status", ignore.case = TRUE)
    ) %>%
    mutate(across(-c(Participant_id, Timepoint), as.character))
}
sym_b1 <- symptom_cols2(batch1)
sym_b2 <- symptom_cols2(batch2)
sym_b3 <- symptom_cols2(batch3)
sym_list    <- list(sym_b1, sym_b2, sym_b3)
sym_all_raw <- bind_rows(sym_list[!vapply(sym_list, is.null, logical(1))])

if (nrow(sym_all_raw) == 0) {
  sym_all_1row <- tibble(Participant_id = character(), Timepoint = character(), Sx_present = logical())
} else {
  sym_bool <- sym_all_raw %>%
    mutate(across(where(is.character), ~ str_to_lower(str_trim(.x)))) %>%
    mutate(
      across(
        -c(Participant_id, Timepoint),
        ~ {
          x <- .x
          num <- suppressWarnings(readr::parse_number(x))
          dplyr::case_when(
            !is.na(num) ~ num > 0,
            x %in% c("ja","yes","y","true","1","present","pos","positive") ~ TRUE,
            x %in% c("nee","no","n","false","0","","absent","neg","negative", NA_character_) ~ FALSE,
            TRUE ~ FALSE
          )
        }
      ),
      Timepoint = canon_tp(Timepoint)
    )
  sym_cols <- setdiff(names(sym_bool), c("Participant_id","Timepoint"))
  sym_all_1row <- sym_bool %>%
    group_by(Participant_id, Timepoint) %>%
    summarise(Sx_present = any(rowSums(across(all_of(sym_cols), ~ .x)) > 0), .groups = "drop") %>%
    distinct()
}
write_csv(sym_all_1row %>% arrange(Participant_id, Timepoint), file.path(dbg_dir, "sym_all.csv"))

# SnS_status only as BACKUP for B1–B2 if S&S columns missing for that PID×TP
extract_sns_status <- function(df, batch_id) {
  nm <- names(df)
  id_col <- nm[grepl("^participant[_ ]?id$", nm, ignore.case = TRUE)][1]
  tp_col <- nm[grepl("^time\\s*point$|^timepoint$", nm, ignore.case = TRUE)][1]
  if (is.na(id_col) || is.na(tp_col)) return(empty_sns_frame())
  df_std <- df %>%
    rename(Participant_id = all_of(id_col), Timepoint = all_of(tp_col)) %>%
    mutate(Participant_id = as.character(Participant_id),
           Timepoint      = canon_tp(Timepoint))
  status_col <- nm[grepl("s\\s*&\\s*s.*status|sns.*status|signs.*symptoms.*status", nm, ignore.case = TRUE)][1]
  if (is.na(status_col)) return(empty_sns_frame())
  df_std %>%
    transmute(
      Participant_id, Timepoint,
      SnS_status = suppressWarnings(as.integer(readr::parse_number(.data[[status_col]]))),
      has_SnS    = !is.na(SnS_status),
      Batch      = batch_id
    ) %>%
    filter(!is.na(SnS_status)) %>%
    group_by(Participant_id, Timepoint, Batch) %>%
    summarise(
      SnS_status = max(SnS_status, na.rm = TRUE),
      has_SnS    = any(has_SnS, na.rm = TRUE),
      .groups = "drop"
    )
}
sns_12 <- bind_rows(extract_sns_status(batch1, 1L),
                    extract_sns_status(batch2, 2L))
if (ncol(sns_12) == 0) sns_12 <- empty_sns_frame()
write_csv(sns_12 %>% arrange(Participant_id, Timepoint), file.path(dbg_dir, "sns_keys.csv"))

# ---------- 4) Classification with POPULATION-FIRST rules ----------
meta_plus <- metadata %>%
  left_join(sns_12, by = c("Participant_id","Timepoint","Batch"),
            relationship = "many-to-many") %>%
  left_join(sym_all_1row %>% dplyr::rename(Sx_present_sym = Sx_present),
            by = c("Participant_id","Timepoint"),
            relationship = "many-to-many")
if (!"Sx_present_sym" %in% names(meta_plus)) meta_plus$Sx_present_sym <- NA

classified <- meta_plus %>%
  mutate(
    Beoord_clean = if ("Beoord" %in% names(.)) str_to_lower(str_trim(Beoord)) else NA_character_,
    
    # CFU: explicit > / ≥ / "meer dan 100000" => positive; ranges not positive
    cfu_cat      = if ("CFU_Count" %in% names(.)) cfu_bucket(CFU_Count) else NA_character_,
    cfu_recorded = if ("CFU_Count" %in% names(.)) (!is.na(CFU_Count) & str_trim(as.character(CFU_Count)) != "") else FALSE,
    
    beoord_cat = case_when(
      str_detect(nz_chr(Beoord_clean), "\\+\\+\\+") ~ "+++",
      str_detect(nz_chr(Beoord_clean), "\\+\\+")    ~ "++",
      str_detect(nz_chr(Beoord_clean), "\\+")       ~ "+",
      TRUE ~ NA_character_
    ),
    
    # Culture positivity per rule
    culture_pos = case_when(
      cfu_recorded                        ~ (cfu_cat == ">=1e5"),
      !cfu_recorded & !is.na(beoord_cat) ~ (beoord_cat == "+++"),
      TRUE                               ~ NA
    ),
    
    # Population-derived S&S (authoritative when determinable, ALL batches)
    sx_present_pop = if ("Population" %in% names(.)) population_to_sns(Population) else as.logical(NA),
    
    # FINAL S&S decision (Population first)
    Sx_present_final = case_when(
      !is.na(sx_present_pop)                               ~ sx_present_pop,
      Batch %in% c(1L,2L) & !is.na(Sx_present_sym)         ~ Sx_present_sym,
      Batch %in% c(1L,2L) &  is.na(Sx_present_sym) &
        !is.na(SnS_status)                                  ~ case_when(
          SnS_status == 2L ~ TRUE,
          SnS_status == 0L ~ FALSE,
          TRUE             ~ as.logical(NA) # 1 or other -> unknown
        ),
      TRUE                                                  ~ as.logical(NA)
    ),
    
    # FINAL status
    Infection_Status = case_when(
      culture_pos == TRUE  &  Sx_present_final %in% TRUE   ~ "UTI",
      culture_pos == TRUE  &  Sx_present_final %in% FALSE  ~ "ASB",
      culture_pos == TRUE  &  is.na(Sx_present_final)      ~ "Culture-positive, S&S unknown",
      culture_pos == FALSE                                  ~ "Negative",
      TRUE                                                   ~ "None"
    )
  )

# ---------- 5) Episode-level collapse (unique PID×TP×Batch) ----------
episode_tbl <- classified %>%
  group_by(Participant_id, Timepoint, Batch) %>%
  summarise(
    culture_pos_epi   = any(culture_pos, na.rm = TRUE),
    cfu_recorded_any  = any(ifelse(is.na(cfu_recorded), FALSE, cfu_recorded)),
    cfu_ge_1e5_any    = any(cfu_cat == ">=1e5", na.rm = TRUE),
    beoord_plus3_any  = any(beoord_cat == "+++", na.rm = TRUE),
    
    # tri-state collapse for S&S at episode level
    Sx_present_any    = collapse_tristate(Sx_present_final),
    
    # Source labelling (Population first)
    from_pop          = any(!is.na(sx_present_pop)),
    from_cols         = any(!is.na(Sx_present_sym)),
    from_status       = any(!is.na(SnS_status)),
    
    Sx_source_epi     = case_when(
      from_pop    ~ "Population",
      from_cols   ~ "S&S columns",
      from_status ~ "SnS_status fallback",
      TRUE        ~ NA_character_
    ),
    
    raw_CFU_examples  = paste(head(unique(nz_chr(CFU_Count)), 5), collapse=" | "),
    raw_BEO_examples  = paste(head(unique(nz_chr(Beoord)),    5), collapse=" | "),
    raw_Pop_examples  = paste(head(unique(nz_chr(Population)),5), collapse=" | "),
    .groups = "drop"
  ) %>%
  select(-from_pop, -from_cols, -from_status) %>%
  mutate(
    Infection_Status = case_when(
      culture_pos_epi &  Sx_present_any %in% TRUE  ~ "UTI",
      culture_pos_epi &  Sx_present_any %in% FALSE ~ "ASB",
      culture_pos_epi &  is.na(Sx_present_any)     ~ "Culture-positive, S&S unknown",
      !culture_pos_epi                             ~ "Negative",
      TRUE                                         ~ "None"
    )
  )

# ---------- 6) Plotting + summaries (episode-level, excluding "None") ----------
tp_norm2 <- function(x) {
  tp_chr <- as.character(x)
  is_uricult <- str_detect(tp_chr, regex("uricult", ignore_case = TRUE))
  tp_num <- suppressWarnings(as.integer(str_extract(tp_chr, "\\d+")))
  tp_lab <- case_when(
    is_uricult ~ "Uricult",
    !is.na(tp_num) ~ paste0("T", tp_num),
    TRUE ~ "Unscheduled"
  )
  tibble(
    tp_num = tp_num,
    tp_lab = factor(tp_lab,
                    levels = c(paste0("T", sort(unique(tp_num[!is.na(tp_num)]))),
                               "Uricult","Unscheduled"))
  )
}
episode_tbl <- episode_tbl %>% bind_cols(tp_norm2(.$Timepoint))

episode_plot <- episode_tbl %>%
  filter(Infection_Status %in% c("Negative","ASB","UTI","Culture-positive, S&S unknown"))

status_order <- c("UTI","ASB","Culture-positive, S&S unknown","Negative")
status_levels_story <- c("Negative","ASB","UTI","Culture-positive, S&S unknown")

status_by_tp <- episode_plot %>%
  count(tp_lab, Infection_Status, name = "n") %>%
  group_by(tp_lab) %>% mutate(pct = 100*n/sum(n)) %>% ungroup() %>%
  arrange(tp_lab, desc(n))

trajectory_wide <- episode_plot %>%
  select(Participant_id, tp_lab, Infection_Status) %>%
  arrange(Participant_id, tp_lab) %>%
  tidyr::pivot_wider(names_from = tp_lab, values_from = Infection_Status, values_fn = list)

from_to <- episode_plot %>%
  filter(!is.na(tp_num)) %>%
  arrange(Participant_id, tp_num) %>%
  group_by(Participant_id) %>%
  mutate(Next_Status = dplyr::lead(Infection_Status)) %>%
  ungroup() %>%
  filter(!is.na(Next_Status)) %>%
  count(From = Infection_Status, To = Next_Status, name = "n") %>%
  complete(From = status_order, To = status_order, fill = list(n = 0)) %>%
  mutate(From = factor(From, levels = status_order),
         To   = factor(To,   levels = status_order))

p_status <- ggplot(status_by_tp, aes(tp_lab, n, fill = Infection_Status)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_stack(vjust = 0.5), size = 3) +
  labs(title = "Status distribution per timepoint (episode-level)", x = "Timepoint", y = "N episodes")

pid_levels <- episode_plot %>%
  summarise(first_tp = suppressWarnings(min(tp_num, na.rm = TRUE)), .by = Participant_id) %>%
  mutate(first_tp = ifelse(is.infinite(first_tp), Inf, first_tp)) %>%
  arrange(desc(first_tp), Participant_id) %>% pull(Participant_id)

plot_df <- episode_plot %>% mutate(Participant_id = factor(Participant_id, levels = pid_levels))
p_traj <- ggplot(plot_df, aes(tp_lab, Participant_id, fill = Infection_Status)) +
  geom_tile(color = "white") +
  scale_y_discrete(drop = FALSE) +
  labs(title = "Within-person infection status across time (episode-level)", x = "Timepoint", y = "Participant") +
  theme(axis.text.y = element_text(size = 6))

# ---- Robust transitions plot helper (fixes previous "Transitions plot error") ----
status_cols <- c(
  Negative = "#56B4E9",
  ASB      = "#E69F00",
  UTI      = "#D55E00",
  `Culture-positive, S&S unknown` = "#009E73"
)

make_transitions_plot <- function(from_to, status_levels_story, status_cols) {
  ft <- from_to %>%
    mutate(
      From = factor(as.character(From), levels = status_levels_story),
      To   = factor(as.character(To),   levels = status_levels_story)
    )
  
  if (!nrow(ft) || sum(ft$n, na.rm = TRUE) == 0) {
    return(
      ggplot() + theme_void() +
        ggtitle("No transitions between consecutive timepoints to display")
    )
  }
  
  if (requireNamespace("ggalluvial", quietly = TRUE)) {
    p_alluv <- try({
      ggplot(ft, aes(y = n, axis1 = From, axis2 = To)) +
        ggalluvial::geom_alluvium(aes(fill = From), width = 0, alpha = 0.9) +
        ggalluvial::geom_stratum(width = 0.15, fill = "grey85", colour = "grey40") +
        ggplot2::geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
        scale_fill_manual(values = status_cols, drop = FALSE) +
        labs(title = "Transitions between consecutive timepoints (episodes)",
             x = NULL, y = "Count", fill = "From status") +
        theme_minimal(base_size = 11)
    }, silent = TRUE)
    if (!inherits(p_alluv, "try-error")) return(p_alluv)
  }
  
  ggplot(ft, aes(From, To, fill = n)) +
    geom_tile() +
    geom_text(aes(label = n)) +
    scale_fill_continuous(name = "Count") +
    labs(title = "Transitions between consecutive timepoints (episodes, heatmap)",
         x = "From", y = "To") +
    theme_minimal(base_size = 11)
}
p_transitions_flow <- make_transitions_plot(from_to, status_levels_story, status_cols)

# Assembly metrics by status and timepoint (episode-level)
assembly_norm <- assembly_df %>%
  transmute(Participant_id, Timepoint, assembler, num_contigs, total_bases, gc_content) %>%
  bind_cols(tp_norm2(.$Timepoint))

assembly_by_status <- episode_plot %>%
  select(Participant_id, tp_lab, Infection_Status) %>%
  distinct() %>%
  inner_join(assembly_norm, by = c("Participant_id","tp_lab"),
             relationship = "many-to-many") %>%
  group_by(tp_lab, Infection_Status) %>%
  summarise(
    n_assemblies   = n(),
    mean_contigs   = round(mean(num_contigs, na.rm = TRUE), 2),
    median_contigs = median(num_contigs, na.rm = TRUE),
    mean_size_mb   = round(mean(total_bases, na.rm = TRUE)/1e6, 3),
    median_gc      = round(median(gc_content, na.rm = TRUE), 2),
    asm_mix        = paste(sort(table(assembler)), collapse = "|"),
    .groups = "drop"
  ) %>% arrange(tp_lab, Infection_Status)
write_csv(assembly_by_status, file.path(out_dir, "assembly_metrics_by_status.csv"))

p_asm_contigs <- ggplot(
  inner_join(episode_plot %>% select(Participant_id, tp_lab, Infection_Status) %>% distinct(),
             assembly_norm, by = c("Participant_id","tp_lab"), relationship = "many-to-many"),
  aes(Infection_Status, num_contigs)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  facet_wrap(~ tp_lab, scales = "free_x") +
  labs(title = "Assembly contig counts by status and timepoint (episode-level)",
       x = "Status", y = "# contigs")

# ---------- 7) STORYBOARD PLOTS ----------
base_tbl <- episode_plot %>%
  mutate(Infection_Status = factor(Infection_Status, levels = c("Negative","ASB","UTI","Culture-positive, S&S unknown"), exclude = NULL)) %>%
  select(Participant_id, tp_lab, tp_num, Infection_Status)

pid_levels2 <- base_tbl %>% distinct(Participant_id) %>% arrange(Participant_id) %>% pull(Participant_id)
tp_levels2  <- levels(base_tbl$tp_lab)

# A) Microtiles grid
micro <- base_tbl %>%
  mutate(
    pid_f = forcats::fct_relevel(Participant_id, pid_levels2),
    tp_f  = factor(tp_lab, levels = tp_levels2),
    x     = as.integer(pid_f),
    y     = as.integer(tp_f)
  ) %>%
  group_by(Participant_id, tp_lab) %>%
  arrange(Infection_Status, .by_group = TRUE) %>%
  mutate(
    N    = n(), k = row_number(),
    xmin = x - 0.5 + (k - 1) / N, xmax = x - 0.5 + (k) / N,
    ymin = y - 0.5, ymax = y + 0.5
  ) %>% ungroup()
micro_tot <- micro %>% count(Participant_id, tp_lab, x, y, name = "N")

p_micro <- ggplot() +
  geom_rect(data = micro,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Infection_Status),
            colour = NA) +
  geom_segment(data = micro,
               aes(x = xmin, xend = xmin, y = ymin, yend = ymax,
                   group = interaction(Participant_id, tp_lab)),
               inherit.aes = FALSE, linewidth = 0.15, colour = "white") +
  geom_rect(data = micro_tot,
            aes(xmin = x - 0.5, xmax = x + 0.5, ymin = y - 0.5, ymax = y + 0.5),
            fill = NA, colour = "grey50", linewidth = 0.25) +
  geom_text(data = dplyr::filter(micro_tot, N > 1),
            aes(x = x, y = y, label = N), size = 3, vjust = 0.5, colour = "black") +
  scale_fill_manual(values = status_cols, drop = FALSE) +
  scale_x_continuous(breaks = seq_along(pid_levels2), labels = pid_levels2, expand = expansion(add = 0)) +
  scale_y_continuous(breaks = seq_along(tp_levels2),  labels = tp_levels2,  expand = expansion(add = 0)) +
  labs(title = "Within-person infection status across time (episodes)",
       subtitle = "Each Participant×Timepoint cell is subdivided into one slice per episode. Numbers show total readings when >1.",
       x = "Participant", y = "Timepoint", fill = "Status") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# B) “Changers only”
changers <- base_tbl %>%
  distinct(Participant_id, tp_lab, Infection_Status) %>%
  group_by(Participant_id) %>%
  mutate(n_statuses_over_time = n_distinct(Infection_Status)) %>%
  ungroup() %>%
  filter(n_statuses_over_time > 1)
p_changers <- if (nrow(changers)) {
  ggplot(changers, aes(tp_lab, Participant_id, fill = Infection_Status)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = status_cols, drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(title = sprintf("Participants with status change over time (n = %d)", n_distinct(changers$Participant_id)),
         x = "Timepoint", y = "Participant", fill = "Status") +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 6))
} else {
  ggplot() + theme_void() + ggtitle("No participants changed status across time")
}

# ---------- 8) Save plots ----------
msg("✓ Saving plots to %s …", normalizePath(plot_dir))
try(ggsave(file.path(plot_dir, "status_by_timepoint.png"),      p_status,      width = 8,  height = 5,  dpi = 300), silent = TRUE)
try(ggsave(file.path(plot_dir, "trajectories_heatmap.png"),     p_traj,        width = 10, height = 12, dpi = 300), silent = TRUE)
try(ggsave(file.path(plot_dir, "transitions_alluvial_or_heatmap.png"), p_transitions_flow, width = 12, height = 8,  dpi = 300), silent = TRUE)
try(ggsave(file.path(plot_dir, "assembly_contigs_boxplot.png"), p_asm_contigs, width = 10, height = 6,  dpi = 300), silent = TRUE)
# Extra storytelling
# 1) Waterfall
waterfall <- tibble::tibble(
  step = c("All episodes",
           "Culture positive (strict CFU/Beoord)",
           "S&S determinable",
           "S&S present",
           "Final UTIs"),
  n = c(
    nrow(episode_tbl),
    sum(episode_tbl$culture_pos_epi, na.rm = TRUE),
    sum(!is.na(episode_tbl$Sx_present_any)),
    sum(episode_tbl$Sx_present_any %in% TRUE, na.rm = TRUE),
    sum(episode_tbl$Infection_Status == "UTI")
  )
)
readr::write_csv(waterfall, file.path(out_dir, "waterfall_counts.csv"))
p_waterfall <- ggplot(waterfall, aes(x = factor(step, levels = step), y = n)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.3, size = 3) +
  labs(title = "How rules narrow episodes down to UTIs",
       x = NULL, y = "Episodes") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
safe_save_plot(file.path(plot_dir, "waterfall_counts.png"), p_waterfall, width = 9, height = 5, dpi = 300)

# 2) Status by timepoint × batch
status_by_tp_batch <- episode_plot %>%
  count(Batch, tp_lab, Infection_Status, name = "n") %>%
  group_by(Batch, tp_lab) %>% mutate(pct = 100*n/sum(n)) %>% ungroup()
readr::write_csv(status_by_tp_batch, file.path(out_dir, "status_by_timepoint_batch.csv"))
p_status_by_tp_batch <- ggplot(status_by_tp_batch, aes(tp_lab, n, fill = Infection_Status)) +
  geom_col() +
  facet_wrap(~ Batch, nrow = 1, scales = "free_y") +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            position = position_stack(vjust = 0.5), size = 3) +
  labs(title = "Status distribution per timepoint, by batch",
       x = "Timepoint", y = "Episodes") +
  theme_minimal(base_size = 11)
safe_save_plot(file.path(plot_dir, "status_by_timepoint_by_batch.png"),
               p_status_by_tp_batch, width = 12, height = 4, dpi = 300)

# Storyboard exports
safe_save_plot(file.path(plot_dir, "story_microtiles.png"),       p_micro,            width = 16, height = 9,  dpi = 300)
safe_save_plot(file.path(plot_dir, "story_changers.png"),         p_changers,         width = 10, height = 10, dpi = 300)
safe_save_plot(file.path(plot_dir, "story_transitions_flow.png"), p_transitions_flow, width = 12, height = 8,  dpi = 300)

# PDFs
safe_pdf_begin(file.path(plot_dir, "overview_plots.pdf"), width = 10, height = 6)
print(p_status); print(p_traj); print(p_transitions_flow); print(p_asm_contigs)
dev.off()

safe_pdf_begin(file.path(plot_dir, "overview_status_story.pdf"), width = 16, height = 9)
print(p_micro); print(p_changers); print(p_transitions_flow); 
dev.off()

# ---------- 9) CSV exports / QA (episode-level, excluding "None") ----------
write_csv(status_by_tp,    file.path(out_dir, "status_by_timepoint.csv"))
write_csv(trajectory_wide, file.path(out_dir, "participant_trajectories_wide.csv"))

# UTI audits (episode-level)
uti_epi <- episode_tbl %>%
  filter(Infection_Status == "UTI") %>%
  mutate(culture_driver = case_when(
    cfu_recorded_any & cfu_ge_1e5_any    ~ "CFU > 100k (explicit)",
    !cfu_recorded_any & beoord_plus3_any ~ "Beoord +++ (CFU missing)",
    TRUE                                 ~ "Other/unclear"
  ))
write_csv(
  uti_epi %>% select(Participant_id, Timepoint, Batch, tp_lab, Infection_Status,
                     culture_pos_epi, culture_driver,
                     cfu_recorded_any, cfu_ge_1e5_any, beoord_plus3_any,
                     Sx_present_any, Sx_source_epi,
                     raw_CFU_examples, raw_BEO_examples, raw_Pop_examples),
  file.path(dbg_dir, "audit_uti_episode_review.csv")
)
write_csv(
  uti_epi %>% count(Batch, culture_driver, Sx_source_epi, name="n") %>% arrange(desc(n)),
  file.path(dbg_dir, "audit_uti_breakdown_by_source.csv")
)

# Extra QA: UTIs where Population says "no S&S" (should be 0 with Population-first)
uti_pop_conflict <- episode_tbl %>%
  inner_join(classified %>% distinct(Participant_id, Timepoint, Batch, sx_present_pop),
             by = c("Participant_id","Timepoint","Batch")) %>%
  filter(Infection_Status == "UTI", sx_present_pop == FALSE)
write_csv(uti_pop_conflict, file.path(dbg_dir, "audit_utis_population_conflict.csv"))

# S&S disagreement audit (B1–B2 only; where BOTH sources exist)
sns_disagree <- classified %>%
  filter(Batch %in% c(1L,2L)) %>%
  filter(!is.na(Sx_present_sym), !is.na(SnS_status)) %>%
  transmute(Participant_id, Timepoint, Batch,
            Sx_from_columns = Sx_present_sym,
            Sx_from_status  = case_when(
              SnS_status == 2L ~ TRUE,
              SnS_status == 0L ~ FALSE,
              TRUE             ~ as.logical(NA)
            ),
            disagree = Sx_from_columns != Sx_from_status) %>%
  arrange(desc(disagree), Participant_id, Timepoint)
write_csv(sns_disagree, file.path(dbg_dir, "audit_sns_disagreement.csv"))

# Hypothetical: flips if we trusted SnS_status over Population/columns (B1–B2 only)
flip_if_use_status <- episode_tbl %>%
  filter(Batch %in% c(1L,2L)) %>%
  left_join(classified %>% select(Participant_id, Timepoint, Batch, Sx_present_sym, SnS_status) %>% distinct(),
            by = c("Participant_id","Timepoint","Batch")) %>%
  mutate(Sx_from_status = case_when(
    SnS_status == 2L ~ TRUE,
    SnS_status == 0L ~ FALSE,
    TRUE             ~ as.logical(NA)
  ),
  Sx_from_cols   = Sx_present_sym) %>%
  filter(!is.na(Sx_from_status), !is.na(Sx_from_cols)) %>%
  transmute(Participant_id, Timepoint, Batch,
            Infection_Status_current = Infection_Status,
            Infection_Status_if_status = case_when(
              culture_pos_epi &  Sx_from_status %in% TRUE  ~ "UTI",
              culture_pos_epi &  Sx_from_status %in% FALSE ~ "ASB",
              culture_pos_epi &  is.na(Sx_from_status)     ~ "Culture-positive, S&S unknown",
              !culture_pos_epi                             ~ "Negative",
              TRUE                                         ~ "None"
            ),
            would_flip = Infection_Status_current != Infection_Status_if_status)
write_csv(flip_if_use_status, file.path(dbg_dir, "audit_hypothetical_use_status.csv"))

# Collapsed maps (episode-level)
status_map <- episode_plot %>%
  transmute(Participant_id, Timepoint = as.character(tp_lab), Infection_Status) %>% distinct()
write_csv(status_map, file.path(out_dir, "status_map_pid_tp_status.csv"))
write_csv(status_map, file.path(out_dir, "status_map.csv"))

status_map_full <- episode_plot %>%
  transmute(Participant_id, Timepoint = as.character(tp_lab), tp_lab, Infection_Status) %>%
  arrange(Participant_id, tp_lab)
write_csv(status_map_full, file.path(out_dir, "status_map_full.csv"))

# Additional QA tables you printed before
qa_utis_by_tp <- episode_tbl %>%
  filter(Infection_Status == "UTI") %>%
  mutate(culture_driver = case_when(
    cfu_ge_1e5_any    ~ "CFU>100k (explicit)",
    beoord_plus3_any  ~ "Beoord +++ (CFU missing)",
    TRUE              ~ "Other/unclear"
  )) %>%
  count(tp_lab, Sx_source_epi, culture_driver, name = "n") %>%
  arrange(tp_lab, desc(n))
readr::write_csv(qa_utis_by_tp, file.path(dbg_dir, "qa_utis_by_tp.csv"))

qa_cpos_by_tp <- episode_tbl %>%
  count(tp_lab, culture_pos_epi, Sx_present_any, name = "n") %>%
  arrange(tp_lab, desc(n))
readr::write_csv(qa_cpos_by_tp, file.path(dbg_dir, "qa_culturepos_by_tp.csv"))

# Population strings by TP (top 10 per TP)
tmp_tp <- tp_norm2(metadata$Timepoint)
pop_by_tp <- metadata %>%
  bind_cols(tmp_tp) %>%
  mutate(pop_l = tolower(nz_chr(Population))) %>%
  count(tp_lab, pop_l, sort = TRUE) %>%
  group_by(tp_lab) %>% slice_head(n = 10) %>% ungroup()
readr::write_csv(pop_by_tp, file.path(dbg_dir, "qa_population_strings_by_tp.csv"))

# Tiny README to sit next to outputs
readme_lines <- c(
  "Classification summary:",
  "- Culture POSITIVE only if CFU text explicitly says > or ≥ 100,000 (or Beoord +++ when CFU missing).",
  "- Ranges (e.g., 10k–100k; 100000–200000) are NOT positive.",
  "- S&S is Population-first (pt zonder UWI = FALSE; pt met UWI voor/na meting = FALSE at this TP; pt met UWI = TRUE).",
  "- If Population is missing/indeterminate: use explicit symptom columns; if missing (B1–B2), SnS_status fallback (2=TRUE, 0=FALSE, else unknown).",
  paste0("- Single numeric >100k alone counts as positive CFU? CFU_ACCEPT_SINGLE_NUMERIC_GT = ",
         CFU_ACCEPT_SINGLE_NUMERIC_GT),
  "",
  "Sanity notes:",
  "- With these data, scheduled TPs mostly say 'pt zonder UWI' or 'pt met UWI voor/na meting', so S&S=FALSE there;",
  "  therefore UTIs appear only at Uricult, where Population says 'pt met UWI' and culture is positive."
)
writeLines(readme_lines, file.path(out_dir, "README_rules.txt"))

msg("✓ Done. Plots in %s , CSVs in %s , debug CSVs in %s",
    normalizePath(plot_dir), normalizePath(out_dir), normalizePath(dbg_dir))

# Final sanity: list what's in debug right now
dbg_files <- list.files(dbg_dir, full.names = TRUE)
if (length(dbg_files)) {
  info <- file.info(dbg_files)
  out <- data.frame(file = normalizePath(rownames(info)),
                    size_kb = round(info$size/1024, 1),
                    mtime = info$mtime, row.names = NULL)
  print(out); msg("Listed %d files from debug.", nrow(out))
} else {
  msg("No files found in: %s", normalizePath(dbg_dir))
}
