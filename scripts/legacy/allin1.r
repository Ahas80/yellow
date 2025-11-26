#!/usr/bin/env Rscript
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
#!/usr/bin/env Rscript
# =============================================================
# 02_gene_presence_analysis.R
# -------------------------------------------------------------
#  * Expects `assembly_df` in the workspace  – OR reads the CSV
#  * Runs Abricate (VFDB) on every FASTA (cached to /results/abricate)
#  * Produces presence/absence matrices, prevalence tables & plots
#  * Optional nucmer trajectories for participants with ≥2 time-points
# =============================================================
## ---------------- 1 ·  libraries  -------------------------------------------
suppressPackageStartupMessages({
  library(dplyr);  library(tidyr);  library(readr)
  library(purrr);  library(furrr);  library(stringr)
  library(ggplot2); library(forcats); library(glue)
})

## ---------------- 2 ·  load assembly_df  ------------------------------------
if (!exists("assembly_df")) {
  if (file.exists("assembly_metadata.csv")) {
    assembly_df <- read_csv("assembly_metadata.csv", show_col_types = FALSE)
    message("✓  loaded assembly_df from assembly_metadata.csv  (",
            nrow(assembly_df), " rows)")
  } else {
    stop("assembly_df not found – run 01_prepare_assembly_metadata.R first")
  }
}
## ---- 2b · ensure full_path exists  -----------------------------
if (!"full_path" %in% names(assembly_df)) {
  if (!"file_name" %in% names(assembly_df))
    stop("assembly_df lacks both full_path and file_name columns.")
  
  fasta_dir <- "ont-yellow-routine-fastas"   # use repo-relative path
  if (!dir.exists(fasta_dir))
    stop("FASTA dir not found: ", fasta_dir, " (adjust in 02_gene_presence_analysis.R)")
  assembly_df <- assembly_df %>%
    mutate(full_path = file.path(fasta_dir, file_name),
           found     = file.exists(full_path))
  
  if (any(!assembly_df$found)) {
    bad <- assembly_df %>% filter(!found) %>% pull(full_path)
    stop("Missing FASTA files:\n", paste(bad, collapse = "\n"))
  }
  assembly_df <- select(assembly_df, -found)  # clean-up helper col
  message("✓  full_path column rebuilt from file_name (", nrow(assembly_df), " rows)")
}
## ---- 2c · normalize timepoints on assemblies (tp_lab + tp_num) --------------
# Drop-in replacement for tp_norm() in 02_gene_presence_analysis.R
# --- replace your current 2c block ---
## ---- 2c · normalize timepoints on assemblies (tp_lab + tp_num) --------------
tp_norm <- function(x) {
  x <- as.character(x)
  is_uricult <- stringr::str_detect(x, stringr::regex("uricult", ignore_case = TRUE))
  tp_num <- as.integer(stringr::str_extract(x, "\\d+"))
  tp_num[is_uricult] <- NA_integer_
  tp_lab <- dplyr::case_when(
    is_uricult ~ "Uricult",
    !is.na(tp_num) ~ paste0("T", tp_num),
    TRUE ~ "Unscheduled"
  )
  tp_lab <- factor(tp_lab,
                   levels = c(paste0("T", sort(unique(tp_num[!is.na(tp_num)]))),
                              "Uricult","Unscheduled"))
  tibble::tibble(tp_lab = tp_lab, tp_num = tp_num)
}

# Overwrite clean tp_* columns (avoid ... suffixes)
tp <- tp_norm(assembly_df$Timepoint)
assembly_df <- assembly_df %>%
  dplyr::mutate(tp_lab = tp$tp_lab, tp_num = tp$tp_num)

# Drop any leftover suffixed duplicates from previous runs
dup_tp <- grep("^tp_(lab|num)\\.\\.\\.[0-9]+$", names(assembly_df), value = TRUE)
if (length(dup_tp)) {
  assembly_df <- dplyr::select(assembly_df, -tidyselect::all_of(dup_tp))
}

# (Do NOT touch vf_pa_all here; it's not built yet)



# Drop ANY pre-existing tp_* columns (including ...##) before adding fresh ones


## ---------------- 3 ·  global output dirs  ----------------------------------
dir.create("results/abricate",   recursive = TRUE, showWarnings = FALSE)
dir.create("results/plots",      recursive = TRUE, showWarnings = FALSE)

## ---------------- 4 ·  Abricate-VFDB scan  -----------------------------------
run_abr_cached <- function(fasta, db = "vfdb", min_cov = 70, min_id = 70) {
  cache <- file.path("results/abricate", paste0(basename(fasta), ".vfdb.tsv"))
  if (file.exists(cache))
    return(readr::read_tsv(cache, show_col_types = FALSE, progress = FALSE))
  
  cmd <- glue::glue("abricate --quiet --db {db} --mincov {min_cov} --minid {min_id} {shQuote(fasta)} > {shQuote(cache)}")
  exit <- system(cmd)
  if (exit != 0) warning("Abricate non-zero exit (", basename(fasta), ")")
  readr::read_tsv(cache, show_col_types = FALSE, progress = FALSE)
}

future::plan(future::multisession, workers = max(1, parallel::detectCores() - 1))
on.exit(future::plan(sequential), add = TRUE)
safe_abr <- purrr::safely(run_abr_cached, otherwise = NULL, quiet = TRUE)

vf_hits_all <- assembly_df %>%
  dplyr::mutate(vfdb = furrr::future_map(full_path, ~ safe_abr(.x)$result, .progress = TRUE)) %>%
  dplyr::filter(purrr::map_lgl(vfdb, ~ !is.null(.x) && NROW(.x) > 0)) %>%
  tidyr::unnest(vfdb)

if (!"tp_lab" %in% names(vf_hits_all)) {
  tp_col <- grep("^tp_lab", names(vf_hits_all), value = TRUE)[1]
  if (length(tp_col)) vf_hits_all <- dplyr::rename(vf_hits_all, tp_lab = !!tp_col)
}

# Be resilient to different gene column names from abricate
gene_col <- intersect(c("GENE","GENE_NAME","NAME","PRODUCT","GENE SYMBOL"), names(vf_hits_all))[1]
if (is.na(gene_col)) stop("No gene name column found in Abricate output.")
vf_hits_all <- vf_hits_all %>% dplyr::rename(GENE = dplyr::all_of(gene_col))

saveRDS(vf_hits_all, "results/vf_hits_all.rds")
message("✓  vf_hits_all saved  (", nrow(vf_hits_all), " rows)")

## ---------------- 5 ·  presence/absence matrix ------------------------------
vf_pa_all <- vf_hits_all %>%
  dplyr::distinct(Participant_id, tp_lab, GENE) %>%   # keep many-to-many at sample level
  dplyr::mutate(present = 1) %>%
  tidyr::pivot_wider(names_from = GENE, values_from = present, values_fill = 0)

readr::write_csv(vf_pa_all, "results/vf_pa_all.csv")

## ---------------- 6 ·  gene-level prevalence --------------------------------
tbl_gene <- vf_hits_all %>%
  dplyr::distinct(Participant_id, GENE) %>%
  dplyr::count(GENE, name = "n_participants") %>%
  dplyr::arrange(dplyr::desc(n_participants))

readr::write_csv(tbl_gene, "results/stats_gene_level.csv")

## ---------------- 7 ·  quick cohort plots -----------------------------------
dir.create("results/plots", showWarnings = FALSE)

top25 <- tbl_gene %>% slice_max(n_participants, n = 25) %>%
  mutate(GENE = fct_reorder(GENE, n_participants))

ggplot(top25, aes(GENE, n_participants)) +
  geom_col(fill = "steelblue") + coord_flip() +
  labs(title = "Top 25 VFDB genes (entire cohort)",
       y = "Participants", x = NULL) +
  theme_minimal(base_size = 11)
ggsave("results/plots/core_bar_top25_all.png", width = 6, height = 6)

ggplot(tbl_gene, aes(n_participants)) +
  geom_histogram(binwidth = 1, fill = "grey70") +
  labs(title = "VF gene prevalence distribution",
       x = "# participants", y = "Gene count") +
  theme_minimal(base_size = 11)
ggsave("results/plots/core_histogram_all.png", width = 5, height = 4)

## ---------------- 8 ·  per-participant Heat-map & UpSet ---------------------
library(ComplexUpset)

variable_genes <- tbl_gene %>%
  dplyr::filter(dplyr::between(n_participants, 1, max(n_participants) - 1)) %>%
  dplyr::pull(GENE)

vf_pa_var <- vf_pa_all %>%
  dplyr::select(Participant_id, tp_lab, tidyselect::all_of(variable_genes))

plot_heatmap <- function(pid) {
  vf_pa_var %>%
    dplyr::filter(Participant_id == pid) %>%
    tidyr::pivot_longer(-c(Participant_id, tp_lab),
                        names_to = "GENE", values_to = "present") %>%
    ggplot2::ggplot(ggplot2::aes(GENE, tp_lab, fill = factor(present))) +
    ggplot2::geom_tile(color = "grey80") +
    ggplot2::scale_fill_manual(values = c(`0` = "white", `1` = "steelblue")) +
    ggplot2::labs(title = paste("Variable genes – participant", pid),
                  x = NULL, y = "Time-point", fill = "Present") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   legend.position = "none")
}

get_upset_df <- function(pid) {
  vf_pa_all %>%
    dplyr::filter(Participant_id == pid) %>%
    tidyr::pivot_longer(-c(Participant_id, tp_lab),
                        names_to = "GENE", values_to = "present") %>%
    dplyr::filter(present == 1) %>%
    tidyr::unite(GENE_TP, GENE, tp_lab, sep = "_") %>%
    dplyr::mutate(val = TRUE) %>%
    tidyr::pivot_wider(names_from = GENE_TP, values_from = val, values_fill = FALSE)
}

dir.create("results/plots/participants", showWarnings = FALSE)
purrr::walk(unique(vf_pa_all$Participant_id), \(pid) {
  ggplot2::ggsave(glue::glue("results/plots/participants/heatmap_{pid}.png"),
                  plot_heatmap(pid), width = 12, height = 6, dpi = 300)
  
  df <- get_upset_df(pid)
  ComplexUpset::upset(df,
                      intersect = setdiff(names(df), "Participant_id"),
                      min_size = 1, name = as.character(pid))
  ggplot2::ggsave(glue::glue("results/plots/participants/upset_{pid}.png"),
                  width = 12, height = 8, dpi = 300)
})

## ---------------- 9 ·  nucmer trajectories (≥2 numeric TPs) ------------------
ids_multi <- vf_pa_all %>%
  dplyr::distinct(Participant_id, tp_lab) %>%
  dplyr::mutate(tp_num = readr::parse_number(as.character(tp_lab))) %>%
  dplyr::filter(!is.na(tp_num)) %>%
  dplyr::count(Participant_id, name = "n_tp") %>%
  dplyr::filter(n_tp >= 2) %>%
  dplyr::pull(Participant_id)

if (length(ids_multi)) {
  have_mummer <- all(nzchar(Sys.which(c("nucmer", "delta-filter", "dnadiff"))))
  if (!have_mummer) {
    message("↪  nucmer/dnadiff not found on PATH; skipping trajectories.")
  } else {
    message("→  nucmer on ", length(ids_multi), " multi-TP participants")
    assembly_long <- assembly_df %>%
      dplyr::filter(Participant_id %in% ids_multi) %>%
      dplyr::filter(!is.na(tp_num))   # numeric T0/T1/... only
    
    pair_tbl <- assembly_long %>%
      dplyr::group_by(Participant_id, assembler) %>%
      dplyr::arrange(tp_num, .by_group = TRUE) %>%
      dplyr::mutate(path_A = full_path,
                    path_B = dplyr::lead(full_path),
                    tp_A   = tp_lab,
                    tp_B   = dplyr::lead(tp_lab),
                    outdir = glue::glue("results/nucmer/{Participant_id}_{assembler}_{tp_A}_vs_{tp_B}")) %>%
      dplyr::filter(!is.na(path_B)) %>%
      dplyr::ungroup()
    
    run_nucmer <- function(a, b, od) {
      dir.create(od, recursive = TRUE, showWarnings = FALSE)
      pref <- file.path(od, "run")
      system(glue::glue("nucmer --mum --prefix={shQuote(pref)} {shQuote(a)} {shQuote(b)}"))
      system(glue::glue("delta-filter -1 {pref}.delta > {pref}.1delta"))
      system(glue::glue("dnadiff -p {pref}_dd -d {pref}.1delta"))
      rpt <- glue::glue("{pref}_dd.report")
      if (!file.exists(rpt)) return(tibble::tibble(AvgIdentity = NA_real_, TotalSnpCnt = NA_real_))
      L <- readLines(rpt)
      grab <- \(k) { m <- grep(k, L, value = TRUE); if (length(m)) as.numeric(stringr::str_extract(m[1], "\\d+\\.?\\d*")) else NA_real_ }
      tibble::tibble(AvgIdentity = grab("AvgIdentity"), TotalSnpCnt = grab("TotalSNPs"))
    }
    
    pair_df <- pair_tbl %>%
      dplyr::mutate(res = furrr::future_pmap(list(path_A, path_B, outdir), run_nucmer, .progress = TRUE)) %>%
      tidyr::unnest(res)
    
    dir.create("results/plots/pairwise_identity", showWarnings = FALSE)
    
    plot_pair <- function(pid, asm) {
      dat <- pair_df %>%
        dplyr::filter(Participant_id == pid, assembler == asm) %>%
        dplyr::arrange(tp_num)                 # keep temporal order
      
      if (!nrow(dat)) return(NULL)
      
      g <- ggplot2::ggplot(dat, ggplot2::aes(x = tp_A, y = AvgIdentity, colour = TotalSnpCnt))
      
      if (nrow(dat) > 1) {
        g <- g + ggplot2::geom_line(ggplot2::aes(group = 1), linewidth = 0.8)
      }
      
      g + ggplot2::geom_point(size = 3) +
        ggplot2::scale_colour_viridis_c(option = "D", name = "SNP count") +
        ggplot2::labs(title = glue::glue("Pairwise identity – P{pid} ({asm})"),
                      x = "Time-point A", y = "Avg identity (%)") +
        ggplot2::theme_minimal(base_size = 10)
    }
    
    pw <- pair_df %>% dplyr::distinct(Participant_id, assembler)
    pw %>% purrr::pwalk(\(Participant_id, assembler) {
      g <- plot_pair(Participant_id, assembler)
      if (!is.null(g))
        ggplot2::ggsave(glue::glue("results/plots/pairwise_identity/P{Participant_id}_{assembler}.png"),
                        g, width = 6, height = 4, dpi = 300)
    })
  }
}
#!/usr/bin/env Rscript
# =============================================================
# 03_plotting.R (refactored)
# -------------------------------------------------------------
# Reads:
#   results/vf_hits_all.rds
#   results/vf_pa_all.csv
#   results/stats_gene_level.csv
#   results/status_map.csv                [optional]
#   results/diff_genes_UTI_vs_ASB_fisher.csv   [optional]
#
# Writes:
#   results/stats_*csv
#   results/stats_descriptive.md          [if knitr available]
#   results/plots/*.png (and subfolders)
#   results/persistence_*.csv
#   results/permanova_UTI_vs_ASB.txt      [if created earlier]
#
# Adds richer plots:
#   - Volcano (UTI vs ASB enrichment)
#   - Status-stratified prevalence heatmap
#   - PCoA (Jaccard) by status
#   - Gene–gene co-occurrence heatmap
#   - Persistence vs transient bars per participant
#   - UpSet improvements:
#       * Persistence per participant (sets = timepoints)
#       * Status-stratified (sets = UTI/ASB/Negative/None)
#       * Family-level (sets = gene families across participants)
# =============================================================

## ---------------- 0 · setup & helpers ---------------------------------------
suppressPackageStartupMessages({
  library(dplyr);    library(tidyr);    library(readr)
  library(ggplot2);  library(forcats);  library(glue)
  library(purrr);    library(stringr)
})

has_pkg <- function(p) { requireNamespace(p, quietly = TRUE) }
use_pkg <- function(p) { if (has_pkg(p)) suppressPackageStartupMessages(library(p, character.only = TRUE)) }

plots_dir <- "results/plots"
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_minimal(base_size = 11))
theme_update(legend.title = element_text(hjust = 0),
             legend.text  = element_text(hjust = 0))

tp_norm <- function(x) {
  tp_chr     <- as.character(x)
  is_uricult <- str_detect(tp_chr, regex("uricult", ignore_case = TRUE))
  tp_num     <- suppressWarnings(as.integer(str_extract(tp_chr, "\\d+")))
  tp_lab     <- dplyr::case_when(
    is_uricult     ~ "Uricult",
    !is.na(tp_num) ~ paste0("T", tp_num),
    TRUE           ~ "Unscheduled"
  )
  tp_levels <- c(paste0("T", sort(unique(tp_num[!is.na(tp_num)]))),
                 "Uricult", "Unscheduled")
  tibble::tibble(tp_lab = factor(tp_lab, levels = tp_levels),
                 tp_num = tp_num)
}

ensure_tp_lab <- function(df) {
  if ("tp_lab" %in% names(df)) return(df)
  if (!"Timepoint" %in% names(df)) stop("Need 'tp_lab' or 'Timepoint' in data frame.")
  bind_cols(df, tp_norm(df$Timepoint))
}

safe_ggsave <- function(path, plot = last_plot(), width = 7, height = 5, dpi = 300) {
  try(ggsave(path, plot = plot, width = width, height = height, dpi = dpi), silent = TRUE)
}

## ---------------- 1 · load core data ----------------------------------------
vf_hits_path <- "results/vf_hits_all.rds"
vf_pa_path   <- "results/vf_pa_all.csv"
tbl_gene_path<- "results/stats_gene_level.csv"

if (!file.exists(vf_hits_path) || !file.exists(vf_pa_path) || !file.exists(tbl_gene_path)) {
  stop("Missing required inputs. Make sure 02_gene_presence_analysis.R has run.\n",
       "- results/vf_hits_all.rds\n- results/vf_pa_all.csv\n- results/stats_gene_level.csv")
}

vf_hits_all <- readRDS(vf_hits_path)
vf_pa_all   <- read_csv(vf_pa_path, show_col_types = FALSE)
tbl_gene    <- read_csv(tbl_gene_path, show_col_types = FALSE)

vf_pa_all   <- ensure_tp_lab(vf_pa_all)

## ---------------- 2 · sample/participant/cheatsheet -------------------------
# DISTINCT genes per sample (Participant_id, Timepoint)
genes_per_sample <- vf_hits_all %>%
  distinct(Participant_id, Timepoint, GENE) %>%
  count(Participant_id, Timepoint, name = "n_genes") %>%
  bind_cols(tp_norm(.$Timepoint))

participant_tbl <- genes_per_sample %>%
  filter(!is.na(tp_num)) %>%
  arrange(Participant_id, tp_num) %>%
  group_by(Participant_id) %>%
  summarise(
    n_timepoints = n(),
    mean_genes   = mean(n_genes),
    sd_genes     = sd(n_genes),
    min_genes    = min(n_genes),
    max_genes    = max(n_genes),
    genes_T0     = dplyr::first(n_genes),
    genes_last   = dplyr::last(n_genes),
    delta        = genes_last - genes_T0,
    .groups      = "drop"
  )

cohort_tbl <- participant_tbl %>%
  summarise(
    participants       = n(),
    timepoints_total   = sum(n_timepoints),
    genes_median_T0    = median(genes_T0),
    genes_median_last  = median(genes_last),
    mean_delta         = mean(delta),
    sd_delta           = sd(delta)
  )

write_csv(genes_per_sample, "results/stats_sample_level.csv")
write_csv(participant_tbl,  "results/stats_participant_level.csv")
write_csv(cohort_tbl,       "results/stats_cohort_level.csv")
message("✓  CSV summaries written to results/")

# Markdown cheat-sheet (if knitr is available)
if (has_pkg("knitr")) {
  md <- c(
    "# Descriptive statistics\n\n",
    "*Generated on:* ", as.character(Sys.Date()), "\n\n",
    "## Cohort overview\n",
    knitr::kable(cohort_tbl, format = "markdown"), "\n\n",
    "## Participant-level summary\n",
    knitr::kable(participant_tbl, format = "markdown"), "\n\n",
    "## First 10 samples\n",
    knitr::kable(head(genes_per_sample, 10), format = "markdown")
  )
  writeLines(md, "results/stats_descriptive.md")
}

## ---------------- 3 · core/accessory quick plots ----------------------------
dir.create(plots_dir, showWarnings = FALSE)

# Top-25 genes
top25 <- tbl_gene %>% slice_max(n_participants, n = 25) %>%
  mutate(GENE = fct_reorder(GENE, n_participants))
ggplot(top25, aes(GENE, n_participants)) +
  geom_col(fill = "steelblue") + coord_flip() +
  labs(title = "Top 25 VFDB genes (entire cohort)", y = "Participants", x = NULL)
safe_ggsave(file.path(plots_dir, "core_bar_top25_all.png"), width = 6, height = 6)

# Prevalence histogram
ggplot(tbl_gene, aes(n_participants)) +
  geom_histogram(binwidth = 1, fill = "grey70") +
  labs(title = "VF gene prevalence distribution", x = "# participants", y = "Gene count")
safe_ggsave(file.path(plots_dir, "core_histogram_all.png"), width = 5, height = 4)

# Richness by timepoint
ggplot(genes_per_sample, aes(tp_lab, n_genes)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  labs(title = "Per-sample VF gene richness by timepoint", x = "Timepoint", y = "# VF genes (distinct)")
safe_ggsave(file.path(plots_dir, "richness_by_timepoint.png"), width = 7, height = 4.5, dpi = 300)

# Trajectories (numeric)
traj_df <- genes_per_sample %>% filter(!is.na(tp_num))
if (nrow(traj_df) > 0) {
  ggplot(traj_df, aes(tp_num, n_genes, group = Participant_id)) +
    geom_line(alpha = 0.35) +
    geom_point(size = 1.6) +
    scale_x_continuous(breaks = sort(unique(traj_df$tp_num))) +
    labs(title = "Within-participant trajectories of VF gene richness", x = "Numeric timepoint", y = "# VF genes (distinct)")
  safe_ggsave(file.path(plots_dir, "richness_trajectories_numeric.png"), width = 7, height = 4.5, dpi = 300)
} else {
  message("↪  No numeric timepoints found; skipping trajectory plot.")
}

## ---------------- 4 · UpSet improvements ------------------------------------
if (has_pkg("ComplexUpset")) {
  use_pkg("ComplexUpset")
  # 4A) Persistence per participant (sets = timepoints; elements = genes)
  dir.create(file.path(plots_dir, "persistence"), showWarnings = FALSE)
  
  persist_upset_for_pid <- function(pid) {
    df <- vf_hits_all %>%
      filter(Participant_id == pid) %>%
      mutate(tp_lab = tp_norm(Timepoint)$tp_lab) %>%   # compute tp_lab once
      distinct(GENE, tp_lab) %>%
      mutate(val = TRUE) %>%
      pivot_wider(names_from = tp_lab, values_from = val, values_fill = FALSE) %>%
      distinct(GENE, .keep_all = TRUE)
    
    tp_sets <- setdiff(names(df), "GENE")
    if (length(tp_sets) < 2) return(NULL)
    
    ComplexUpset::upset(
      df, intersect = tp_sets, min_size = 1,
      name = paste("P", pid, "genes"),
      base_annotations = list('Intersection size' = intersection_size(text = list(size = 3)))
    )
  }
  
  persist_stats <- purrr::map_dfr(unique(vf_hits_all$Participant_id), function(pid) {
    df <- vf_hits_all %>%
      filter(Participant_id == pid) %>%
      mutate(tp_lab = tp_norm(Timepoint)$tp_lab) %>%
      distinct(GENE, tp_lab) %>%
      mutate(val = TRUE) %>%
      pivot_wider(names_from = tp_lab, values_from = val, values_fill = FALSE) %>%
      distinct(GENE, .keep_all = TRUE)
    
    tp_sets <- setdiff(names(df), "GENE")
    if (length(tp_sets) < 2) return(NULL)
    
    # ✅ across() inside mutate + select(...) is valid; rowSums() on logicals works
    df <- df %>%
      mutate(
        present_n = rowSums(select(., dplyr::all_of(tp_sets))),
        persist   = present_n == length(tp_sets)
      )
    
    readr::write_csv(df %>% arrange(desc(persist), desc(present_n)),
                     glue::glue("results/persistence_P{pid}.csv"))
    
    p <- persist_upset_for_pid(pid)
    if (!is.null(p)) safe_ggsave(glue::glue("{plots_dir}/persistence/upset_P{pid}.png"),
                                 p, width = 10, height = 6, dpi = 300)
    
    tibble::tibble(
      Participant_id   = pid,
      n_timepoints     = length(tp_sets),
      persistent_genes = sum(df$persist),
      transient_genes  = sum(!df$persist)
    )
  })
  
  if (nrow(persist_stats)) readr::write_csv(persist_stats, "results/persistence_summary.csv")
  
  # 4B) Status-stratified UpSet (sets = Infection_Status)
  if (file.exists("results/status_map.csv")) {
    status_map <- read_csv("results/status_map.csv", show_col_types = FALSE) %>%
      ensure_tp_lab() %>%
      select(Participant_id, tp_lab, Infection_Status) %>%
      distinct()
    gene_status_df <- vf_pa_all %>%
      pivot_longer(-c(Participant_id, tp_lab), names_to = "GENE", values_to = "present") %>%
      filter(present > 0) %>% select(-present) %>%
      left_join(status_map, by = c("Participant_id","tp_lab")) %>%
      filter(!is.na(Infection_Status)) %>%
      distinct(GENE, Infection_Status) %>%
      mutate(val = TRUE) %>%
      pivot_wider(names_from = Infection_Status, values_from = val, values_fill = FALSE)
    
    # optional prevalence trimming (5–95%)
    prev <- tbl_gene %>%
      mutate(p = n_participants / n_distinct(vf_pa_all$Participant_id)) %>%
      filter(p >= 0.05, p <= 0.95) %>% pull(GENE)
    gene_status_df <- gene_status_df %>% filter(GENE %in% prev)
    status_sets <- setdiff(names(gene_status_df), "GENE")
    if (length(status_sets) >= 2) {
      p <- ComplexUpset::upset(
        gene_status_df, intersect = status_sets, min_size = 1,
        name = "Genes by status",
        base_annotations = list('Intersection size' = intersection_size(text = list(size = 3)))
      )
      safe_ggsave(file.path(plots_dir, "upset_genes_by_status.png"), p, width = 10, height = 6, dpi = 300)
    }
    
    # Optional: limit to significant UTI vs ASB genes
    if (file.exists("results/diff_genes_UTI_vs_ASB_fisher.csv")) {
      sig <- read_csv("results/diff_genes_UTI_vs_ASB_fisher.csv", show_col_types = FALSE) %>%
        filter(p_adj <= 0.05) %>% pull(GENE)
      sig_df <- gene_status_df %>% filter(GENE %in% sig)
      if (nrow(sig_df)) {
        p <- ComplexUpset::upset(
          sig_df, intersect = status_sets, min_size = 1,
          name = "Sig. genes (UTI vs ASB)",
          base_annotations = list('Intersection size' = intersection_size(text = list(size = 3)))
        )
        safe_ggsave(file.path(plots_dir, "upset_sig_genes_by_status.png"), p, width = 10, height = 6, dpi = 300)
      }
    }
  }
  
  # 4C) Family-level UpSet (sets = gene families; elements = participants)
  get_family <- function(g) {
    m <- stringr::str_match(g, "^([A-Za-z]+)")[,2]
    ifelse(is.na(m), g, m)
  }
  fam_df <- vf_hits_all %>%
    mutate(FAMILY = get_family(GENE)) %>%
    distinct(Participant_id, FAMILY)
  fam_wide <- fam_df %>%
    mutate(val = TRUE) %>%
    pivot_wider(names_from = FAMILY, values_from = val, values_fill = FALSE) %>%
    distinct(Participant_id, .keep_all = TRUE)
  if (nrow(fam_wide)) {
    fam_counts <- colSums(fam_wide[, -1, drop = FALSE])
    n_p <- nrow(fam_wide)
    keep <- names(fam_counts)[fam_counts >= 0.05*n_p & fam_counts <= 0.95*n_p]
    fam_mat <- select(fam_wide, Participant_id, all_of(keep))
    if (ncol(fam_mat) > 2) {
      p <- ComplexUpset::upset(
        fam_mat, intersect = names(fam_mat)[-1], min_size = 1,
        name = "Families across participants",
        base_annotations = list('Intersection size' = intersection_size(text = list(size = 3)))
      )
      safe_ggsave(file.path(plots_dir, "upset_families_by_participant.png"), p, width = 11, height = 6, dpi = 300)
    }
  }
} else {
  message("↪  ComplexUpset not installed; skipping UpSet figures.")
}

## ---------------- 5 · Volcano: UTI vs ASB enrichment -------------------------
if (file.exists("results/diff_genes_UTI_vs_ASB_fisher.csv")) {
  if (has_pkg("ggrepel")) use_pkg("ggrepel")
  de <- read_csv("results/diff_genes_UTI_vs_ASB_fisher.csv", show_col_types = FALSE) %>%
    mutate(
      log2OR = log2(OR),
      neglog10FDR = -log10(p_adj),
      sig = p_adj <= 0.05,
      total_pos = coalesce(`UTI_TRUE`, 0) + coalesce(`ASB_TRUE`, 0)
    )
  top_labs <- bind_rows(
    de %>% arrange(p_adj) %>% slice_head(n = 20),
    de %>% arrange(desc(abs(log2OR))) %>% slice_head(n = 20)
  ) %>% distinct(GENE)
  
  g <- ggplot(de, aes(log2OR, neglog10FDR)) +
    geom_point(aes(alpha = sig, size = total_pos)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    { if (has_pkg("ggrepel")) ggrepel::geom_text_repel(data = semi_join(de, top_labs, by = "GENE"),
                                                       aes(label = GENE), max.overlaps = 40, size = 3) else NULL } +
    labs(title = "UTI vs ASB: per-gene enrichment",
         x = "log2(odds ratio)  (UTI / ASB)",
         y = "-log10(FDR)")
  safe_ggsave(file.path(plots_dir, "volcano_UTI_vs_ASB.png"), g, width = 7, height = 5, dpi = 300)
} else {
  message("↪  diff_genes_UTI_vs_ASB_fisher.csv not found; skipping volcano.")
}

## ---------------- 6 · Status-stratified prevalence heatmap -------------------
if (file.exists("results/status_map.csv")) {
  status_map <- read_csv("results/status_map.csv", show_col_types = FALSE)
  status_map <- status_map %>%
    mutate(tp_lab = if (!"tp_lab" %in% names(.)) tp_norm(Timepoint)$tp_lab else tp_lab) %>%
    select(Participant_id, tp_lab, Infection_Status) %>% distinct()
  
  prev_long <- vf_pa_all %>%
    pivot_longer(-c(Participant_id, tp_lab), names_to = "GENE", values_to = "present") %>%
    left_join(status_map, by = c("Participant_id","tp_lab")) %>%
    filter(!is.na(Infection_Status)) %>%
    group_by(GENE, Infection_Status) %>%
    summarise(prev = mean(present > 0) * 100, .groups = "drop")
  
  rank40 <- prev_long %>%
    pivot_wider(names_from = Infection_Status, values_from = prev) %>%
    mutate(delta_UTI_ASB = coalesce(UTI, 0) - coalesce(ASB, 0)) %>%
    arrange(desc(abs(delta_UTI_ASB))) %>%
    slice_head(n = 40) %>% pull(GENE)
  
  heat_df <- prev_long %>%
    filter(GENE %in% rank40) %>%
    mutate(
      Infection_Status = factor(Infection_Status, levels = c("UTI","ASB","Negative","None")),
      GENE = factor(GENE, levels = rev(rank40))
    )
  
  g <- ggplot(heat_df, aes(Infection_Status, GENE, fill = prev)) +
    geom_tile() +
    scale_fill_gradient(name = "% present", low = "white", high = "steelblue") +
    labs(title = "Gene prevalence by clinical status (top 40 by |UTI–ASB|)") +
    theme(axis.text.y = element_text(size = 6))
  safe_ggsave(file.path(plots_dir, "heatmap_prevalence_by_status_top40.png"), g, width = 7.5, height = 10, dpi = 300)
} else {
  message("↪  status_map.csv not found; skipping status heatmap.")
}

## ---------------- 7 · PCoA (Jaccard) by status -------------------------------
if (has_pkg("vegan")) {
  use_pkg("vegan"); if (has_pkg("ggrepel")) use_pkg("ggrepel")
  if (file.exists("results/status_map.csv")) {
    status_map <- read_csv("results/status_map.csv", show_col_types = FALSE) %>%
      mutate(tp_lab = if (!"tp_lab" %in% names(.)) tp_norm(Timepoint)$tp_lab else tp_lab) %>%
      select(Participant_id, tp_lab, Infection_Status) %>% distinct()
    ord_data <- vf_pa_all %>% left_join(status_map, by = c("Participant_id","tp_lab"))
  } else {
    ord_data <- vf_pa_all %>% mutate(Infection_Status = NA_character_)
  }
  X <- ord_data %>% select(-Participant_id, -tp_lab, -Infection_Status) %>% as.matrix()
  X <- X > 0
  if (nrow(X) >= 3) {
    dJ <- vegan::vegdist(X, method = "jaccard", binary = TRUE)
    pcoa <- stats::cmdscale(dJ, k = 2, eig = TRUE)
    ord <- tibble::tibble(PC1 = pcoa$points[,1], PC2 = pcoa$points[,2]) %>%
      bind_cols(ord_data %>% select(Participant_id, tp_lab, Infection_Status))
    g <- ggplot(ord, aes(PC1, PC2, color = Infection_Status)) +
      geom_point(alpha = 0.85, size = 2) +
      labs(title = "PCoA (Jaccard) of VF presence/absence",
           subtitle = "Colored by Infection_Status (if available)")
    safe_ggsave(file.path(plots_dir, "pcoa_jaccard_status.png"), g, width = 7, height = 5, dpi = 300)
    # Save ordination coordinates for reproducibility
    readr::write_csv(ord, file.path(plots_dir, "pcoa_coordinates.csv"))
  } else message("↪  Too few samples for PCoA.")
} else {
  message("↪  vegan not installed; skipping PCoA.")
}

## ---------------- 8 · Gene–gene co-occurrence heatmap ------------------------
# choose variable genes (10–90% prev), top 60 by variance proxy p*(1-p)
prev_tbl <- vf_pa_all %>%
  select(-Participant_id, -tp_lab) %>%
  summarise(across(everything(), ~ mean(. > 0))) %>%
  pivot_longer(everything(), names_to = "GENE", values_to = "prev")

sel_genes <- prev_tbl %>%
  filter(prev >= 0.10, prev <= 0.90) %>%
  mutate(score = prev * (1 - prev)) %>%
  arrange(desc(score)) %>%
  slice_head(n = 60) %>% pull(GENE)

if (length(sel_genes) >= 4) {
  B <- vf_pa_all %>% select(all_of(sel_genes))
  # Prefer Jaccard via vegan if available; else fallback to binary distance
  if (has_pkg("vegan")) {
    use_pkg("vegan")
    D <- vegan::vegdist(t(B > 0), method = "jaccard", binary = TRUE)
    Sim <- as.matrix(1 - D)
  } else {
    D <- stats::dist(t(B > 0), method = "binary")
    Sim <- 1 - as.matrix(D)  # not true Jaccard, but acceptable fallback
  }
  hc  <- hclust(as.dist(1 - Sim), method = "average")
  ord <- hc$labels[hc$order]
  mat <- Sim[ord, ord]
  df  <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(df) <- c("Gene1","Gene2","Jaccard")
  
  g <- ggplot(df, aes(Gene1, Gene2, fill = Jaccard)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = "Gene–gene co-occurrence (Jaccard similarity)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6))
  safe_ggsave(file.path(plots_dir, "heatmap_gene_cooccurrence_top60.png"), g, width = 8, height = 8, dpi = 300)
} else {
  message("↪  Not enough variable genes for co-occurrence heatmap.")
}

## ---------------- 9 · Persistence vs transient summary -----------------------
persist_tbl <- vf_hits_all %>%
  distinct(Participant_id, GENE, tp_lab = tp_norm(Timepoint)$tp_lab) %>%
  group_by(Participant_id, GENE) %>%
  summarise(n_tp = n_distinct(tp_lab), .groups = "drop")

tp_counts <- vf_hits_all %>%
  distinct(Participant_id, tp_lab = tp_norm(Timepoint)$tp_lab) %>%
  count(Participant_id, name = "n_tp") %>% filter(n_tp >= 2)

stab <- persist_tbl %>%
  semi_join(tp_counts, by = "Participant_id") %>%
  group_by(Participant_id) %>%
  summarise(
    genes_total      = n(),
    genes_persistent = sum(n_tp == max(n_tp)),
    genes_transient  = sum(n_tp <  max(n_tp)),
    frac_persistent  = round(100 * genes_persistent / genes_total, 1),
    .groups = "drop"
  ) %>% arrange(desc(frac_persistent))

write_csv(stab, "results/persistence_per_participant.csv")

g <- ggplot(stab, aes(reorder(Participant_id, frac_persistent), frac_persistent)) +
  geom_col() + coord_flip() +
  labs(title = "Share of persistent genes per participant", x = "Participant", y = "% persistent genes")
safe_ggsave(file.path(plots_dir, "bar_persistence_share_by_participant.png"), g, width = 7, height = 8, dpi = 300)

## ---------------- 10 · Done --------------------------------------------------
message("✓  All plots complete – check results/plots/ and CSVs in results/.")
#!/usr/bin/env Rscript
# =============================================================
# 04_gene_breakdown.R  
# -------------------------------------------------------------
# Reads : results/vf_hits_all.rds   (long VFDB-hit table from 02)
# Writes :
#   results/annotated_gene_table.csv
#   results/per_sample_category_counts.csv
#   results/nitrate_presence_matrix.csv
#   results/plots/nitrate_upset.png
#
# Part 2  adds:
#   • status-stratified prevalence & Fisher tests for focus genes
#   • ST-adjusted enrichment via per-gene logistic regression (glm) + volcano
#
# Changes included:
#  • Imports broom (for tidying glm )
#  • Fix: build starter gene_map *after* loading vf_hits_all
#  • Prevent: the canonical tribble never overwrites an existing CSV map
#    (we define gene_map_default here, and only fall back to it if no CSV exists)
#  • Status-stratified prevalence & Fisher tests for focus genes
#  • ST-adjusted enrichment via per-gene logistic regression (glm)
#  • Volcano plot of ST-adjusted results
# =============================================================

## ---------- 1 · libraries, dirs, theme --------------------------------------
suppressPackageStartupMessages({
  library(dplyr);   library(tidyr);      library(readr)
  library(tibble);  library(ggplot2);    library(ComplexUpset)
  library(stringr); library(purrr);      library(forcats)
  library(broom)   # for tidying glm (used in Part 2)
})

dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)

# ggplot2 3.5+ friendly legend alignment
theme_set(theme_minimal(base_size = 11))
theme_update(legend.title = element_text(hjust = 0),
             legend.text  = element_text(hjust = 0))

## ---------- 2 · helpers ------------------------------------------------------
# Uniform timepoint labels (‘T0’, ‘T1’, …, ‘Uricult’, ‘Unscheduled’)
tp_norm <- function(x) {
  tp_chr     <- as.character(x)
  is_uricult <- str_detect(tp_chr, regex("uricult", ignore_case = TRUE))
  tp_num     <- suppressWarnings(as.integer(str_extract(tp_chr, "\\d+")))
  tp_lab     <- dplyr::case_when(
    is_uricult     ~ "Uricult",
    !is.na(tp_num) ~ paste0("T", tp_num),
    TRUE           ~ "Unscheduled"
  )
  tp_levels <- c(paste0("T", sort(unique(tp_num[!is.na(tp_num)]))),
                 "Uricult", "Unscheduled")
  tibble(tp_lab = factor(tp_lab, levels = tp_levels), tp_num = tp_num)
}

# Canonicalise gene strings (case/space/punct insensitive)
canon <- function(x) stringr::str_to_lower(gsub("[^a-z0-9]+", "", x))

# Regex alias rules (first match wins)
alias_rules <- tibble::tribble(
  ~pattern,                                   ~Gene,
  "^bla[-_ ]?ctx[-_ ]?m[-_ ]?15(\\.[0-9]+)?$",  "blaCTX-M-15",
  "^bla[-_ ]?ctx[-_ ]?m[-_ ]?[0-9]+[a-zA-Z]?$", "blaCTX-M",
  "^bla[-_ ]?oxa[-_ ]?48$",                     "blaOXA-48",
  "^bla[-_ ]?oxa[-_ ]?181$",                    "blaOXA-181",
  "^bla[-_ ]?oxa[-_ ]?232$",                    "blaOXA-232",
  "^bla[-_ ]?oxa[-_ ]?244$",                    "blaOXA-244",
  "^aac\\(6'\\)[-_ ]?ib[-_ ]?cr[0-9]*$",        "aac(6')-Ib-cr",
  "^aac\\(6'\\)[-_ ]?ib($|[-_ ].*)",            "aac(6')-Ib",
  "^aac\\(3\\).*$",                              "aac(3)",
  "^qnrA[0-9]*$",                               "qnrA",
  "^qnrB[0-9]*$",                               "qnrB",
  "^qnrS[0-9]*$",                               "qnrS",
  "^dfrA[-_ ]?1(\\b|\\D|$)",                    "dfrA1",
  "^dfrA[-_ ]?17(\\b|\\D|$)",                   "dfrA17",
  "^dfrA[0-9]+.*$",                              "dfrA",
  "^mcr[-_ ]?([1-9]|10)(\\.[0-9]+)?$",          "mcr-\\1",
  "^gyrA(\\b|[_-].*)$",                         "gyrA",
  "^parC(\\b|[_-].*)$",                         "parC",
  "^sfa[/_-]?focde$",                           "sfa/focDE",
  "^afa[/_-]?drabc$",                           "afa/draBC",
  "^kpsm[-_ ]?ii$",                              "kpsM II",
  "^kpsmt[-_ ]?ii$",                             "kpsMT II",
  "^repfia.*$",                                  "repFIA",
  "^repfii.*$",                                  "repFII",
  "^is1(\\b|[_-].*)$",                           "IS1",
  "^is3(\\b|[_-].*)$",                           "IS3",
  "^is91(\\b|[_-].*)$",                          "IS91",
  "^iscr1(\\b|[_-].*)$",                         "ISCR1",
  "^tn21(\\b|[_-].*)$",                          "Tn21"
)

apply_aliases <- function(df, gene_col = "Gene") {
  x <- as.character(df[[gene_col]])
  seen <- rep(FALSE, length(x))
  for (i in seq_len(nrow(alias_rules))) {
    pat  <- alias_rules$pattern[i]
    repl <- alias_rules$Gene[i]
    idx  <- !seen & grepl(pat, x, ignore.case = TRUE, perl = TRUE)
    if (any(idx, na.rm = TRUE)) {
      x[idx] <- gsub(pat, repl, x[idx], ignore.case = TRUE, perl = TRUE)
      seen[idx] <- TRUE
    }
  }
  df[[gene_col]] <- x
  df
}

## ---------- 3 · load VFDB hits (BEFORE gene_map) ------------------------------
vf_hits_all <- readRDS("results/vf_hits_all.rds")   # from script 02

# Ensure we have a 'Gene' column; add tp_lab only if missing
hits <- vf_hits_all
gene_col <- intersect(c("GENE","gene","Gene"), names(hits))[1]
if (is.na(gene_col)) stop("Could not find a gene column among: GENE/gene/Gene")
names(hits)[names(hits) == gene_col] <- "Gene"

if (!"tp_lab" %in% names(hits)) {
  hits <- dplyr::bind_cols(hits, tp_norm(hits$Timepoint))
}

# trim whitespace in Gene + apply regex aliases
hits$Gene <- trimws(hits$Gene)
hits <- apply_aliases(hits, "Gene")

# Convenience: add a canonicalised gene column
hits_nc <- hits %>% mutate(Gene_clean = canon(Gene))

## ---------- 4 · gene_map (CSV > starter > default) ---------------------------
gene_map_path <- "results/gene_map.csv"

# (A) Your canonical annotations — preserved as the *default* (won't overwrite CSV)
gene_map_default <- tibble::tribble(
  ~Gene,          ~Category,                ~Subcategory,
  # ---- Virulence (focus + core) ----
  "afa/draBC",    "Virulence factor",       "Adhesin (Dr-binding fimbriae)",
  "bmaE",         "Virulence factor",       "Adhesin (M hemagglutinin)",
  "cvaC",         "Virulence factor",       "Bacteriocin (Colicin V production)",
  "f17A",         "Virulence factor",       "Adhesin (F17 fimbriae)",
  "fimG",         "Virulence factor",       "Type 1 fimbriae (assembly)",
  "fimH",         "Virulence factor",       "Adhesin (Type 1 fimbriae)",
  "focG",         "Virulence factor",       "Adhesin (F1C fimbriae)",
  "hlyA",         "Virulence factor",       "Toxin (alpha-hemolysin)",
  "ibeA",         "Virulence factor",       "Invasin (brain endothelium)",
  "iha",          "Virulence factor",       "Adhesin (iron-regulated)",
  "iroN",         "Virulence factor",       "Siderophore receptor (iron acquisition)",
  "kpsM II",      "Virulence factor",       "Capsule synthesis (Group II)",
  "kpsMT II",     "Virulence factor",       "Capsule transport (Group II)",
  "malX",         "Virulence factor",       "PAI-associated PTS",
  "neuA",         "Virulence factor",       "Capsule (sialic acid metabolism)",
  "ompT",         "Virulence factor",       "Outer membrane protease",
  "papA",         "Virulence factor",       "Adhesin (P fimbriae)",
  "papC",         "Virulence factor",       "P fimbriae assembly (usher)",
  "papEF",        "Virulence factor",       "P fimbriae accessory",
  "pic",          "Virulence factor",       "Protease (colonization)",
  "sfa/focDE",    "Virulence factor",       "Adhesin (S/F1C fimbriae)",
  "sfaS",         "Virulence factor",       "Adhesin (S fimbriae)",
  "sat",          "Virulence factor",       "Toxin (autotransporter)",
  "tcpC",         "Virulence factor",       "Immune modulation (TLR antagonist)",
  "traT",         "Virulence factor",       "Serum resistance (surface exclusion)",
  "tsh",          "Virulence factor",       "Adhesin (temperature-sensitive hemagglutinin)",
  "usp",          "Virulence factor",       "Uropathogenic-specific protein",
  "vat",          "Virulence factor",       "Toxin (vacuolating autotransporter)",
  "yfcV",         "Virulence factor",       "Adhesin (fimbrial-like protein)",
  "chuA",         "Virulence factor",       "Heme uptake receptor",
  "fyuA",         "Virulence factor",       "Yersiniabactin receptor",
  "iucD",         "Virulence factor",       "Aerobactin synthesis",
  # ---- Antibiotic resistance ----
  "blaOXA-48",    "Antibiotic resistance",  "Beta-lactamase (OXA-48-like)",
  "blaOXA-181",   "Antibiotic resistance",  "Beta-lactamase (OXA-48 variant)",
  "blaOXA-232",   "Antibiotic resistance",  "Beta-lactamase (OXA-48 variant)",
  "blaOXA-244",   "Antibiotic resistance",  "Beta-lactamase (OXA-48 variant)",
  "blaCTX-M-15",  "Antibiotic resistance",  "ESBL (cephalosporins)",
  "blaCTX-M",     "Antibiotic resistance",  "ESBL (CTX-M family)",
  "blaTEM",       "Antibiotic resistance",  "Beta-lactamase",
  "blaSHV",       "Antibiotic resistance",  "Beta-lactamase",
  "aac(6')-Ib-cr","Antibiotic resistance",  "PMQR / aminoglycoside",
  "aac(6')-Ib",   "Antibiotic resistance",  "Aminoglycoside acetyltransferase",
  "aac(3)",       "Antibiotic resistance",  "Aminoglycoside acetyltransferase (AAC3 family)",
  "aadA",         "Antibiotic resistance",  "Aminoglycoside adenyltransferase",
  "aph(3')-Ia",   "Antibiotic resistance",  "Aminoglycoside phosphotransferase",
  "aph(3')",      "Antibiotic resistance",  "Aminoglycoside phosphotransferase (family)",
  "mph(A)",       "Antibiotic resistance",  "Macrolide phosphotransferase",
  "dfrA1",        "Antibiotic resistance",  "Trimethoprim (DHFR A1)",
  "dfrA17",       "Antibiotic resistance",  "Trimethoprim (DHFR A17)",
  "dfrA",         "Antibiotic resistance",  "Trimethoprim (DHFR family)",
  "sul1",         "Antibiotic resistance",  "Sulfonamide",
  "sul2",         "Antibiotic resistance",  "Sulfonamide",
  "sul3",         "Antibiotic resistance",  "Sulfonamide",
  "tetA",         "Antibiotic resistance",  "Tetracycline efflux",
  "tetB",         "Antibiotic resistance",  "Tetracycline efflux",
  "catA1",        "Antibiotic resistance",  "Chloramphenicol acetyltransferase",
  "floR",         "Antibiotic resistance",  "Phenicols",
  "qnrA",         "Antibiotic resistance",  "PMQR (qnrA)",
  "qnrB",         "Antibiotic resistance",  "PMQR (qnrB)",
  "qnrS",         "Antibiotic resistance",  "PMQR (qnrS)",
  "gyrA",         "Antibiotic resistance",  "Quinolone (QRDR site)",
  "parC",         "Antibiotic resistance",  "Quinolone (QRDR site)",
  "mdf(A)",       "Antibiotic resistance",  "Multidrug efflux (chromosomal)",
  "mcr-1",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-2",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-3",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-4",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-5",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-6",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-7",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-8",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-9",        "Antibiotic resistance",  "Colistin (mobile)",
  "mcr-10",       "Antibiotic resistance",  "Colistin (mobile)",
  # ---- Nitrate metabolism ----
  "fnr",          "Nitrate metabolism",     "Regulator",
  "narG",         "Nitrate metabolism",     "Nar",
  "narH",         "Nitrate metabolism",     "Nar",
  "narI",         "Nitrate metabolism",     "Nar",
  "narJ",         "Nitrate metabolism",     "Nar",
  "narK",         "Nitrate metabolism",     "Transport",
  "narL",         "Nitrate metabolism",     "Regulator",
  "narP",         "Nitrate metabolism",     "Regulator",
  "narQ",         "Nitrate metabolism",     "Regulator",
  "narX",         "Nitrate metabolism",     "Regulator",
  "napA",         "Nitrate metabolism",     "Nap",
  "napB",         "Nitrate metabolism",     "Nap",
  "napC",         "Nitrate metabolism",     "Nap",
  "napD",         "Nitrate metabolism",     "Nap",
  "napF",         "Nitrate metabolism",     "Nap",
  "napG",         "Nitrate metabolism",     "Nap",
  "napH",         "Nitrate metabolism",     "Nap",
  "nasA",         "Nitrate metabolism",     "Nas",
  "nasB",         "Nitrate metabolism",     "Nas",
  "nasC",         "Nitrate metabolism",     "Nas",
  "nasG",         "Nitrate metabolism",     "Nas",
  "nasH",         "Nitrate metabolism",     "Nas",
  "nirB",         "Nitrate metabolism",     "Nitrite reductase",
  "nirD",         "Nitrate metabolism",     "Nitrite reductase",
  "nrfA",         "Nitrate metabolism",     "Nitrite reductase",
  "nrfB",         "Nitrate metabolism",     "Nitrite reductase",
  "norR",         "Nitrate metabolism",     "NO detox regulator",
  "norV",         "Nitrate metabolism",     "NO detox",
  "norW",         "Nitrate metabolism",     "NO detox",
  # ---- Mobile / plasmid ----
  "IS1",          "Mobile element",         "Transposase",
  "IS3",          "Mobile element",         "Transposase",
  "IS91",         "Mobile element",         "Transposase",
  "ISCR1",        "Mobile element",         "ISCR element",
  "Tn21",         "Mobile element",         "Transposon",
  "ccdB",         "Mobile element",         "Addiction toxin",
  "hok",          "Mobile element",         "Addiction toxin",
  "repFIA",       "Mobile element",         "Plasmid replicon (IncFIA)",
  "repFII",       "Mobile element",         "Plasmid replicon (IncFII)"
)

# (B) CSV takes precedence; otherwise, build a starter map from this dataset
if (file.exists(gene_map_path)) {
  gene_map <- readr::read_csv(gene_map_path, show_col_types = FALSE) %>%
    mutate(
      Gene        = as.character(Gene),
      Category    = dplyr::coalesce(as.character(Category), "Unassigned"),
      Subcategory = dplyr::coalesce(as.character(Subcategory), Category)
    )
  # Optional (disabled): augment missing entries from defaults without overwriting
  # missing_defaults <- anti_join(gene_map_default, gene_map, by = "Gene")
  # if (nrow(missing_defaults)) {
  #   gene_map <- bind_rows(gene_map, missing_defaults)
  #   readr::write_csv(gene_map, gene_map_path)  # update with additions only
  # }
} else {
  message("↪  No results/gene_map.csv found; creating a starter map…")
  genes_seen <- hits %>%
    rename(Gene = dplyr::any_of("Gene")) %>%
    mutate(Gene = trimws(as.character(Gene))) %>%
    distinct(Gene) %>% filter(!is.na(Gene), Gene != "")
  
  # Heuristic for uncatalogued genes; prefer defaults where available
  starter <- genes_seen %>%
    left_join(gene_map_default, by = "Gene") %>%
    mutate(
      Category = dplyr::coalesce(
        Category,
        case_when(
          str_detect(Gene, regex("^(fim|fml|pil|foc|sfa|pap|afa|dra|cfa)", TRUE)) ~ "Adhesion/Fimbriae",
          str_detect(Gene, regex("^(kps|kfi|neu|ugd|rmpA|caps|wzx|wzy)", TRUE))    ~ "Capsule/Surface",
          str_detect(Gene, regex("^(iut|iuc|iro|irp|fyuA|chu|fep|ent|fec|ybt)", TRUE)) ~ "Iron acquisition/Siderophores",
          str_detect(Gene, regex("^(hly|cnf|sat|vat|cdt|astA|subAB|stx|lt|st)", TRUE)) ~ "Toxins",
          str_detect(Gene, regex("^(omp|iss|ibe|tra|usp|malX)", TRUE))            ~ "Invasion/Evasion",
          str_detect(Gene, regex("^(bla|qnr|aac|aph|aad|erm|cat|tet|sul|dfr|mcr|gyrA|parC)", TRUE)) ~ "AMR",
          TRUE ~ "Unassigned"
        )
      ),
      Subcategory = dplyr::coalesce(Subcategory, Category)
    )
  
  gene_map <- starter
  readr::write_csv(gene_map, gene_map_path)
  message("✓  Starter mapping written to ", gene_map_path,
          " — review/edit categories there any time and re-run.")
}

## ---------- 5 · annotate & per-sample category counts ------------------------
annotated <- hits %>%
  left_join(gene_map, by = "Gene") %>%
  mutate(
    Category    = dplyr::coalesce(Category,    "Unassigned"),
    Subcategory = dplyr::coalesce(Subcategory, "Unassigned")
  )

readr::write_csv(annotated, "results/annotated_gene_table.csv")

cat_summary <- annotated %>%
  group_by(Participant_id, Timepoint, Category) %>%
  summarise(n_genes = n_distinct(Gene), .groups = "drop") %>%
  pivot_wider(names_from = Category, values_from = n_genes, values_fill = 0)
readr::write_csv(cat_summary, "results/per_sample_category_counts.csv")

## ---------- 6 · nitrate-system snapshot --------------------------------------
nitrate_sys <- list(
  Nar = c("narG","narH","narJ","narI"),
  Nap = c("napF","napD","napA","napG","napH","napB","napC"),
  Nas = c("nasA","nasB","nasG","nasH","nasC")
)

nitrate_long <- hits_nc %>%
  filter(Gene_clean %in% canon(unlist(nitrate_sys))) %>%
  mutate(System = case_when(
    Gene_clean %in% canon(nitrate_sys$Nar) ~ "Nar",
    Gene_clean %in% canon(nitrate_sys$Nap) ~ "Nap",
    Gene_clean %in% canon(nitrate_sys$Nas) ~ "Nas"
  )) %>%
  distinct(Participant_id, tp_lab, System)

if (nrow(nitrate_long)) {
  nitrate_mat <- nitrate_long %>%
    mutate(present = 1L) %>%
    pivot_wider(names_from = System, values_from = present, values_fill = 0L)
  
  for (col in c("Nar","Nap","Nas"))
    if (!col %in% names(nitrate_mat)) nitrate_mat[[col]] <- 0L
  
  readr::write_csv(nitrate_mat, "results/nitrate_presence_matrix.csv")
  
  if (any(colSums(nitrate_mat[, c("Nar","Nap","Nas"), drop = FALSE]) > 0)) {
    ups_df <- nitrate_mat %>%
      unite(Sample, Participant_id, tp_lab, sep = "_") %>%
      mutate(across(c(Nar, Nap, Nas), as.logical))
    p <- ComplexUpset::upset(
      ups_df, intersect = c("Nar", "Nap", "Nas"),
      name = "N-system\npresence", min_size = 1
    )
    ggsave("results/plots/nitrate_upset.png", p, width = 6, height = 4, dpi = 300)
  }
} else {
  message("↪  No Nar/Nap/Nas genes found – skipping nitrate analysis")
}

## ---------- 7 · focus genes (exact + family regex) ---------------------------
# exact list (kept as in your script)
focus_exact <- c(
  # Virulence
  "afa/draBC","bmaE","cvaC","f17A","fimG","fimH","focG","hlyA","ibeA","iha",
  "iroN","kpsM II","kpsMT II","malX","neuA","ompT","papA","papC","papEF","pic",
  "sfa/focDE","sat","tcpC","traT","tsh","usp","vat","yfcV","chuA","fyuA","iucD",
  # AMR
  "blaOXA-48","blaOXA-181","blaOXA-244","blaOXA-232","blaCTX-M-15","blaTEM","blaSHV",
  "aadA","aac(3)","aph(3')","mph(A)","dfrA","sul1","sul2","sul3","tetA","tetB","catA1",
  "floR","qnrA","qnrB","qnrS","aac(6')-Ib-cr","gyrA","parC","mdf(A)",
  "mcr-1","mcr-2","mcr-3","mcr-4","mcr-5","mcr-6","mcr-7","mcr-8","mcr-9","mcr-10"
)

# family regex (operate on canonical strings; first match wins)
family_map <- list(
  "qnrA" = "^qnra", "qnrB" = "^qnrb", "qnrS" = "^qnrs",
  "aac(3)" = "^aac3", "aph(3')" = "^aph3",
  "dfrA" = "^dfra",
  "tetA" = "^teta", "tetB" = "^tetb",
  "mcr-1" = "^mcr1(\\d*)$", "mcr-2" = "^mcr2(\\d*)$", "mcr-3" = "^mcr3(\\d*)$",
  "mcr-4" = "^mcr4(\\d*)$", "mcr-5" = "^mcr5(\\d*)$", "mcr-6" = "^mcr6(\\d*)$",
  "mcr-7" = "^mcr7(\\d*)$", "mcr-8" = "^mcr8(\\d*)$", "mcr-9" = "^mcr9(\\d*)$",
  "mcr-10"= "^mcr10(\\d*)$"
  # keep exact-only for blaCTX-M-15, blaOXA-48/181/232/244, aac(6')-Ib-cr, gyrA, parC
)

focus_exact_canon <- stats::setNames(focus_exact, canon(focus_exact))

map_focus <- function(gene_clean) {
  res <- rep(NA_character_, length(gene_clean))
  # exact matches
  m <- match(gene_clean, names(focus_exact_canon))
  hit_exact <- !is.na(m)
  res[hit_exact] <- unname(focus_exact_canon[m[hit_exact]])
  # family fallbacks (fill only remaining NAs)
  remaining <- which(is.na(res))
  if (length(remaining)) {
    for (label in names(family_map)) {
      pat <- family_map[[label]]
      idx <- remaining[ grepl(pat, gene_clean[remaining], ignore.case = TRUE, perl = TRUE) ]
      if (length(idx)) {
        res[idx] <- label
        remaining <- setdiff(remaining, idx)
        if (!length(remaining)) break
      }
    }
  }
  res
}

# apply mapping → per-sample presence (one row per sample-gene)
hits_focus <- hits_nc %>%
  mutate(FocusKey = map_focus(Gene_clean)) %>%
  filter(!is.na(FocusKey))

focus_pa <- hits_focus %>%
  distinct(Participant_id, tp_lab, FocusKey) %>%
  mutate(present = 1L)


## ---------- 8 · status map join & prevalence / Fisher -------------------------
if (file.exists("results/status_map.csv")) {
  status_map <- readr::read_csv("results/status_map.csv", show_col_types = FALSE)
  
  if (!"tp_lab" %in% names(status_map)) {
    if (!"Timepoint" %in% names(status_map)) {
      stop("status_map.csv must contain 'tp_lab' or 'Timepoint'.")
    }
    status_map <- dplyr::bind_cols(status_map, tp_norm(status_map$Timepoint))
  }
  
  status_map <- status_map %>%
    dplyr::select(Participant_id, tp_lab, Infection_Status) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      Infection_Status = forcats::fct_relevel(as.factor(Infection_Status), "ASB", "UTI")
    )
  
  # sample × gene matrix (0/1)
  focus_wide <- focus_pa %>%
    tidyr::pivot_wider(names_from = FocusKey, values_from = present, values_fill = 0L)
  
  samples_focus <- status_map %>%
    dplyr::left_join(focus_wide, by = c("Participant_id","tp_lab"))
  
  focus_cols <- setdiff(names(samples_focus), c("Participant_id","tp_lab","Infection_Status"))
  if (!length(focus_cols)) {
    message("↪  No focus-gene columns found; skipping status-stratified analysis.")
  } else {
    samples_focus <- samples_focus %>%
      dplyr::mutate(across(all_of(focus_cols), ~ tidyr::replace_na(., 0L)))
    
    vt <- samples_focus %>%
      tidyr::pivot_longer(all_of(focus_cols), names_to = "FocusKey", values_to = "present") %>%
      dplyr::mutate(present = present > 0) %>%
      dplyr::filter(Infection_Status %in% c("UTI","ASB"))
    
    # Prevalence by status (+ CSV + plot)
    prev_tbl <- vt %>%
      dplyr::group_by(FocusKey, Infection_Status) %>%
      dplyr::summarise(n_pos = sum(present), n_total = dplyr::n(),
                       prevalence = n_pos / n_total, .groups = "drop")
    
    readr::write_csv(prev_tbl, "results/vf_focus_prevalence_by_status.csv")
    
    p_prev <- ggplot2::ggplot(
      prev_tbl,
      ggplot2::aes(forcats::fct_reorder(FocusKey, prevalence), prevalence, fill = Infection_Status)
    ) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      ggplot2::labs(title = "Focus-gene prevalence by status", x = NULL, y = "Prevalence") +
      ggplot2::theme_minimal(base_size = 11)
    
    ggplot2::ggsave("results/plots/vf_focus_prevalence_by_status.png",
                    p_prev, width = 8, height = 6, dpi = 300)
    
    # Fisher's exact test per gene
    gene_enrichment_focus <- vt %>%
      dplyr::count(FocusKey, Infection_Status, present, name = "n") %>%
      tidyr::complete(FocusKey,
                      Infection_Status = c("UTI","ASB"),
                      present = c(FALSE, TRUE),
                      fill = list(n = 0L)) %>%
      tidyr::pivot_wider(names_from = c(Infection_Status, present), values_from = n) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        p  = fisher.test(matrix(c(`UTI_TRUE`, `UTI_FALSE`, `ASB_TRUE`, `ASB_FALSE`), nrow = 2))$p.value,
        OR = ((`UTI_TRUE` + 0.5)/(`UTI_FALSE` + 0.5)) /
          ((`ASB_TRUE` + 0.5)/(`ASB_FALSE` + 0.5))
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(p_adj = p.adjust(p, method = "BH")) %>%
      dplyr::arrange(p_adj)
    
    readr::write_csv(gene_enrichment_focus, "results/diff_focus_genes_UTI_vs_ASB_fisher.csv")
    ## ---------- 9 · ST-adjusted enrichment (per-gene GLM) --------------------
    # Try to add ST from MLST table; analysis still proceeds without it.
    st_attached <- FALSE
    if (file.exists("results/mlst/mlst_matrix.csv")) {
      mlst_tbl <- suppressMessages(readr::read_csv("results/mlst/mlst_matrix.csv", show_col_types = FALSE))
      if (!"tp_lab" %in% names(mlst_tbl) && "Timepoint" %in% names(mlst_tbl)) {
        mlst_tbl <- dplyr::bind_cols(mlst_tbl, tp_norm(mlst_tbl$Timepoint))
      }
      if (all(c("Participant_id","tp_lab","ST") %in% names(mlst_tbl))) {
        mlst_st <- mlst_tbl %>%
          dplyr::select(Participant_id, tp_lab, ST) %>%
          dplyr::mutate(ST = as.character(ST)) %>%
          dplyr::distinct()
        samples_focus <- samples_focus %>%
          dplyr::left_join(mlst_st, by = c("Participant_id","tp_lab"))
        st_attached <- TRUE
      } else {
        message("↪  MLST table lacks Participant_id/tp_lab/ST — skipping ST covariate.")
      }
    } else {
      message("↪  results/mlst/mlst_matrix.csv not found — skipping ST covariate.")
    }
    
    # Long table with ST included (if available)
    vt_st <- samples_focus %>%
      tidyr::pivot_longer(all_of(focus_cols), names_to = "FocusKey", values_to = "present") %>%
      dplyr::mutate(
        present = as.integer(present > 0),
        Infection_Status = forcats::fct_relevel(as.factor(Infection_Status), "ASB","UTI"),
        ST = if (st_attached) forcats::fct_na_value_to_level(as.factor(ST), level = "ST_NA") else factor("ST_missing")
      ) %>%
      dplyr::filter(Infection_Status %in% c("ASB","UTI"))
    
    # helper: robust per-gene logistic fit (handles separation)
    safe_logistic_fit <- function(dat, st_attached) {
      fml <- if (st_attached) present ~ Infection_Status + ST else present ~ Infection_Status
      
      fit <- try(suppressWarnings(stats::glm(fml, data = dat, family = stats::binomial(),
                                             control = list(maxit = 50))), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        td <- broom::tidy(fit)
        ok <- any(td$term == "Infection_StatusUTI") &&
          is.finite(td$estimate[td$term == "Infection_StatusUTI"]) &&
          is.finite(td$std.error[td$term == "Infection_StatusUTI"])
        if (ok) return(list(method = "glm", tidy = td))
      }
      
      if (requireNamespace("logistf", quietly = TRUE)) {
        fit2 <- try(logistf::logistf(fml, data = dat), silent = TRUE)
        if (!inherits(fit2, "try-error")) {
          if ("Infection_StatusUTI" %in% names(fit2$coef)) {
            ti <- "Infection_StatusUTI"
            td2 <- tibble::tibble(
              term      = ti,
              estimate  = unname(fit2$coef[ti]),
              std.error = unname(fit2$se[ti]),
              statistic = NA_real_,
              p.value   = unname(fit2$prob[ti])
            )
            return(list(method = "logistf", tidy = td2))
          }
        }
      }
      
      if (requireNamespace("brglm2", quietly = TRUE)) {
        fit3 <- try(suppressWarnings(stats::glm(fml, data = dat, family = stats::binomial(),
                                                method = brglm2::brglmFit)), silent = TRUE)
        if (!inherits(fit3, "try-error")) {
          td3 <- broom::tidy(fit3)
          if (any(td3$term == "Infection_StatusUTI")) {
            return(list(method = "brglm2", tidy = td3))
          }
        }
      }
      
      message("⚠️  GLM failed (and no robust fallback available) for a gene; returning NA.")
      list(method = "failed",
           tidy = tibble::tibble(term = "Infection_StatusUTI",
                                 estimate = NA_real_, std.error = NA_real_,
                                 statistic = NA_real_, p.value = NA_real_))
    }
    
    # per-gene model: present ~ Infection_Status (+ ST if present)
    glm_by_ST <- vt_st %>%
      dplyr::group_by(FocusKey) %>%
      dplyr::group_modify(function(dat, key) {
        if (dplyr::n_distinct(dat$Infection_Status) < 2 || dplyr::n_distinct(dat$present) < 2) {
          return(tibble::tibble(term = "Infection_StatusUTI",
                                estimate = NA_real_, std.error = NA_real_,
                                statistic = NA_real_, p.value = NA_real_, n = nrow(dat),
                                method = "insufficient-variation"))
        }
        res <- safe_logistic_fit(dat, st_attached)
        res$tidy %>%
          dplyr::filter(term == "Infection_StatusUTI") %>%
          dplyr::mutate(n = nrow(dat), method = res$method)
      }) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        OR        = exp(estimate),
        conf.low  = exp(estimate - 1.96*std.error),
        conf.high = exp(estimate + 1.96*std.error),
        log2OR    = log2(OR),
        p_adj     = p.adjust(p.value, method = "BH")
      ) %>%
      dplyr::arrange(p_adj)
    
    readr::write_csv(glm_by_ST, "results/diff_focus_genes_UTI_vs_ASB_glm_by_ST.csv")
    
    # Volcano (build first, then save; don't add ggsave with '+')
    volcano_df <- glm_by_ST %>% dplyr::filter(is.finite(log2OR), is.finite(p_adj))
    
    p_volc <- ggplot2::ggplot(volcano_df, ggplot2::aes(x = log2OR, y = -log10(p_adj))) +
      ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
      ggplot2::geom_point(ggplot2::aes(color = p_adj < 0.05), alpha = 0.8, size = 2) +
      ggplot2::scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "firebrick")) +
      ggplot2::labs(
        title = "ST-adjusted enrichment (UTI vs ASB)",
        x = "log2(OR) for UTI (vs ASB)",
        y = expression(-log[10](BH~adjusted~p)),
        color = "BH < 0.05"
      ) +
      ggplot2::theme_minimal(base_size = 11)
    
    ggplot2::ggsave("results/plots/vf_focus_volcano_glm_by_ST.png",
                    p_volc, width = 7, height = 5, dpi = 300)
  } # end else (focus_cols found)
} else {
  message("↪  results/status_map.csv not found — skipping status-stratified and ST-adjusted analyses.")
}
readr::write_csv(annotated, "results/annotated_gene_table.csv")
message("✓ wrote results/annotated_gene_table.csv (", nrow(annotated), " rows)")

#!/usr/bin/env Rscript
# =============================================================
#  05_gene_overview_plots.R
#  • bar-plot of gene prevalence
#  • heat-map of variable genes (all labels saved to PNG & PDF)
# =============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

# ---------------------------------------------------------------- paths -------
csv_in   <- "results/vf_pa_all.csv"
plot_dir <- "results/plots"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------- 1 · load -------
vf <- read_csv(csv_in, show_col_types = FALSE)

meta_cols  <- c("Participant_id", "Timepoint")
gene_cols  <- setdiff(names(vf), meta_cols)

# ----------------------------------------------------- 2 · tidy matrix -------
vf[gene_cols] <- vf[gene_cols] |>
  mutate(across(everything(), ~ replace_na(as.numeric(.x), 0)))

vf <- vf |>
  group_by(across(all_of(meta_cols))) |>
  summarise(across(all_of(gene_cols), max), .groups = "drop") |>
  arrange(across(all_of(meta_cols)))

# strictly-numeric presence/absence matrix
mat <- vf |>
  unite(row_id, Participant_id, Timepoint, sep = "_", remove = FALSE) |>
  column_to_rownames("row_id") |>
  select(all_of(gene_cols)) |>
  as.matrix()

# ------------------------------------------------ 3 · core vs variable --------
core_genes     <- names(which(colSums(mat) == nrow(mat)))
write_lines(core_genes, "results/core_gene_list.txt")

variable_genes <- setdiff(colnames(mat), core_genes)
write_lines(variable_genes, "results/variable_gene_list.txt")


cat("\n===== Core genes (present in *every* isolate) =====\n",
    paste(core_genes, collapse = ", "), "\n\n")
cat("===== Variable genes (present in < 100 % isolates) =====\n",
    paste(variable_genes, collapse = ", "), "\n\n")

# -------------------------------------- 4 · prevalence bar-plot --------------
gene_prev <- colSums(mat) |>
  enframe(name = "gene", value = "n_iso") |>
  mutate(prevalence = n_iso / nrow(mat)) |>
  arrange(desc(prevalence))

topN <- 40
prev_plot <- ggplot(slice_max(gene_prev, prevalence, n = topN),
                    aes(reorder(gene, prevalence), prevalence)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Isolates",
       title = paste("Top", topN, "genes by prevalence")) +
  theme_minimal(base_size = 10)

bar_file <- file.path(plot_dir, "gene_prevalence_bar.png")
ggsave(bar_file, prev_plot, width = 6, height = 8, dpi = 300)

# ------------------------------ 5 · variable-gene heat-map --------------------
var_mat <- mat[, variable_genes, drop = FALSE]
var_mat <- var_mat[, colSums(var_mat) > 0, drop = FALSE]   # safety

ann_row <- vf |>
  unite(row_id, Participant_id, Timepoint, sep = "_", remove = FALSE) |>
  select(row_id, Participant_id, Timepoint) |>
  column_to_rownames("row_id")

heat_file_png <- file.path(plot_dir, "variable_gene_heatmap.png")
heat_file_pdf <- file.path(plot_dir, "variable_gene_heatmap.pdf")  # wide PDF

# ---- PNG (bitmap) -----------------------------------------------------------
pheatmap(var_mat,
         cluster_rows   = FALSE,
         cluster_cols   = TRUE,
         show_rownames  = FALSE,
         fontsize_col   = 5,
         annotation_row = ann_row,
         filename       = heat_file_png)

# ---- PDF (vector: every label stays crisp) ----------------------------------
pdf(heat_file_pdf,
    width  = 0.18 * ncol(var_mat) + 2,   # ≈5 mm / gene + margins
    height = 8)
pheatmap(var_mat,
         cluster_rows   = FALSE,
         cluster_cols   = TRUE,
         show_rownames  = FALSE,
         fontsize_col   = 6,
         annotation_row = ann_row)
dev.off()

# --------------------------- 6 · optional on-screen preview -------------------
if (interactive()) {
  message("Showing heat-map on screen …")
  pheatmap(var_mat,
           cluster_rows   = FALSE,
           cluster_cols   = TRUE,
           show_rownames  = FALSE,
           fontsize_col   = 5,
           annotation_row = ann_row,
           main           = "Variable genes across participants & time-points")
}

# ------------------------------------------------------------------ done ------
message("\n✓  Outputs written:\n  • ", bar_file,
        "\n  • ", heat_file_png,
        "\n  • ", heat_file_pdf, "\n")
#!/usr/bin/env Rscript
# =============================================================
# 06_MLST.R  –  Multi-locus sequence typing for “Yellow RoUTIne”
# Works with mlst ≥ 2.19   ( --quiet --csv --legacy )
# Last updated: 8 Sep 2025  (QC + provenance + persistence)
# =============================================================

# ---- 0 · configuration ------------------------------------------------------
THREADS   <- pmin(4L, max(1L, parallel::detectCores() - 1L))   # PubMLST ≤ 4
FASTA_DIR <- "ont-yellow-routine-fastas"  # use repo-relative path
SCHEME    <- "ecoli"                      # PubMLST scheme, or "auto"
DEBUG     <- FALSE                        # TRUE → echo every mlst command

if (DEBUG) { options(warn = 1, error = recover); msg <- message } else {
  msg <- function(...) message(format(Sys.time(), "%H:%M:%S "), ...)
}

suppressPackageStartupMessages({
  library(dplyr);  library(readr);  library(tidyr)
  library(purrr);  library(furrr);  library(processx)
  library(stringr);library(fs);     library(scales); library(rlang)
  library(ggplot2)
})

msg("06_MLST.R starting")

plots_dir <- "results/plots"
dir_create(plots_dir, recurse = TRUE)

# ---------- 1 · read assembly metadata ---------------------------------------
assembly_df <- read_csv("assembly_metadata.csv", show_col_types = FALSE)
assembly_cols <- names(assembly_df)

if (!"full_path" %in% names(assembly_df)) {
  assembly_df <- assembly_df %>%
    mutate(full_path = file.path(FASTA_DIR, file_name)) %>%
    mutate(found     = file.exists(full_path))
  if (any(!assembly_df$found))
    stop("Missing FASTA files:\n",
         paste(assembly_df$full_path[!assembly_df$found], collapse = "\n"))
  assembly_df <- select(assembly_df, -found)
  msg("✓ full_path rebuilt (", nrow(assembly_df), " rows)")
}

# ---------- 2 · locate mlst & sanity-check -----------------------------------
MLST_BIN <- Sys.which("mlst")
if (MLST_BIN == "") stop("`mlst` not in $PATH – install or relink Homebrew")
if (system2(MLST_BIN, "--check") != 0L)
  stop("mlst --check failed – see stderr above")
msg("✓ mlst binary = ", MLST_BIN)

# ---------- 3 · output folders ------------------------------------------------
dir_create("results/mlst/raw",  recurse = TRUE)
dir_create("results/mlst/logs", recurse = TRUE)

# ---------- 4 · helpers -------------------------------------------------------
inspect_fail_logs <- function(n = 5) {
  fails <- fs::dir_ls("results/mlst/logs", glob = "*.log.txt")
  if (!length(fails)) { message("No failing logs found."); return(invisible()) }
  cat("\nShowing first", n, "failing logs:\n")
  for (lf in head(fails, n)) {
    cat("\n======", fs::path_file(lf), "======\n")
    cat(readLines(lf, n = 40), sep = "\n")
  }
  invisible()
}

tp_norm <- function(x) {
  x <- as.character(x)
  # fixed bug: ignore_case (not ignore.case)
  is_uricult <- stringr::str_detect(x, stringr::regex("uricult", ignore_case = TRUE))
  tp_num <- suppressWarnings(as.integer(stringr::str_extract(x, "\\d+")))
  tp_lab <- dplyr::case_when(
    is_uricult ~ "Uricult",
    !is.na(tp_num) ~ paste0("T", tp_num),
    TRUE ~ "Unscheduled"
  )
  tibble::tibble(
    tp_lab = factor(
      tp_lab,
      levels = c(paste0("T", sort(unique(tp_num[!is.na(tp_num)]))), "Uricult", "Unscheduled")
    )
  )
}


run_mlst <- function(fasta, scheme = SCHEME) {
  cache   <- fs::path("results/mlst/raw",  paste0(fs::path_file(fasta), ".mlst.csv"))
  log_out <- fs::path("results/mlst/logs", paste0(fs::path_file(fasta), ".log.txt"))
  
  if (!fs::file_exists(cache)) {
    cmd <- c("--quiet", "--threads", "1",
             "--scheme", scheme, "--csv", "--legacy", fasta)
    if (DEBUG) msg("cmd: ", MLST_BIN, " ", paste(cmd, collapse = " "))
    
    px <- processx::run(MLST_BIN, cmd,
                        echo             = DEBUG,
                        stderr_to_stdout = TRUE,
                        error_on_status  = FALSE)
    
    write_lines(px$stdout, log_out)
    if (px$status != 0) return(tibble(.status = "failed"))
    write_lines(px$stdout, cache)
  }
  
  dat <- read_csv(cache,
                  col_types = cols(.default = col_character()),
                  na        = c("", "?"),
                  progress  = FALSE) %>%
    rename_with(tolower) %>%
    mutate(scheme = scheme, .before = 1)
  
  # Log any parsing problems for this FASTA
  prob <- readr::problems(dat)
  if (nrow(prob)) {
    readr::write_csv(prob,
                     fs::path("results/mlst/logs", paste0(fs::path_file(fasta), ".problems.csv"))
    )
  }
  
  if (nrow(dat) == 0)
    return(tibble(.status = "no-data"))
  
  st_cols <- names(dat)[tolower(names(dat)) == "st"]
  if (!length(st_cols))
    return(tibble(.status = "unexpected-format"))
  
  dplyr::mutate(dat, .status = "ok")
}

# ---------- 5 · run in parallel ----------------------------------------------
plan(multisession, workers = THREADS)

mlst_raw <- assembly_df %>%
  mutate(mlst = future_map(full_path, run_mlst, .progress = TRUE)) %>%
  unnest(mlst)

# ---------- 6 · concise overview + provenance ---------------------------------
summary_tbl <- mlst_raw %>%
  count(.status, name = "n") %>%
  mutate(pct = percent(n / sum(n)))

msg("------ run summary ------")
print(summary_tbl)
if (!any(mlst_raw$.status == "ok"))
  stop("No successful typings – call inspect_fail_logs()")

# Save provenance (session + tool versions)
readr::write_lines(capture.output(sessionInfo()), "results/mlst/sessionInfo.txt")
suppressWarnings(
  readr::write_lines(system2(MLST_BIN, "--version", stdout = TRUE),
                     "results/mlst/mlst_version.txt")
)
suppressWarnings(
  readr::write_lines(system2("blastn", "-version", stdout = TRUE),
                     "results/mlst/blastn_version.txt")
)

# ---------- 7 · final tidy-up -------------------------------------------------
st_col <- names(mlst_raw)[tolower(names(mlst_raw)) == "st"]
if (!length(st_col))
  stop("Could not find an ST column – names are:\n",
       paste(names(mlst_raw), collapse = ", "))

mlst_tbl <- mlst_raw %>%
  filter(.status == "ok") %>%
  select(-.status) %>%
  mutate(ST = .data[[st_col[1]]]) %>%
  select(-all_of(st_col))

# Ensure we have a stable sample identifier
if (!"Isolate_ID" %in% names(mlst_tbl)) {
  if ("file_name" %in% names(mlst_tbl)) {
    mlst_tbl <- mlst_tbl %>%
      mutate(Isolate_ID = tools::file_path_sans_ext(basename(.data[["file_name"]])))
  } else if ("full_path" %in% names(mlst_tbl)) {
    mlst_tbl <- mlst_tbl %>%
      mutate(Isolate_ID = tools::file_path_sans_ext(basename(.data[["full_path"]])))
  } else {
    mlst_tbl <- mlst_tbl %>% mutate(Isolate_ID = row_number())
  }
}
mlst_tbl <- mlst_tbl %>% relocate(Isolate_ID, ST, .before = 1)

# ---------- 7b · MLST locus set + QC flags -----------------------------------
# Locus columns = everything not in assembly metadata or obvious meta fields
meta_cols <- unique(c("scheme","ST","Isolate_ID","file_name","full_path", assembly_cols))
locus_cols <- setdiff(names(mlst_tbl), meta_cols)

if (!length(locus_cols)) {
  warning("No MLST locus columns were detected; check input formatting.")
} else {
  write_lines(locus_cols, "results/mlst/mlst_locus_list.txt")
}

mlst_tbl <- mlst_tbl %>%
  mutate(across(all_of(locus_cols), as.character)) %>%
  rowwise() %>%
  mutate(
    n_loci_typed    = if (length(locus_cols))
      sum(!is.na(c_across(all_of(locus_cols))) &
            c_across(all_of(locus_cols)) != "" &
            !grepl("^(\\?|0)$", c_across(all_of(locus_cols))))
    else 0L,
    prop_loci_typed = if (length(locus_cols)) n_loci_typed / length(locus_cols) else NA_real_,
    mlst_complete   = if (length(locus_cols)) n_loci_typed == length(locus_cols) else NA,
    has_new_allele  = if (length(locus_cols)) any(grepl("NEW|\\*$", c_across(all_of(locus_cols)), ignore.case = TRUE)) else NA,
    ambiguous_call  = if (length(locus_cols)) any(grepl("[,;/]", c_across(all_of(locus_cols)))) else NA
  ) %>%
  ungroup()

mlst_qc <- mlst_tbl %>%
  count(mlst_complete, has_new_allele, ambiguous_call, name = "n")
write_csv(mlst_qc, "results/mlst/mlst_qc_summary.csv")

# ---------- 7c · Top STs (useful for downstream modeling) ---------------------
top_STs <- mlst_tbl %>%
  mutate(ST = as.character(ST)) %>%
  count(ST, sort = TRUE, name = "n_isolates") %>%
  mutate(pct = n_isolates / sum(n_isolates))
write_csv(top_STs, "results/mlst/top_STs.csv")

# ---------- 7d · Within-host ST persistence (if Participant_id/Timepoint) -----
has_pid <- "Participant_id" %in% names(mlst_tbl)
has_tp  <- "Timepoint"      %in% names(mlst_tbl)

if (has_pid && has_tp) {
  st_long <- mlst_tbl %>%
    select(Participant_id, Timepoint, ST, mlst_complete) %>%
    bind_cols(tp_norm(.$Timepoint)) %>%
    distinct(Participant_id, tp_lab, ST, mlst_complete)
  
  st_persist <- st_long %>%
    group_by(Participant_id) %>%
    summarise(
      n_tp        = n(),
      n_ST        = n_distinct(ST),
      dominant_ST = names(sort(table(ST), decreasing = TRUE))[1],
      frac_domin  = as.numeric(max(table(ST)))/n_tp,
      any_incomplete = any(!mlst_complete),
      .groups = "drop"
    ) %>%
    arrange(desc(frac_domin))
  write_csv(st_persist, "results/mlst/ST_persistence_by_participant.csv")
  
  # quick visual
  g <- ggplot(st_persist, aes(reorder(Participant_id, frac_domin), frac_domin)) +
    geom_col() + coord_flip() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Participant", y = "Fraction with dominant ST",
         title = "Within-host ST persistence")
  ggsave(file.path(plots_dir, "st_persistence_by_participant.png"),
         g, width = 7, height = 8, dpi = 300)
} else {
  msg("↪  Skipping ST persistence: Participant_id/Timepoint not found in MLST table.")
}

# ---------- 8 · write primary outputs -----------------------------------------
write_tsv(mlst_tbl,    "results/mlst/mlst_all.tsv")
write_csv(mlst_tbl,    "results/mlst/mlst_matrix.csv")
write_csv(summary_tbl, "results/mlst/log_summary.csv")

msg("✓ wrote ", nrow(mlst_tbl), " assemblies → results/mlst/")
msg("06_MLST.R finished ✅")
#!/usr/bin/env Rscript
# =============================================================
# 07_explore_mlst.R – quick descriptive stats & QC
# Requires the output of 06_MLST.R in results/mlst/
# =============================================================

## ---- 0 · libraries ---------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(fs)
  library(stringr)
  library(tidyr)
  library(scales)
})

## ---- 1 · file paths --------------------------------------------------------
mlst_file    <- "results/mlst/mlst_all.tsv"
meta_file    <- "assembly_metadata.csv"
plot_file    <- "results/mlst/top20_STs.pdf"
out_freq_csv <- "results/mlst/ST_frequencies.csv"

dir.create("results/debug", recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(plot_file), recursive = TRUE, showWarnings = FALSE)

## ---- 2 · load data ---------------------------------------------------------
if (!file.exists(mlst_file)) stop("Missing required file: ", mlst_file)
mlst <- read_tsv(mlst_file, show_col_types = FALSE)

# Ensure an ST column exists (tolerate different namings)
st_col <- names(mlst)[tolower(names(mlst)) == "st"]
if (!length(st_col)) {
  alt <- names(mlst)[grepl("(?i)^sequence[_ ]?type$", names(mlst))]
  if (length(alt)) st_col <- alt[1]
}
if (!length(st_col)) stop("Could not find an ST column in ", mlst_file,
                          ". Columns: ", paste(names(mlst), collapse = ", "))
if (!"ST" %in% names(mlst)) mlst$ST <- mlst[[st_col[1]]]
mlst <- mlst %>% mutate(ST = as.character(ST))

# Optional clinical metadata
if (file.exists(meta_file)) {
  meta <- read_csv(meta_file, show_col_types = FALSE)
} else {
  message("↪  ", meta_file, " not found; skipping meta-joins.")
}

## ---- 3 · collapse to one row per Isolate_ID (for counts & clarity) ----------
# Always keep a single representative row per Isolate_ID:
# preference: complete MLST -> more typed loci -> first deterministically.
mlst_in <- mlst %>%
  mutate(
    n_loci_typed  = dplyr::coalesce(.data[["n_loci_typed"]], -Inf),
    mlst_complete = dplyr::coalesce(.data[["mlst_complete"]], FALSE)
  ) %>%
  group_by(Isolate_ID) %>%
  arrange(desc(mlst_complete), desc(n_loci_typed), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  # drop helper cols if they weren't present in the raw table
  select(-any_of(c("n_loci_typed","mlst_complete")))

## ---- 3b · isolate-level ST frequency table ----------------------------------
mlst_for_counts <- mlst_in %>%
  mutate(ST = as.character(ST)) %>%
  filter(!is.na(ST), ST != "")

st_freq <- mlst_for_counts %>%
  count(ST, sort = TRUE) %>%
  mutate(pct = scales::percent(n / sum(n)))

write_csv(st_freq, out_freq_csv)
message("✓ wrote isolate-level ST frequency table → ", out_freq_csv)

## ---- 4 · bar-plot of top 20 STs (unique isolates) ---------------------------
top20 <- st_freq %>% slice_max(n, n = 20, with_ties = FALSE)

p <- ggplot(top20, aes(x = reorder(factor(ST), n), y = n)) +
  geom_col() +
  coord_flip() +
  labs(x = "Sequence type (ST)", y = "Isolates",
       title = "Top 20 E. coli STs - Yellow RoUTIne (unique isolates)") + # ASCII hyphen
  theme_classic(base_size = 12)

if (capabilities("cairo")) {
  ggsave(plot_file, p, width = 7, height = 5, device = cairo_pdf)
} else {
  ggsave(plot_file, p, width = 7, height = 5)
}
message("✓ saved plot → ", plot_file)

## ---- 5 · join ST back onto clinical meta (robust keys) ----------------------
if (exists("meta")) {
  iso_col <- names(meta)[grepl("(?i)^isolate[_ ]?id$", names(meta))][1]
  if (!is.na(iso_col) && "Isolate_ID" %in% names(mlst_in)) {
    
    # 5a) Diagnose duplicates in raw inputs (pre-collapse)
    mlst_key_diag <- mlst %>%
      group_by(Isolate_ID) %>%
      summarise(
        n_rows = n(),
        n_ST   = n_distinct(ST, na.rm = TRUE),
        ST_set = paste(sort(unique(na.omit(ST))), collapse = ","),
        .groups = "drop"
      )
    write_csv(mlst_key_diag, "results/debug/mlst_isolate_summary.csv")
    
    meta_key_diag <- meta %>%
      mutate(Isolate_ID = .data[[iso_col]]) %>%
      group_by(Isolate_ID) %>%
      summarise(
        n_rows = n(),
        n_pid  = dplyr::n_distinct(dplyr::coalesce(as.character(Participant_id), "NA")),
        n_tp   = dplyr::n_distinct(dplyr::coalesce(as.character(Timepoint),      "NA")),
        .groups = "drop"
      )
    write_csv(meta_key_diag, "results/debug/meta_isolate_summary.csv")
    
    # 5b) Decide join key
    can_use_isolate_only <- all(meta_key_diag$n_pid <= 1 & meta_key_diag$n_tp <= 1)
    
    # 5c) Build metadata with one row per Isolate_ID if using Isolate_ID-only
    if (can_use_isolate_only) {
      meta_1row <- meta %>%
        arrange(.data[[iso_col]], desc(!is.na(.data[["Timepoint"]]))) %>%
        distinct(.data[[iso_col]], .keep_all = TRUE)
      
      mlst_with_meta <- mlst_in %>%
        left_join(
          meta_1row,
          by = setNames("Isolate_ID", iso_col),
          relationship = "many-to-one"
        )
      
      branch <- "Isolate_ID only"
    } else {
      # Composite key path (requires both sides to have pid & timepoint)
      has_pid_tp_meta <- all(c("Participant_id","Timepoint") %in% names(meta))
      has_pid_tp_mlst <- all(c("Participant_id","Timepoint") %in% names(mlst_in))
      if (has_pid_tp_meta && has_pid_tp_mlst) {
        mlst_with_meta <- mlst_in %>%
          left_join(
            meta,
            by = c("Isolate_ID" = iso_col,
                   "Participant_id" = "Participant_id",
                   "Timepoint"      = "Timepoint"),
            relationship = "many-to-one"
          )
        branch <- "Composite (Isolate_ID + Participant_id + Timepoint)"
      } else {
        stop(
          "Isolate_ID maps to multiple participants/timepoints in metadata, ",
          "but MLST table lacks Participant_id/Timepoint. ",
          "Either (a) add those columns upstream, or (b) fix duplicated Isolate_IDs in meta."
        )
      }
    }
    
    write_csv(mlst_with_meta, "results/mlst/mlst_with_meta.csv")
    message("✓ wrote MLST + meta join → results/mlst/mlst_with_meta.csv")
    
    ## ---- Post-join audits -----------------------------------------------------
    writeLines(paste0("[07_explore_mlst] join branch: ", branch),
               "results/debug/07_join_branch.txt")
    
    # Isolates with >1 ST upstream (pre-collapse)
    multi_st <- mlst_key_diag %>% filter(n_ST > 1)
    write_csv(multi_st, "results/debug/mlst_multi_ST_isolates.csv")
    
    # After join, check duplicates and missing meta linkage
    dup_after <- mlst_with_meta %>% count(Isolate_ID) %>% filter(n > 1)
    write_csv(dup_after, "results/debug/join_duplicates_after.csv")
    
    miss_after <- mlst_with_meta %>%
      filter(is.na(.data[[iso_col]])) %>%        # no matching meta row
      select(Isolate_ID) %>% distinct()
    write_csv(miss_after, "results/debug/join_missing_meta.csv")
    
    message(sprintf(
      "Join OK: %d MLST rows → %d unique Isolate_ID; meta rows kept=%d; after-join dupIDs=%d, missing meta=%d",
      nrow(mlst), dplyr::n_distinct(mlst$Isolate_ID),
      dplyr::n_distinct(meta[[iso_col]]),
      nrow(dup_after), nrow(miss_after)
    ))
    
    ## ---- Audit: which row kept per Isolate_ID --------------------------------
    collapse_report <- mlst %>%
      mutate(
        ST            = as.character(ST),
        n_loci_typed  = dplyr::coalesce(.data[["n_loci_typed"]], -Inf),
        mlst_complete = dplyr::coalesce(.data[["mlst_complete"]], FALSE),
        file_pref     = dplyr::coalesce(.data[["full_path"]], .data[["file_name"]])
      ) %>%
      group_by(Isolate_ID) %>%
      arrange(desc(mlst_complete), desc(n_loci_typed), .by_group = TRUE) %>%
      summarise(
        n_candidates    = dplyr::n(),
        ST_set          = paste(sort(unique(na.omit(ST))), collapse = ","),
        n_ST            = dplyr::n_distinct(ST, na.rm = TRUE),
        chosen_ST       = dplyr::first(ST),
        chosen_complete = dplyr::first(mlst_complete),
        chosen_n_loci   = dplyr::first(n_loci_typed),
        chosen_file     = dplyr::first(file_pref),
        .groups = "drop"
      ) %>%
      mutate(
        reason = dplyr::case_when(
          n_candidates == 1      ~ "single_row",
          chosen_complete        ~ "mlst_complete_preferred",
          TRUE                   ~ "most_loci_preferred"
        )
      )
    
    readr::write_csv(collapse_report, "results/debug/mlst_collapse_choices.csv")
    message("✓ wrote collapse audit → results/debug/mlst_collapse_choices.csv")
    
    ## ---- Optional: flag multi-ST isolates in the joined table ----------------
    flag_multi <- mlst_key_diag %>% transmute(Isolate_ID, multi_ST_flag = n_ST > 1)
    mlst_with_meta <- mlst_with_meta %>%
      left_join(flag_multi, by = "Isolate_ID") %>%
      mutate(multi_ST_flag = dplyr::coalesce(multi_ST_flag, FALSE))
    write_csv(mlst_with_meta, "results/mlst/mlst_with_meta.csv")
  }
}

# End of script
#!/usr/bin/env Rscript
# =============================================================
# 08_core_vs_plasmid_STs.R – chromosomal STs vs plasmid pMLST
# Parallel version (furrr) – 9 Jul 2025
# =============================================================

suppressPackageStartupMessages({
  library(dplyr);  library(readr);  library(tidyr)
  library(purrr);   library(furrr);  library(fs)
  library(stringr); library(scales); library(processx)
})

## ---- 0 · paths -------------------------------------------------------------
fa_dir        <- "ont-yellow-routine-fastas"
core_mlst_tsv <- "results/mlst/mlst_all.tsv"

out_core_csv  <- "results/mlst/ST_core_freq.csv"
out_pmlst_long<- "results/mlst/pMLST_hits_long.csv"
out_pmlst_wide<- "results/mlst/plasmid_types_per_isolate.csv"

## ---- 1 · chromosomal ST frequencies ---------------------------------------
core_tbl <- read_tsv(core_mlst_tsv, show_col_types = FALSE)

core_tbl %>%
  count(ST, sort = TRUE) %>%
  mutate(pct = percent(n / sum(n))) %>%
  write_csv(out_core_csv)

message("✓ chromosomal ST table → ", out_core_csv)

## ---- 2 · pMLST setup -------------------------------------------------------
MLST_BIN <- Sys.which("mlst")
if (MLST_BIN == "") stop("mlst not found on PATH")

pSchemes    <- c("ac", "f", "hi1", "hi2", "i1", "n")
fasta_files <- fs::dir_ls(fa_dir, glob = "*.fasta")

run_one <- function(fasta, scheme) {
  out <- fs::path("results/mlst/raw",
                  paste0(fs::path_file(fasta), ".", scheme, ".csv"))
  if (fs::file_exists(out))                       # cached
    return(read_csv(out, col_types = cols(.default = "c")))
  
  res <- processx::run(MLST_BIN,
                       c("--quiet", "--threads", "1",
                         "--scheme", scheme, "--csv", "--legacy", fasta),
                       error_on_status = FALSE)
  
  if (res$status != 0 || !str_detect(res$stdout, ",")) return(NULL)
  
  dat <- read_csv(I(res$stdout), col_types = cols(.default = "c")) |>
    rename_with(tolower)
  write_csv(dat, out)                             # cache
  dat
}

## ---- 3 · ETA quick test ----------------------------------------------------
test_time <- system.time(run_one(fasta_files[1], pSchemes[1]))[["elapsed"]]
est_min   <- test_time * length(fasta_files) * length(pSchemes) /
  (60 * 4)                   # 4 parallel workers below
message(sprintf("~%.1f min expected on 4 workers (%.1f s/test)",
                est_min, test_time))

## ---- 4 · parallel fan-out --------------------------------------------------
plan(multisession, workers = 4)          # ← adjust cores here

grid <- tidyr::expand_grid(file = fasta_files, scheme = pSchemes)

pmlst_hits <- grid |>
  mutate(dat = furrr::future_pmap(
    list(file, scheme), run_one,
    .progress = TRUE)) |>
  unnest(dat)                            # drops NULL automatically

if (!nrow(pmlst_hits)) {
  message("No plasmid MLST alleles detected.")
  quit(status = 0)
}

## ---- 5 · write long & wide --------------------------------------------------
write_csv(pmlst_hits, out_pmlst_long)
message("✓ pMLST long table     → ", out_pmlst_long)

pmlst_hits %>%
  select(Isolate_ID = isolate_id, p_scheme = scheme, st) %>%
  pivot_wider(names_from = p_scheme,
              values_from = st,
              values_fill = NA_character_,
              names_prefix = "pST_") %>%
  arrange(Isolate_ID) %>%
  write_csv(out_pmlst_wide)

message("✓ per-isolate pST wide → ", out_pmlst_wide)

## ---- 6 · console top-5 summary --------------------------------------------
cat("\n---------  pMLST summary  ---------------------------------\n")
pmlst_hits %>%
  count(scheme, st, sort = TRUE) %>%
  group_by(scheme) %>%
  slice_max(n, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  print(n = 30)

cat("\nFinished 08_core_vs_plasmid_STs.R ✅\n")
#!/usr/bin/env Rscript
# =============================================================
# 09_inc_plasmid_network.R  –  Inc typing via ABRicate + networks
# Yellow RoUTIne  ·  9 Jul 2025
# =============================================================

# ---- configuration ------------------------------------------
FASTA_DIR   <- "ont-yellow-routine-fastas"      # same as 06_MLST.R
THREADS     <- max(1L, parallel::detectCores() - 1L)
ABR_BIN     <- Sys.which("abricate")
DB          <- "plasmidfinder"                  # Inc typing
CACHE_DIR   <- "results/abricate_plasmid"
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

# output files -------------------------------------------------
csv_long    <- "results/plasmidfinder_hits_long.csv"
csv_matrix  <- "results/plasmidfinder_presence_absence.csv"
net_pdf1    <- "results/replicon_cooccurrence.pdf"
net_pdf2    <- "results/ST_vs_replicon_network.pdf"

# ---- libraries ----------------------------------------------
suppressPackageStartupMessages({
  library(dplyr);  library(readr);  library(tidyr)
  library(stringr);library(furrr);  library(fs)
  library(igraph); library(ggraph); library(ggplot2)
})

if (ABR_BIN == "") stop("ABRicate not found on PATH")

# ---- 1 · run (cached) ---------------------------------------
run_abricate <- function(fasta) {
  out <- fs::path(CACHE_DIR, paste0(fs::path_file(fasta), ".tab"))
  if (fs::file_exists(out)) return(out)
  
  cmd <- c("--quiet", "--threads", "1", "--db", DB, fasta)
  res <- processx::run(ABR_BIN, cmd, echo = FALSE, error_on_status = FALSE)
  if (res$status != 0)
    warning("ABRicate failed on ", fasta, "; skipping") else
      write_lines(res$stdout, out)
  out
}

fasta_files <- fs::dir_ls(FASTA_DIR, glob = "*.fasta")
plan(multisession, workers = THREADS)
tab_files <- future_map_chr(fasta_files, run_abricate, .progress = TRUE)

## ---- 2 · read + tidy long ---------------------------------------------------
read_one_tab <- function(f) {
  read_tsv(f,
           comment = "",                 #  <-  DON'T drop the header!
           col_types = cols(.default = "c"),
           show_col_types = FALSE) |>
    mutate(file = fs::path_file(f))
}


hits_long <- purrr::map_dfr(tab_files, read_one_tab) |>
  mutate(
    identity = suppressWarnings(parse_number(`%IDENTITY`)),
    coverage = suppressWarnings(parse_number(`%COVERAGE`)),
    Isolate_ID = stringr::str_extract(file, "[0-9A-Za-z]+-[0-9]+")
  ) |>
  select(Isolate_ID, accession = ACCESSION, identity, coverage,
         SEQUENCE, GENE, everything())          # keep a few informative cols


write_csv(hits_long, csv_long)
message("✓ hits (long)  → ", csv_long)

## ---- 3 · presence/absence matrix -------------------------------------------
matrix <- hits_long |>
  filter(!is.na(accession)) |>
  distinct(Isolate_ID, accession) |>
  mutate(value = 1L) |>
  tidyr::pivot_wider(names_from = accession,
                     values_from = value,
                     values_fill = 0)

write_csv(matrix, csv_matrix)
message("✓ presence/absence → ", csv_matrix)

## ---- 4 · replicon co-occurrence --------------------------------------------
# frequency trim to reduce clutter (e.g., keep replicons in ≥5 isolates)
num_mat <- matrix %>% select(-Isolate_ID) %>% as.matrix()
keep    <- which(colSums(num_mat) >= 5)
if (!length(keep)) {
  message("↪  No replicons meet the frequency threshold (>=5); skipping graphs.")
  quit(status = 0)
}
mat <- num_mat[, keep, drop = FALSE]

coocc_edges <- t(mat) %*% mat            # 2) multiply once


# … the rest of the script (ggraph code, bipartite graph, etc.) stays unchanged …
# co-occur counts
diag(coocc_edges) <- 0                    # drop self loops

edges_tbl <- as.data.frame(as.table(coocc_edges)) %>%
  filter(Freq > 0) %>%
  rename(from = Var1, to = Var2, weight = Freq)

g1 <- graph_from_data_frame(edges_tbl, directed = FALSE)
pdf(net_pdf1, width = 7, height = 7)
ggraph(g1, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = .6) +
  geom_node_point(size = 4) +
  geom_node_text(aes(label = name), vjust = 1.5) +
  theme_void() +
  ggtitle("Replicon co-occurrence (Inc types)")
dev.off()
message("✓ network plot 1 → ", net_pdf1)

# ---- 5 · ST vs replicon bipartite graph ---------------------
mlst <- read_tsv("results/mlst/mlst_all.tsv", show_col_types = FALSE) %>%
  select(Isolate_ID, ST)

bipartite_edges <- hits_long %>%
  left_join(mlst, by = "Isolate_ID") %>%
  distinct(ST, accession)

g2 <- graph_from_data_frame(bipartite_edges, directed = FALSE)
V(g2)$type <- str_detect(V(g2)$name, "^\\d+$")  # TRUE for ST nodes

pdf(net_pdf2, width = 8, height = 6)
ggraph(g2, layout = "fr") +
  geom_edge_link(alpha = .4) +
  geom_node_point(aes(color = type), size = 4) +
  scale_color_manual(values = c("steelblue", "tomato"),
                     labels = c("Replicon", "ST"),
                     name = "") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_void() +
  ggtitle("Chromosomal ST vs Inc replicon network")
dev.off()
message("✓ network plot 2 → ", net_pdf2)

message("09_inc_plasmid_network.R finished ✅")
#!/usr/bin/env Rscript
# =============================================================
# 10_replicon_heatmap.R   –   Inc-type presence/absence heat-map
# Makes:   results/replicon_heatmap.pdf
# =============================================================

## ---- 0 · paths -------------------------------------------------------------
mat_file <- "results/plasmidfinder_presence_absence.csv"
mlst_file <- "results/mlst/mlst_all.tsv"       # for ST annotation
out_pdf   <- "results/replicon_heatmap.pdf"

## ---- 1 · libraries ---------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr);  library(readr);  library(tidyr)
  library(RColorBrewer)                           # colour palette
})

### prefer ComplexHeatmap; fall back to pheatmap; else stop with message --------
use_complex <- requireNamespace("ComplexHeatmap", quietly = TRUE)
use_pheat   <- requireNamespace("pheatmap", quietly = TRUE)
if (!use_complex && !use_pheat) {
  stop("Neither ComplexHeatmap nor pheatmap is installed. Install one of them to continue.")
}

## ---- 2 · load presence/absence matrix --------------------------------------
mat <- read_csv(mat_file, show_col_types = FALSE)
rownames(mat) <- mat$Isolate_ID
mat <- mat %>% select(-Isolate_ID) %>% as.matrix()

## ---- 3 · add ST side-bar ---------------------------------------------------
mlst <- read_tsv(mlst_file, show_col_types = FALSE) %>%
  select(Isolate_ID, ST)
st_vec <- mlst$ST[match(rownames(mat), mlst$Isolate_ID)]
st_colours <- setNames(
  colorRampPalette(brewer.pal(8, "Set3"))(length(unique(st_vec))),
  sort(unique(st_vec))
)

## ---- 4 · tidy up (optional filtering) --------------------------------------
# keep only replicons present in ≥5 isolates for a cleaner figure
keep_cols <- which(colSums(mat) >= 5)
heat <- mat[, keep_cols, drop = FALSE]

## ---- 5 · produce the heat-map ---------------------------------------------
if (use_complex) {
  library(ComplexHeatmap)
  ha <- HeatmapAnnotation(ST = st_vec,
                          col = list(ST = st_colours),
                          show_annotation_name = FALSE,
                          annotation_width = unit(4, "mm"))
  
  pdf(out_pdf, width = 10, height = 8)
  Heatmap(heat,
          name = "Inc\npresent",
          col = c("white", "black"),
          show_row_names = FALSE,
          top_annotation = ha,
          cluster_rows = TRUE,
          cluster_columns = TRUE)
  dev.off()
  
} else if (use_pheat) {         # ---- pheatmap fallback ----------------------
  library(pheatmap)
  ann <- data.frame(ST = st_vec, row.names = rownames(heat))
  pdf(out_pdf, width = 10, height = 8)
  pheatmap(heat,
           color = c("white", "black"),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_row = ann,
           annotation_colors = list(ST = st_colours),
           show_rownames = FALSE,
           border_color = NA)
  dev.off()
}

# Optional PNG output for convenience (if ComplexHeatmap path was used)
png_sub <- sub("\\.pdf$", ".png", out_pdf)
if (file.exists(out_pdf)) {
  try({
    gr <- grDevices::png(png_sub, width = 1200, height = 1000, res = 150)
    grDevices::dev.off()  # placeholder in case of automated conversion step
  }, silent = TRUE)
}

message("✓ wrote heat-map  →  ", out_pdf)
#!/usr/bin/env Rscript
# =============================================================
# 10_replicon_heatmap.R   –   Inc-type presence/absence heat-map
# Makes:   results/replicon_heatmap.pdf
# =============================================================

## ---- 0 · paths -------------------------------------------------------------
mat_file <- "results/plasmidfinder_presence_absence.csv"
mlst_file <- "results/mlst/mlst_all.tsv"       # for ST annotation
out_pdf   <- "results/replicon_heatmap.pdf"

## ---- 1 · libraries ---------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr);  library(readr);  library(tidyr)
  library(RColorBrewer)                           # colour palette
})

### prefer ComplexHeatmap; fall back to pheatmap; else stop with message --------
use_complex <- requireNamespace("ComplexHeatmap", quietly = TRUE)
use_pheat   <- requireNamespace("pheatmap", quietly = TRUE)
if (!use_complex && !use_pheat) {
  stop("Neither ComplexHeatmap nor pheatmap is installed. Install one of them to continue.")
}

## ---- 2 · load presence/absence matrix --------------------------------------
mat <- read_csv(mat_file, show_col_types = FALSE)
rownames(mat) <- mat$Isolate_ID
mat <- mat %>% select(-Isolate_ID) %>% as.matrix()

## ---- 3 · add ST side-bar ---------------------------------------------------
mlst <- read_tsv(mlst_file, show_col_types = FALSE) %>%
  select(Isolate_ID, ST)
st_vec <- mlst$ST[match(rownames(mat), mlst$Isolate_ID)]
st_colours <- setNames(
  colorRampPalette(brewer.pal(8, "Set3"))(length(unique(st_vec))),
  sort(unique(st_vec))
)

## ---- 4 · tidy up (optional filtering) --------------------------------------
# keep only replicons present in ≥5 isolates for a cleaner figure
keep_cols <- which(colSums(mat) >= 5)
heat <- mat[, keep_cols, drop = FALSE]

## ---- 5 · produce the heat-map ---------------------------------------------
if (use_complex) {
  library(ComplexHeatmap)
  ha <- HeatmapAnnotation(ST = st_vec,
                          col = list(ST = st_colours),
                          show_annotation_name = FALSE,
                          annotation_width = unit(4, "mm"))
  
  pdf(out_pdf, width = 10, height = 8)
  Heatmap(heat,
          name = "Inc\npresent",
          col = c("white", "black"),
          show_row_names = FALSE,
          top_annotation = ha,
          cluster_rows = TRUE,
          cluster_columns = TRUE)
  dev.off()
  
} else if (use_pheat) {         # ---- pheatmap fallback ----------------------
  library(pheatmap)
  ann <- data.frame(ST = st_vec, row.names = rownames(heat))
  pdf(out_pdf, width = 10, height = 8)
  pheatmap(heat,
           color = c("white", "black"),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_row = ann,
           annotation_colors = list(ST = st_colours),
           show_rownames = FALSE,
           border_color = NA)
  dev.off()
}

# Optional PNG output for convenience (if ComplexHeatmap path was used)
png_sub <- sub("\\.pdf$", ".png", out_pdf)
if (file.exists(out_pdf)) {
  try({
    gr <- grDevices::png(png_sub, width = 1200, height = 1000, res = 150)
    grDevices::dev.off()  # placeholder in case of automated conversion step
  }, silent = TRUE)
}

message("✓ wrote heat-map  →  ", out_pdf)
#!/usr/bin/env Rscript
# =============================================================
# 11_within_person_strain_compare.R
# -------------------------------------------------------------
# Purpose:
#   Within-participant strain comparison of E. coli across time,
#   extending the Yellow RoUTIne pipeline. Compares consecutive
#   timepoints per participant and quantifies similarity via:
#     • MLST equality (ST, majority vote per PID×TP)
#     • VF gene presence/absence Jaccard distance
#     • Plasmid replicon Jaccard distance (if mappable)
#     • ANI / SNPs from nucmer/dnadiff (if available)
#
# Inputs (expected from prior scripts):
#   - results/vf_pa_all.csv                      (from 02_gene_presence_analysis.R)
#   - results/mlst/mlst_matrix.csv               (from 06_MLST.R)
#   - results/status_map.csv                     (from 01_prepare_assembly_metadata.R; optional)
#   - results/plasmidfinder_presence_absence.csv (from 09_inc_plasmid_network.R; optional)
#   - results/nucmer/*/run_dd.report             (if nucmer trajectories ran in 02)
#
# Outputs:
#   - results/within_person/pairwise_summary.csv
#   - results/within_person/persistence_by_participant.csv
#   - results/within_person/gene_changes_long.csv
#   - results/within_person/benchmark_within_vs_between.csv
#   - results/within_person/nucmer_pairwise_metrics.csv         (if parsed)
#   - results/within_person/debug/*.csv                         (conflicts, missing, parsing)
#   - results/within_person/params.json
#   - results/plots/within_person/*.png                         (distributions & trajectories)
#
# Key assumptions (edit thresholds via CLI if needed):
#   1) Timepoints are normalized to labels like "T0","T1",…,"Uricult" (as in prior scripts).
#   2) "Persistent (same strain)" default rule:
#        - Same ST AND
#        - VF Jaccard distance <= 0.05 AND
#        - If ANI/SNPs available: ANI >= 99.90% (or SNPs <= 20)
#      Otherwise "Replacement". Missing/discordant metrics → "Uncertain".
#   3) If multiple isolates exist within the same Participant×Timepoint,
#      VF presence/absence collapses by any-present; ST uses majority vote
#      (conflicts logged to /results/within_person/debug/).
#
# CLI (optional):
#   --pids=P01,P02        (restrict to these participants)
#   --tp_rule=all         ("all" = consecutive by factor order; reserved for future)
#   --vf_jaccard_max=0.05 (threshold for persistence)
#   --ani_min=99.9        (threshold for persistence if ANI present)
#   --snp_max=20          (threshold for persistence if SNP counts present)
#
# Style & behavior:
#   • Only needed libraries; startup messages suppressed.
#   • Creates results/ and subfolders; logs major steps with message().
#   • Fails gracefully with informative errors if required inputs missing.
# =============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(purrr); library(ggplot2); library(forcats); library(jsonlite)
})

# ------------------------------ 0 · helpers & dirs ----------------------------
base_dir <- normalizePath(getwd(), mustWork = TRUE)
out_dir  <- file.path(base_dir, "results", "within_person")
plot_dir <- file.path(base_dir, "results", "plots", "within_person")
dbg_dir  <- file.path(base_dir, "results", "within_person", "debug")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dbg_dir,  recursive = TRUE, showWarnings = FALSE)

msg <- function(fmt, ...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), sprintf(fmt, ...), "\n")
  flush.console()
}

safe_write_csv <- function(x, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  ok <- tryCatch({ readr::write_csv(x, tmp, ...); TRUE }, error = function(e) { message(conditionMessage(e)); FALSE })
  if (ok && file.exists(tmp)) file.rename(tmp, path)
  if (file.exists(path)) msg("wrote: %s  (nrow=%d)", normalizePath(path), nrow(x))
}
bench_df <- readr::read_csv(file.path(out_dir, "benchmark_within_vs_between.csv"), show_col_types = FALSE)
wilcox_p <- NA_real_
if (sum(bench_df$kind == "within", na.rm = TRUE) >= 5 &&
    sum(bench_df$kind == "between", na.rm = TRUE) >= 5) {
  wilcox_p <- suppressWarnings(
    wilcox.test(bench_df$vf_jaccard[bench_df$kind == "within"],
                bench_df$vf_jaccard[bench_df$kind == "between"],
                alternative = "less")$p.value
  )
}
subtitle_txt <- if (!is.na(wilcox_p)) paste0("Wilcoxon (Within < Between) p = ", signif(wilcox_p, 3)) else "Wilcoxon p = NA"

pD <- bench_df %>%
  ggplot(aes(kind, vf_jaccard, fill = kind)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3) +
  labs(title = "Benchmark VF distances", subtitle = subtitle_txt, x = NULL, y = "VF Jaccard", fill = "")
safe_save_plot(file.path(plot_dir, "benchmark_within_vs_between.png"), pD, 6, 4.5)
safe_save_plot <- function(filename, plot, width, height, dpi = 300, units = "in") {
  try({
    if (requireNamespace("ragg", quietly = TRUE)) {
      ggplot2::ggsave(filename, plot, device = ragg::agg_png, width = width, height = height, dpi = dpi, units = units)
    } else {
      ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi, units = units)
    }
    msg("saved plot: %s", normalizePath(filename))
  }, silent = TRUE)
}

# Timepoint normalization used throughout the pipeline
tp_norm <- function(x) {
  tp_chr     <- as.character(x)
  is_uricult <- str_detect(tp_chr, regex("uricult", ignore_case = TRUE))
  tp_num     <- suppressWarnings(as.integer(str_extract(tp_chr, "\\d+")))
  tp_lab     <- dplyr::case_when(
    is_uricult     ~ "Uricult",
    !is.na(tp_num) ~ paste0("T", tp_num),
    TRUE           ~ "Unscheduled"
  )
  tp_levels <- c(paste0("T", sort(unique(tp_num[!is.na(tp_num)]))),
                 "Uricult", "Unscheduled")
  tibble::tibble(tp_lab = factor(tp_lab, levels = tp_levels),
                 tp_num = tp_num)
}

# Ensure we have tp_lab (+tp_num) without duplicating columns
ensure_tp <- function(df) {
  out <- df
  if ("tp_lab" %in% names(out)) {
    # Only compute tp_num if missing, using the current tp_lab
    if (!"tp_num" %in% names(out)) {
      tp_chr <- as.character(out$tp_lab)
      # Derive numeric part; Uricult/Unscheduled -> NA
      tp_num <- suppressWarnings(as.integer(stringr::str_extract(tp_chr, "\\d+")))
      out$tp_num <- tp_num
    }
  } else if ("Timepoint" %in% names(out)) {
    out <- bind_cols(out, tp_norm(out$Timepoint))
  } else {
    stop("Data frame lacks both 'tp_lab' and 'Timepoint' fields.")
  }
  out
}

# Jaccard distance for binary named vectors (0/1); returns NA if union = 0
jaccard_distance01 <- function(a, b) {
  a <- as.integer(a > 0); b <- as.integer(b > 0)
  inter <- sum(a & b, na.rm = TRUE)
  uni   <- sum((a | b), na.rm = TRUE)
  if (uni == 0) return(NA_real_)
  1 - (inter / uni)
}

# ------------------------------ 1 · parse CLI & thresholds --------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  hit <- args[grepl(paste0("^--", key, "="), args)]
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}

PID_FILTER <- get_arg("pids", default = NA_character_)
TP_RULE    <- get_arg("tp_rule", default = "all")
VF_JAC_MAX <- suppressWarnings(as.numeric(get_arg("vf_jaccard_max", "0.05")))
ANI_MIN    <- suppressWarnings(as.numeric(get_arg("ani_min", "99.9")))
SNP_MAX    <- suppressWarnings(as.integer(get_arg("snp_max", "20")))
if (is.na(VF_JAC_MAX)) VF_JAC_MAX <- 0.05
if (is.na(ANI_MIN))    ANI_MIN    <- 99.9
if (is.na(SNP_MAX))    SNP_MAX    <- 20L

params <- list(
  pid_filter = if (is.na(PID_FILTER)) "ALL" else PID_FILTER,
  tp_rule = TP_RULE, vf_jaccard_max = VF_JAC_MAX,
  ani_min = ANI_MIN, snp_max = SNP_MAX
)
writeLines(jsonlite::toJSON(params, pretty = TRUE, auto_unbox = TRUE),
           file.path(out_dir, "params.json"))

# ------------------------------ 2 · check required inputs ---------------------
vf_path   <- "results/vf_pa_all.csv"
mlst_path <- "results/mlst/mlst_matrix.csv"
if (!file.exists(vf_path) || !file.exists(mlst_path)) {
  if (!file.exists(vf_path))   msg("ERROR: %s not found (run 02_gene_presence_analysis.R)", vf_path)
  if (!file.exists(mlst_path)) msg("ERROR: %s not found (run 06_MLST.R)", mlst_path)
  stop("[11] Cannot proceed without required inputs.")
}
status_map_path <- "results/status_map.csv"
replicon_path   <- "results/plasmidfinder_presence_absence.csv"

# ------------------------------ 3 · load input tables -------------------------
msg("Loading VF presence/absence …")
vf <- suppressMessages(readr::read_csv(vf_path, show_col_types = FALSE))
vf <- ensure_tp(vf)
vf <- vf %>% mutate(Participant_id = as.character(Participant_id)) %>%
  arrange(Participant_id, tp_lab)

# Identify gene columns; standardize binary & collapse duplicates by any-present
gene_cols <- setdiff(names(vf), c("Participant_id", "tp_lab", "tp_num", "Timepoint"))
if (!length(gene_cols)) stop("No VF gene columns detected in vf_pa_all.csv.")
vf[gene_cols] <- lapply(vf[gene_cols], function(x) {
  x <- suppressWarnings(as.integer(x)); x[is.na(x)] <- 0L; as.integer(x > 0L)
})
vf <- vf %>%
  group_by(Participant_id, tp_lab, tp_num) %>%
  summarise(across(all_of(gene_cols), ~ as.integer(any(. > 0, na.rm = TRUE))), .groups = "drop")

msg("Loading MLST table …")
mlst <- suppressMessages(readr::read_csv(mlst_path, show_col_types = FALSE))
mlst <- ensure_tp(mlst)

# Ensure ST column exists and is character
st_col <- names(mlst)[tolower(names(mlst)) == "st"]
if (!length(st_col)) stop("MLST table lacks an 'ST' column.")
if (!"ST" %in% names(mlst)) mlst$ST <- as.character(mlst[[st_col[1]]]) else mlst$ST <- as.character(mlst$ST)

mlst_key <- mlst %>%
  mutate(Participant_id = as.character(Participant_id)) %>%
  select(Participant_id, tp_lab, ST) %>%
  filter(!is.na(Participant_id), !is.na(tp_lab)) %>%
  group_by(Participant_id, tp_lab) %>%
  summarise(
    ST_majority = {
      s <- na.omit(ST)
      if (!length(s)) NA_character_ else names(sort(table(s), decreasing = TRUE))[1]
    },
    ST_n_unique = n_distinct(ST, na.rm = TRUE),
    .groups = "drop"
  )

conflicts <- mlst_key %>% filter(ST_n_unique > 1)
if (nrow(conflicts)) {
  msg("⚠ Detected %d Participant×Timepoint with multiple ST calls; majority vote will be used, conflicts logged.", nrow(conflicts))
  safe_write_csv(conflicts, file.path(dbg_dir, "mlst_conflicts.csv"))
}

# Optional status map (used only for annotations/plots)
status_map <- NULL
if (file.exists(status_map_path)) {
  status_map <- suppressMessages(readr::read_csv(status_map_path, show_col_types = FALSE))
  status_map <- ensure_tp(status_map)
  status_map <- status_map %>% select(Participant_id, tp_lab, Infection_Status) %>% distinct()
}

# Optional replicon matrix (map to PID×TP via Isolate_ID if possible)
rep_pid_tp <- NULL
if (file.exists(replicon_path)) {
  repl <- suppressMessages(readr::read_csv(replicon_path, show_col_types = FALSE))
  if ("Isolate_ID" %in% names(repl) && "Isolate_ID" %in% names(mlst)) {
    key <- mlst %>%
      select(Isolate_ID, Participant_id, tp_lab) %>%
      distinct() %>%
      filter(!is.na(Isolate_ID), !is.na(Participant_id), !is.na(tp_lab))
    rep_join <- repl %>% mutate(Isolate_ID = as.character(Isolate_ID)) %>%
      left_join(key, by = "Isolate_ID", relationship = "many-to-one") %>%
      filter(!is.na(Participant_id), !is.na(tp_lab))
    rep_cols <- setdiff(names(rep_join), c("Isolate_ID","Participant_id","tp_lab"))
    if (length(rep_cols)) {
      rep_join[rep_cols] <- lapply(rep_join[rep_cols], function(x) { x <- suppressWarnings(as.integer(x)); x[is.na(x)] <- 0L; as.integer(x > 0L) })
      rep_pid_tp <- rep_join %>%
        group_by(Participant_id, tp_lab) %>%
        summarise(across(all_of(rep_cols), ~ as.integer(any(. > 0, na.rm = TRUE))), .groups = "drop")
      msg("✓ Replicon matrix mapped to Participant×Timepoint (n=%d rows).", nrow(rep_pid_tp))
    } else {
      msg("↪ Replicon file found but no replicon columns to summarize.")
    }
  } else {
    msg("↪ Replicon mapping skipped: missing 'Isolate_ID' in one or both tables.")
  }
}

# ------------------------------ 4 · optional nucmer/dnadiff parsing -----------
parse_dnadiff <- function(report_file) {
  L <- readLines(report_file, warn = FALSE)
  grab <- function(key, pattern = "\\d+\\.?\\d*") {
    m <- grep(key, L, value = TRUE)
    if (!length(m)) return(NA_real_)
    as.numeric(stringr::str_extract(m[1], pattern))
  }
  tibble::tibble(
    AvgIdentity = grab("AvgIdentity"),
    TotalSNPs   = grab("TotalSNPs"),
    RefLen      = grab("^TotalBasesRef"),
    QryLen      = grab("^TotalBasesQry")
  )
}
collect_nucmer <- function() {
  root <- file.path("results","nucmer")
  if (!dir.exists(root)) return(tibble())
  dirs <- list.dirs(root, recursive = TRUE, full.names = TRUE)
  rpt  <- file.path(dirs, "run_dd.report")
  rpt  <- rpt[file.exists(rpt)]
  if (!length(rpt)) return(tibble())
  purrr::map_dfr(rpt, function(f) {
    # Expected dir name: results/nucmer/PID_ASM_Tx_vs_Ty
    dn <- basename(dirname(f))
    parts <- strsplit(dn, "_")[[1]]
    if (length(parts) < 4) return(tibble())
    pid <- parts[1]
    assembler <- parts[2]
    ab   <- strsplit(paste(parts[3:length(parts)], collapse = "_"), "_vs_")[[1]]
    if (length(ab) != 2) return(tibble())
    tpA  <- ab[1]; tpB <- ab[2]
    parse_dnadiff(f) %>%
      mutate(Participant_id = pid, assembler = assembler, tp_A = tpA, tp_B = tpB,
             report_path = f)
  })
}
nucmer_tbl <- collect_nucmer()
if (nrow(nucmer_tbl)) {
  safe_write_csv(nucmer_tbl, file.path(out_dir, "nucmer_pairwise_metrics.csv"))
  msg("✓ Parsed %d nucmer/dnadiff reports.", nrow(nucmer_tbl))
} else {
  msg("↪ No nucmer/dnadiff reports found; ANI/SNP metrics will be unavailable.")
}

# ------------------------------ 5 · restrict PIDs if requested ----------------
if (!is.na(PID_FILTER)) {
  keep <- unlist(strsplit(PID_FILTER, "\\s*,\\s*"))
  vf        <- vf        %>% filter(Participant_id %in% keep)
  mlst_key  <- mlst_key  %>% filter(Participant_id %in% keep)
  if (!is.null(status_map)) status_map <- status_map %>% filter(Participant_id %in% keep)
  if (!is.null(rep_pid_tp)) rep_pid_tp <- rep_pid_tp %>% filter(Participant_id %in% keep)
  if (nrow(vf) == 0) stop("After --pids filter, no rows remain.")
}

# ------------------------------ 6 · build sample key table --------------------
# NOTE: avoid duplicating tp_lab — we already have it; only ensure tp_num.
tp_info <- vf %>%
  select(Participant_id, tp_lab, tp_num) %>%
  mutate(tp_num = suppressWarnings(as.integer(tp_num))) %>%
  distinct() %>%
  arrange(Participant_id, as.integer(tp_num), tp_lab)

samples <- tp_info %>%
  left_join(mlst_key %>% select(Participant_id, tp_lab, ST = ST_majority),
            by = c("Participant_id","tp_lab"))

missing_st <- samples %>% filter(is.na(ST))
if (nrow(missing_st)) safe_write_csv(missing_st, file.path(dbg_dir, "missing_ST_pid_tp.csv"))

# ------------------------------ 7 · build consecutive pairs -------------------
pairs <- samples %>%
  group_by(Participant_id) %>%
  arrange(Participant_id,
          # Order: numeric T's by tp_num; then Uricult; then Unscheduled
          is.na(tp_num), tp_num, factor(tp_lab, levels = levels(tp_lab))) %>%
  mutate(tp_next = lead(tp_lab), tp_num_next = lead(tp_num), ST_next = lead(ST)) %>%
  ungroup() %>%
  filter(!is.na(tp_next)) %>%
  transmute(Participant_id,
            tp_A = tp_lab,    tp_B = tp_next,
            tpA_num = tp_num, tpB_num = tp_num_next,
            ST_A = ST, ST_B = ST_next)

if (!nrow(pairs)) {
  msg("No consecutive pairs available. Provide participants with ≥2 timepoints.")
  safe_write_csv(samples, file.path(out_dir, "sample_key.csv"))
  quit(save = "no", status = 0)
}

# ------------------------------ 8 · compute VF Jaccard & gene changes ---------
vf_mat <- vf %>% select(Participant_id, tp_lab, all_of(gene_cols))

get_vf_vec <- function(pid, tp) {
  row <- vf_mat %>% filter(Participant_id == pid, tp_lab == tp)
  if (nrow(row) != 1L) return(setNames(rep(NA_integer_, length(gene_cols)), gene_cols))
  v <- as.integer(row[1, gene_cols, drop = TRUE])
  names(v) <- gene_cols
  v
}

pairs <- pairs %>%
  mutate(
    vfA_vec = purrr::map2(Participant_id, tp_A, \(pid, tp) get_vf_vec(pid, tp)),
    vfB_vec = purrr::map2(Participant_id, tp_B, \(pid, tp) get_vf_vec(pid, tp)),
    vf_jaccard = purrr::map2_dbl(vfA_vec, vfB_vec, jaccard_distance01),
    vf_intersect = purrr::map2(vfA_vec, vfB_vec, \(a,b) names(which(a == 1L & b == 1L))),
    vf_gain     = purrr::map2(vfA_vec, vfB_vec, \(a,b) names(which(a == 0L & b == 1L))),
    vf_loss     = purrr::map2(vfA_vec, vfB_vec, \(a,b) names(which(a == 1L & b == 0L)))
  )

gene_changes_long <- pairs %>%
  select(Participant_id, tp_A, tp_B, vf_gain, vf_loss) %>%
  tidyr::pivot_longer(c(vf_gain, vf_loss), names_to = "change_type", values_to = "GENE") %>%
  tidyr::unnest(GENE) %>%
  mutate(change_type = dplyr::recode(change_type, vf_gain = "gained_in_B", vf_loss = "lost_in_B"))

safe_write_csv(gene_changes_long, file.path(out_dir, "gene_changes_long.csv"))

# ------------------------------ 9 · plasmid replicon Jaccard (optional) -------
rep_cols <- NULL
if (!is.null(rep_pid_tp)) {
  rep_cols <- setdiff(names(rep_pid_tp), c("Participant_id","tp_lab"))
  get_rep_vec <- function(pid, tp) {
    row <- rep_pid_tp %>% filter(Participant_id == pid, tp_lab == tp)
    if (nrow(row) != 1L) return(rep(NA_integer_, length(rep_cols)))
    as.integer(row[1, rep_cols, drop = TRUE])
  }
  pairs <- pairs %>%
    mutate(
      repA_vec = purrr::map2(Participant_id, tp_A, \(pid, tp) get_rep_vec(pid, tp)),
      repB_vec = purrr::map2(Participant_id, tp_B, \(pid, tp) get_rep_vec(pid, tp)),
      replicon_jaccard = if (length(rep_cols)) purrr::map2_dbl(repA_vec, repB_vec, jaccard_distance01) else NA_real_
    )
} else {
  pairs$replicon_jaccard <- NA_real_
}

# ------------------------------ 10 · merge nucmer/dnadiff metrics (optional) ---
if (exists("nucmer_tbl") && nrow(nucmer_tbl)) {
  nucmer_summ <- nucmer_tbl %>%
    group_by(Participant_id, tp_A, tp_B) %>%
    summarise(
      ANI  = suppressWarnings(max(AvgIdentity, na.rm = TRUE)),
      SNPs = suppressWarnings(min(TotalSNPs, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      ANI  = ifelse(is.infinite(ANI),  NA_real_, ANI),
      SNPs = ifelse(is.infinite(SNPs), NA_real_, SNPs)
    )
  pairs <- pairs %>%
    left_join(nucmer_summ, by = c("Participant_id","tp_A","tp_B"))
} else {
  pairs$ANI  <- NA_real_
  pairs$SNPs <- NA_real_
}

# ------------------------------ 11 · classify persistence vs replacement -------
pairs <- pairs %>%
  mutate(
    ST_same = !is.na(ST_A) & !is.na(ST_B) & ST_A == ST_B,
    vf_ok   = !is.na(vf_jaccard) & vf_jaccard <= VF_JAC_MAX,
    ani_ok  = is.na(ANI)  | (ANI  >= ANI_MIN),
    snp_ok  = is.na(SNPs) | (SNPs <= SNP_MAX),
    call    = dplyr::case_when(
      ST_same & vf_ok & ani_ok & snp_ok ~ "Persistent",
      !ST_same                           ~ "Replacement",
      ST_same & !vf_ok                   ~ "Replacement",
      ST_same & vf_ok & (!ani_ok | !snp_ok) ~ "Replacement",
      TRUE                               ~ "Uncertain"
    )
  )

# ------------------------------ 12 · participant-level summaries ---------------
persistence_by_participant <- pairs %>%
  group_by(Participant_id) %>%
  summarise(
    n_pairs = n(),
    n_persistent = sum(call == "Persistent", na.rm = TRUE),
    n_replacement = sum(call == "Replacement", na.rm = TRUE),
    n_uncertain = sum(call == "Uncertain", na.rm = TRUE),
    prop_persistent = ifelse(n_pairs > 0, n_persistent / n_pairs, NA_real_),
    .groups = "drop"
  )

# ------------------------------ 13 · write pairwise summary --------------------
pairwise_summary <- pairs %>%
  select(Participant_id, tp_A, tp_B, tpA_num, tpB_num,
         ST_A, ST_B, ST_same,
         vf_jaccard, replicon_jaccard, ANI, SNPs, call)

safe_write_csv(pairwise_summary, file.path(out_dir, "pairwise_summary.csv"))
safe_write_csv(persistence_by_participant, file.path(out_dir, "persistence_by_participant.csv"))

# ------------------------------ 14 · benchmark within vs between ---------------
# Build a between-participant baseline using VF Jaccard for pairs at the same
# timepoint. We compute explicit pairwise combinations to keep it stable.

vf_for_bench <- vf_mat %>% select(Participant_id, tp_lab, all_of(gene_cols))

# helper: compute all unordered pairwise Jaccard distances within a single tp_lab
.compute_between_for_tp <- function(df_tp) {
  # df_tp has columns: Participant_id, tp_lab, <genes...>
  if (nrow(df_tp) < 2) return(tibble::tibble())
  pids <- df_tp$Participant_id
  tpv  <- unique(as.character(df_tp$tp_lab))
  tpv  <- if (length(tpv)) tpv[1] else NA_character_
  
  mat  <- as.matrix(df_tp[, gene_cols, drop = FALSE])
  idx  <- t(combn(seq_len(nrow(mat)), 2))
  
  jd <- apply(idx, 1, function(ix) {
    jaccard_distance01(mat[ix[1], ], mat[ix[2], ])
  })
  
  tibble::tibble(
    kind       = "between",
    tp_lab     = tpv,
    pid_A      = pids[idx[, 1]],
    pid_B      = pids[idx[, 2]],
    vf_jaccard = as.numeric(jd)
  ) %>% filter(!is.na(vf_jaccard))
}

within_dist <- pairwise_summary %>%
  transmute(kind = "within", vf_jaccard = vf_jaccard) %>%
  filter(!is.na(vf_jaccard))

between_dist_all <- vf_for_bench %>%
  group_by(tp_lab) %>% group_split() %>%
  purrr::map_dfr(.compute_between_for_tp)

# size-match the between distribution to the number of within pairs (if needed)
set.seed(1)
n_within <- nrow(within_dist)
between_dist <- if (nrow(between_dist_all) > n_within) {
  dplyr::slice_sample(between_dist_all, n = n_within) %>% transmute(kind, vf_jaccard)
} else {
  between_dist_all %>% transmute(kind, vf_jaccard)
}

benchmark_within_vs_between <- dplyr::bind_rows(within_dist, between_dist)
safe_write_csv(benchmark_within_vs_between,
               file.path(out_dir, "benchmark_within_vs_between.csv"))

# ------------------------------ 15 · plots: VF Jaccard distributions -----------
p_vf_jacc <- ggplot(benchmark_within_vs_between,
                    aes(x = vf_jaccard, fill = kind)) +
  geom_density(alpha = 0.45, adjust = 1.2) +
  theme_bw(base_size = 12) +
  labs(title = "VF Jaccard distance: within-person vs between-person",
       x = "Jaccard distance (0 = identical, 1 = disjoint)",
       y = "Density", fill = "") +
  guides(fill = guide_legend(override.aes = list(alpha = 0.6)))

# rename existing density save (change your earlier save line)
safe_save_plot(file.path(plot_dir, "vf_jaccard_within_vs_between_density.png"),
               p_vf_jacc, width = 7, height = 4.5, dpi = 300)

## Boxplot/jitter for within vs between with Wilcoxon p-value (uses 'kind' & 'vf_jaccard')
bench_df <- readr::read_csv(file.path(out_dir, "benchmark_within_vs_between.csv"), show_col_types = FALSE)

wilcox_p <- NA_real_
if (sum(bench_df$kind == "within", na.rm = TRUE) >= 5 &&
    sum(bench_df$kind == "between", na.rm = TRUE) >= 5) {
  wilcox_p <- suppressWarnings(
    wilcox.test(bench_df$vf_jaccard[bench_df$kind == "within"],
                bench_df$vf_jaccard[bench_df$kind == "between"],
                alternative = "less")$p.value
  )
}
subtitle_txt <- if (!is.na(wilcox_p)) paste0("Wilcoxon (Within < Between) p = ", signif(wilcox_p, 3)) else "Wilcoxon p = NA"

pD <- bench_df %>%
  ggplot(aes(kind, vf_jaccard, fill = kind)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3) +
  labs(title = "Benchmark VF distances", subtitle = subtitle_txt, x = NULL, y = "VF Jaccard", fill = "")
safe_save_plot(file.path(plot_dir, "benchmark_within_vs_between.png"), pD, 6, 4.5)
# (optional) drop a text file like the initial script did
writeLines(paste0("wilcox_p (Within < Between) = ",
                  ifelse(is.na(wilcox_p), "NA", signif(wilcox_p, 6))),
           file.path(out_dir, "benchmark_wilcoxon.txt"))

# ------------------------------ 16 · plots: call summary -----------------------
p_call <- pairs %>%
  count(call, name = "n") %>%
  mutate(call = factor(call, levels = c("Persistent","Replacement","Uncertain"))) %>%
  ggplot(aes(call, n)) + geom_col() + theme_bw(base_size = 12) +
  labs(title = "Pairwise calls", x = NULL, y = "Count")

safe_save_plot(file.path(plot_dir, "pairwise_calls_bar.png"),
               p_call, width = 5, height = 4, dpi = 300)

# ------------------------------ 17 · optional: scatter VF vs ANI ---------------
if (any(!is.na(pairs$ANI))) {
  p_scatter <- pairs %>%
    ggplot(aes(x = vf_jaccard, y = ANI, shape = call)) +
    geom_point(size = 2) +
    theme_bw(base_size = 12) +
    labs(title = "Within-person pairs: VF Jaccard vs ANI",
         x = "VF Jaccard distance",
         y = "ANI (%)",
         shape = "Call")
  safe_save_plot(file.path(plot_dir, "vf_vs_ani_scatter.png"),
                 p_scatter, width = 6, height = 4.5, dpi = 300)
}

# A) Distribution of VF Jaccard by classification (matches initial script idea)
pA <- pairs %>%
  filter(!is.na(vf_jaccard)) %>%
  ggplot(aes(vf_jaccard, fill = call)) +
  geom_histogram(bins = 30, alpha = 0.8, position = "identity") +
  geom_vline(xintercept = VF_JAC_MAX, linetype = "dashed") +
  labs(title = "Within-person VF Jaccard distances",
       subtitle = paste0("Dashed line = threshold ", VF_JAC_MAX),
       x = "VF Jaccard distance (0 = identical, 1 = disjoint)",
       y = "Pairs", fill = "Call")
safe_save_plot(file.path(plot_dir, "vf_jaccard_by_classification.png"), pA, 8, 5)

# B) Trajectories of VF distance (consecutive TPs) — requires tp_num per TP
# Build a reliable PID×TP -> tp_num lookup (avoid relying on 'samples')
tp_lookup <- vf %>%
  dplyr::select(Participant_id, tp_lab, tp_num) %>%
  dplyr::distinct()

# B) Trajectories of VF distance (consecutive TPs)
traj_df <- pairs %>%
  dplyr::left_join(tp_lookup, by = c("Participant_id","tp_A" = "tp_lab")) %>%
  dplyr::rename(tp_num_A = tp_num) %>%
  dplyr::left_join(tp_lookup, by = c("Participant_id","tp_B" = "tp_lab")) %>%
  dplyr::rename(tp_num_B = tp_num) %>%
  dplyr::filter(!is.na(tp_num_A), !is.na(tp_num_B))

pC <- traj_df %>%
  ggplot2::ggplot(ggplot2::aes(tp_num_A, vf_jaccard, group = Participant_id, color = call)) +
  ggplot2::geom_line(alpha = 0.4) +
  ggplot2::geom_point() +
  ggplot2::scale_x_continuous(breaks = sort(unique(traj_df$tp_num_A))) +
  ggplot2::labs(title = "Trajectories of within-person VF distance (consecutive TPs)",
                x = "Numeric timepoint (A)", y = "VF Jaccard")
safe_save_plot(file.path(plot_dir, "vf_trajectory_consecutive.png"), pC, 8, 5)


# ------------------------------ 18 · write extra debug aids --------------------
# (These are handy when diagnosing mismatches or empty gene change tables.)
safe_write_csv(samples, file.path(out_dir, "sample_key.csv"))
safe_write_csv(pairs,    file.path(out_dir, "pairs_with_metrics.csv"))

msg("All done. Saved benchmark table and plots in results/within_person/ and results/plots/within_person/")

# Tiny README
readme_lines <- c(
  "Within-person strain comparison (Yellow RoUTIne):",
  paste0("- Persistence rule: ST_same & VF_Jaccard ≤ ", VF_JAC_MAX,
         if (exists("nucmer_tbl") && nrow(nucmer_tbl))
           paste0(" & ANI ≥ ", ANI_MIN, "% (or SNPs ≤ ", SNP_MAX, ")")
         else ""),
  "- Missing core metrics or discordant thresholds → 'Uncertain'.",
  "- See pairwise_summary.csv for per-pair metrics; debug/ for conflicts/missing data.",
  "- Optional: replicon Jaccard is computed only if Isolate_ID→Participant×Timepoint mapping is available."
)
writeLines(readme_lines, file.path(out_dir, "README_11_within_person.txt"))

# Verification table
verify_file <- function(p) file.exists(p) && isTRUE(file.size(p) > 0)
must_have <- c(file.path(out_dir, "pairwise_summary.csv"),
               file.path(out_dir, "persistence_by_participant.csv"),
               file.path(out_dir, "gene_changes_long.csv"))
ok_vec <- vapply(must_have, verify_file, logical(1))
verif_tbl <- tibble::tibble(file = must_have, exists_nonempty = ok_vec)
safe_write_csv(verif_tbl, file.path(out_dir, "verification_table.csv"))

# Session info snapshot
sess_txt <- capture.output(sessionInfo())
writeLines(sess_txt, file.path(out_dir, "sessionInfo.txt"))
#!/usr/bin/env Rscript
# =============================================================
# 11_within_person_strain_compare.R
# -------------------------------------------------------------
# Purpose:
#   Within-participant strain comparison of E. coli across time,
#   extending the Yellow RoUTIne pipeline. Compares consecutive
#   timepoints per participant and quantifies similarity via:
#     • MLST equality (ST, majority vote per PID×TP)
#     • VF gene presence/absence Jaccard distance
#     • Plasmid replicon Jaccard distance (if mappable)
#     • ANI / SNPs from nucmer/dnadiff (if available)
#
# Inputs (expected from prior scripts):
#   - results/vf_pa_all.csv                      (from 02_gene_presence_analysis.R)
#   - results/mlst/mlst_matrix.csv               (from 06_MLST.R)
#   - results/status_map.csv                     (from 01_prepare_assembly_metadata.R; optional)
#   - results/plasmidfinder_presence_absence.csv (from 09_inc_plasmid_network.R; optional)
#   - results/nucmer/*/run_dd.report             (if nucmer trajectories ran in 02)
#
# Outputs:
#   - results/within_person/pairwise_summary.csv
#   - results/within_person/persistence_by_participant.csv
#   - results/within_person/gene_changes_long.csv
#   - results/within_person/benchmark_within_vs_between.csv
#   - results/within_person/nucmer_pairwise_metrics.csv         (if parsed)
#   - results/within_person/debug/*.csv                         (conflicts, missing, parsing)
#   - results/within_person/params.json
#   - results/plots/within_person/*.png                         (distributions & trajectories)
#
# Key assumptions (edit thresholds via CLI if needed):
#   1) Timepoints are normalized to labels like "T0","T1",…,"Uricult" (as in prior scripts).
#   2) "Persistent (same strain)" default rule:
#        - Same ST AND
#        - VF Jaccard distance <= 0.05 AND
#        - If ANI/SNPs available: ANI >= 99.90% (or SNPs <= 20)
#      Otherwise "Replacement". Missing/discordant metrics → "Uncertain".
#   3) If multiple isolates exist within the same Participant×Timepoint,
#      VF presence/absence collapses by any-present; ST uses majority vote
#      (conflicts logged to /results/within_person/debug/).
#
# CLI (optional):
#   --pids=P01,P02        (restrict to these participants)
#   --tp_rule=all         ("all" = consecutive by factor order; reserved for future)
#   --vf_jaccard_max=0.05 (threshold for persistence)
#   --ani_min=99.9        (threshold for persistence if ANI present)
#   --snp_max=20          (threshold for persistence if SNP counts present)
#
# Style & behavior:
#   • Only needed libraries; startup messages suppressed.
#   • Creates results/ and subfolders; logs major steps with message().
#   • Fails gracefully with informative errors if required inputs missing.
# =============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(purrr); library(ggplot2); library(forcats); library(jsonlite)
})

# ------------------------------ 0 · helpers & dirs ----------------------------
base_dir <- normalizePath(getwd(), mustWork = TRUE)
out_dir  <- file.path(base_dir, "results", "within_person")
plot_dir <- file.path(base_dir, "results", "plots", "within_person")
dbg_dir  <- file.path(base_dir, "results", "within_person", "debug")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dbg_dir,  recursive = TRUE, showWarnings = FALSE)

msg <- function(fmt, ...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), sprintf(fmt, ...), "\n")
  flush.console()
}

safe_write_csv <- function(x, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  ok <- tryCatch({ readr::write_csv(x, tmp, ...); TRUE }, error = function(e) { message(conditionMessage(e)); FALSE })
  if (ok && file.exists(tmp)) file.rename(tmp, path)
  if (file.exists(path)) msg("wrote: %s  (nrow=%d)", normalizePath(path), nrow(x))
}
bench_df <- readr::read_csv(file.path(out_dir, "benchmark_within_vs_between.csv"), show_col_types = FALSE)
wilcox_p <- NA_real_
if (sum(bench_df$kind == "within", na.rm = TRUE) >= 5 &&
    sum(bench_df$kind == "between", na.rm = TRUE) >= 5) {
  wilcox_p <- suppressWarnings(
    wilcox.test(bench_df$vf_jaccard[bench_df$kind == "within"],
                bench_df$vf_jaccard[bench_df$kind == "between"],
                alternative = "less")$p.value
  )
}
subtitle_txt <- if (!is.na(wilcox_p)) paste0("Wilcoxon (Within < Between) p = ", signif(wilcox_p, 3)) else "Wilcoxon p = NA"

pD <- bench_df %>%
  ggplot(aes(kind, vf_jaccard, fill = kind)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3) +
  labs(title = "Benchmark VF distances", subtitle = subtitle_txt, x = NULL, y = "VF Jaccard", fill = "")
safe_save_plot(file.path(plot_dir, "benchmark_within_vs_between.png"), pD, 6, 4.5)
safe_save_plot <- function(filename, plot, width, height, dpi = 300, units = "in") {
  try({
    if (requireNamespace("ragg", quietly = TRUE)) {
      ggplot2::ggsave(filename, plot, device = ragg::agg_png, width = width, height = height, dpi = dpi, units = units)
    } else {
      ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi, units = units)
    }
    msg("saved plot: %s", normalizePath(filename))
  }, silent = TRUE)
}

# Timepoint normalization used throughout the pipeline
tp_norm <- function(x) {
  tp_chr     <- as.character(x)
  is_uricult <- str_detect(tp_chr, regex("uricult", ignore_case = TRUE))
  tp_num     <- suppressWarnings(as.integer(str_extract(tp_chr, "\\d+")))
  tp_lab     <- dplyr::case_when(
    is_uricult     ~ "Uricult",
    !is.na(tp_num) ~ paste0("T", tp_num),
    TRUE           ~ "Unscheduled"
  )
  tp_levels <- c(paste0("T", sort(unique(tp_num[!is.na(tp_num)]))),
                 "Uricult", "Unscheduled")
  tibble::tibble(tp_lab = factor(tp_lab, levels = tp_levels),
                 tp_num = tp_num)
}

# Ensure we have tp_lab (+tp_num) without duplicating columns
ensure_tp <- function(df) {
  out <- df
  if ("tp_lab" %in% names(out)) {
    # Only compute tp_num if missing, using the current tp_lab
    if (!"tp_num" %in% names(out)) {
      tp_chr <- as.character(out$tp_lab)
      # Derive numeric part; Uricult/Unscheduled -> NA
      tp_num <- suppressWarnings(as.integer(stringr::str_extract(tp_chr, "\\d+")))
      out$tp_num <- tp_num
    }
  } else if ("Timepoint" %in% names(out)) {
    out <- bind_cols(out, tp_norm(out$Timepoint))
  } else {
    stop("Data frame lacks both 'tp_lab' and 'Timepoint' fields.")
  }
  out
}

# Jaccard distance for binary named vectors (0/1); returns NA if union = 0
jaccard_distance01 <- function(a, b) {
  a <- as.integer(a > 0); b <- as.integer(b > 0)
  inter <- sum(a & b, na.rm = TRUE)
  uni   <- sum((a | b), na.rm = TRUE)
  if (uni == 0) return(NA_real_)
  1 - (inter / uni)
}

# ------------------------------ 1 · parse CLI & thresholds --------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  hit <- args[grepl(paste0("^--", key, "="), args)]
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}

PID_FILTER <- get_arg("pids", default = NA_character_)
TP_RULE    <- get_arg("tp_rule", default = "all")
VF_JAC_MAX <- suppressWarnings(as.numeric(get_arg("vf_jaccard_max", "0.05")))
ANI_MIN    <- suppressWarnings(as.numeric(get_arg("ani_min", "99.9")))
SNP_MAX    <- suppressWarnings(as.integer(get_arg("snp_max", "20")))
if (is.na(VF_JAC_MAX)) VF_JAC_MAX <- 0.05
if (is.na(ANI_MIN))    ANI_MIN    <- 99.9
if (is.na(SNP_MAX))    SNP_MAX    <- 20L

params <- list(
  pid_filter = if (is.na(PID_FILTER)) "ALL" else PID_FILTER,
  tp_rule = TP_RULE, vf_jaccard_max = VF_JAC_MAX,
  ani_min = ANI_MIN, snp_max = SNP_MAX
)
writeLines(jsonlite::toJSON(params, pretty = TRUE, auto_unbox = TRUE),
           file.path(out_dir, "params.json"))

# ------------------------------ 2 · check required inputs ---------------------
vf_path   <- "results/vf_pa_all.csv"
mlst_path <- "results/mlst/mlst_matrix.csv"
if (!file.exists(vf_path) || !file.exists(mlst_path)) {
  if (!file.exists(vf_path))   msg("ERROR: %s not found (run 02_gene_presence_analysis.R)", vf_path)
  if (!file.exists(mlst_path)) msg("ERROR: %s not found (run 06_MLST.R)", mlst_path)
  stop("[11] Cannot proceed without required inputs.")
}
status_map_path <- "results/status_map.csv"
replicon_path   <- "results/plasmidfinder_presence_absence.csv"

# ------------------------------ 3 · load input tables -------------------------
msg("Loading VF presence/absence …")
vf <- suppressMessages(readr::read_csv(vf_path, show_col_types = FALSE))
vf <- ensure_tp(vf)
vf <- vf %>% mutate(Participant_id = as.character(Participant_id)) %>%
  arrange(Participant_id, tp_lab)

# Identify gene columns; standardize binary & collapse duplicates by any-present
gene_cols <- setdiff(names(vf), c("Participant_id", "tp_lab", "tp_num", "Timepoint"))
if (!length(gene_cols)) stop("No VF gene columns detected in vf_pa_all.csv.")
vf[gene_cols] <- lapply(vf[gene_cols], function(x) {
  x <- suppressWarnings(as.integer(x)); x[is.na(x)] <- 0L; as.integer(x > 0L)
})
vf <- vf %>%
  group_by(Participant_id, tp_lab, tp_num) %>%
  summarise(across(all_of(gene_cols), ~ as.integer(any(. > 0, na.rm = TRUE))), .groups = "drop")

msg("Loading MLST table …")
mlst <- suppressMessages(readr::read_csv(mlst_path, show_col_types = FALSE))
mlst <- ensure_tp(mlst)

# Ensure ST column exists and is character
st_col <- names(mlst)[tolower(names(mlst)) == "st"]
if (!length(st_col)) stop("MLST table lacks an 'ST' column.")
if (!"ST" %in% names(mlst)) mlst$ST <- as.character(mlst[[st_col[1]]]) else mlst$ST <- as.character(mlst$ST)

mlst_key <- mlst %>%
  mutate(Participant_id = as.character(Participant_id)) %>%
  select(Participant_id, tp_lab, ST) %>%
  filter(!is.na(Participant_id), !is.na(tp_lab)) %>%
  group_by(Participant_id, tp_lab) %>%
  summarise(
    ST_majority = {
      s <- na.omit(ST)
      if (!length(s)) NA_character_ else names(sort(table(s), decreasing = TRUE))[1]
    },
    ST_n_unique = n_distinct(ST, na.rm = TRUE),
    .groups = "drop"
  )

conflicts <- mlst_key %>% filter(ST_n_unique > 1)
if (nrow(conflicts)) {
  msg("⚠ Detected %d Participant×Timepoint with multiple ST calls; majority vote will be used, conflicts logged.", nrow(conflicts))
  safe_write_csv(conflicts, file.path(dbg_dir, "mlst_conflicts.csv"))
}

# Optional status map (used only for annotations/plots)
status_map <- NULL
if (file.exists(status_map_path)) {
  status_map <- suppressMessages(readr::read_csv(status_map_path, show_col_types = FALSE))
  status_map <- ensure_tp(status_map)
  status_map <- status_map %>% select(Participant_id, tp_lab, Infection_Status) %>% distinct()
}

# Optional replicon matrix (map to PID×TP via Isolate_ID if possible)
rep_pid_tp <- NULL
if (file.exists(replicon_path)) {
  repl <- suppressMessages(readr::read_csv(replicon_path, show_col_types = FALSE))
  if ("Isolate_ID" %in% names(repl) && "Isolate_ID" %in% names(mlst)) {
    key <- mlst %>%
      select(Isolate_ID, Participant_id, tp_lab) %>%
      distinct() %>%
      filter(!is.na(Isolate_ID), !is.na(Participant_id), !is.na(tp_lab))
    rep_join <- repl %>% mutate(Isolate_ID = as.character(Isolate_ID)) %>%
      left_join(key, by = "Isolate_ID", relationship = "many-to-one") %>%
      filter(!is.na(Participant_id), !is.na(tp_lab))
    rep_cols <- setdiff(names(rep_join), c("Isolate_ID","Participant_id","tp_lab"))
    if (length(rep_cols)) {
      rep_join[rep_cols] <- lapply(rep_join[rep_cols], function(x) { x <- suppressWarnings(as.integer(x)); x[is.na(x)] <- 0L; as.integer(x > 0L) })
      rep_pid_tp <- rep_join %>%
        group_by(Participant_id, tp_lab) %>%
        summarise(across(all_of(rep_cols), ~ as.integer(any(. > 0, na.rm = TRUE))), .groups = "drop")
      msg("✓ Replicon matrix mapped to Participant×Timepoint (n=%d rows).", nrow(rep_pid_tp))
    } else {
      msg("↪ Replicon file found but no replicon columns to summarize.")
    }
  } else {
    msg("↪ Replicon mapping skipped: missing 'Isolate_ID' in one or both tables.")
  }
}

# ------------------------------ 4 · optional nucmer/dnadiff parsing -----------
parse_dnadiff <- function(report_file) {
  L <- readLines(report_file, warn = FALSE)
  grab <- function(key, pattern = "\\d+\\.?\\d*") {
    m <- grep(key, L, value = TRUE)
    if (!length(m)) return(NA_real_)
    as.numeric(stringr::str_extract(m[1], pattern))
  }
  tibble::tibble(
    AvgIdentity = grab("AvgIdentity"),
    TotalSNPs   = grab("TotalSNPs"),
    RefLen      = grab("^TotalBasesRef"),
    QryLen      = grab("^TotalBasesQry")
  )
}
collect_nucmer <- function() {
  root <- file.path("results","nucmer")
  if (!dir.exists(root)) return(tibble())
  dirs <- list.dirs(root, recursive = TRUE, full.names = TRUE)
  rpt  <- file.path(dirs, "run_dd.report")
  rpt  <- rpt[file.exists(rpt)]
  if (!length(rpt)) return(tibble())
  purrr::map_dfr(rpt, function(f) {
    # Expected dir name: results/nucmer/PID_ASM_Tx_vs_Ty
    dn <- basename(dirname(f))
    parts <- strsplit(dn, "_")[[1]]
    if (length(parts) < 4) return(tibble())
    pid <- parts[1]
    assembler <- parts[2]
    ab   <- strsplit(paste(parts[3:length(parts)], collapse = "_"), "_vs_")[[1]]
    if (length(ab) != 2) return(tibble())
    tpA  <- ab[1]; tpB <- ab[2]
    parse_dnadiff(f) %>%
      mutate(Participant_id = pid, assembler = assembler, tp_A = tpA, tp_B = tpB,
             report_path = f)
  })
}
nucmer_tbl <- collect_nucmer()
if (nrow(nucmer_tbl)) {
  safe_write_csv(nucmer_tbl, file.path(out_dir, "nucmer_pairwise_metrics.csv"))
  msg("✓ Parsed %d nucmer/dnadiff reports.", nrow(nucmer_tbl))
} else {
  msg("↪ No nucmer/dnadiff reports found; ANI/SNP metrics will be unavailable.")
}

# ------------------------------ 5 · restrict PIDs if requested ----------------
if (!is.na(PID_FILTER)) {
  keep <- unlist(strsplit(PID_FILTER, "\\s*,\\s*"))
  vf        <- vf        %>% filter(Participant_id %in% keep)
  mlst_key  <- mlst_key  %>% filter(Participant_id %in% keep)
  if (!is.null(status_map)) status_map <- status_map %>% filter(Participant_id %in% keep)
  if (!is.null(rep_pid_tp)) rep_pid_tp <- rep_pid_tp %>% filter(Participant_id %in% keep)
  if (nrow(vf) == 0) stop("After --pids filter, no rows remain.")
}

# ------------------------------ 6 · build sample key table --------------------
# NOTE: avoid duplicating tp_lab — we already have it; only ensure tp_num.
tp_info <- vf %>%
  select(Participant_id, tp_lab, tp_num) %>%
  mutate(tp_num = suppressWarnings(as.integer(tp_num))) %>%
  distinct() %>%
  arrange(Participant_id, as.integer(tp_num), tp_lab)

samples <- tp_info %>%
  left_join(mlst_key %>% select(Participant_id, tp_lab, ST = ST_majority),
            by = c("Participant_id","tp_lab"))

missing_st <- samples %>% filter(is.na(ST))
if (nrow(missing_st)) safe_write_csv(missing_st, file.path(dbg_dir, "missing_ST_pid_tp.csv"))

# ------------------------------ 7 · build consecutive pairs -------------------
pairs <- samples %>%
  group_by(Participant_id) %>%
  arrange(Participant_id,
          # Order: numeric T's by tp_num; then Uricult; then Unscheduled
          is.na(tp_num), tp_num, factor(tp_lab, levels = levels(tp_lab))) %>%
  mutate(tp_next = lead(tp_lab), tp_num_next = lead(tp_num), ST_next = lead(ST)) %>%
  ungroup() %>%
  filter(!is.na(tp_next)) %>%
  transmute(Participant_id,
            tp_A = tp_lab,    tp_B = tp_next,
            tpA_num = tp_num, tpB_num = tp_num_next,
            ST_A = ST, ST_B = ST_next)

if (!nrow(pairs)) {
  msg("No consecutive pairs available. Provide participants with ≥2 timepoints.")
  safe_write_csv(samples, file.path(out_dir, "sample_key.csv"))
  quit(save = "no", status = 0)
}

# ------------------------------ 8 · compute VF Jaccard & gene changes ---------
vf_mat <- vf %>% select(Participant_id, tp_lab, all_of(gene_cols))

get_vf_vec <- function(pid, tp) {
  row <- vf_mat %>% filter(Participant_id == pid, tp_lab == tp)
  if (nrow(row) != 1L) return(setNames(rep(NA_integer_, length(gene_cols)), gene_cols))
  v <- as.integer(row[1, gene_cols, drop = TRUE])
  names(v) <- gene_cols
  v
}

pairs <- pairs %>%
  mutate(
    vfA_vec = purrr::map2(Participant_id, tp_A, \(pid, tp) get_vf_vec(pid, tp)),
    vfB_vec = purrr::map2(Participant_id, tp_B, \(pid, tp) get_vf_vec(pid, tp)),
    vf_jaccard = purrr::map2_dbl(vfA_vec, vfB_vec, jaccard_distance01),
    vf_intersect = purrr::map2(vfA_vec, vfB_vec, \(a,b) names(which(a == 1L & b == 1L))),
    vf_gain     = purrr::map2(vfA_vec, vfB_vec, \(a,b) names(which(a == 0L & b == 1L))),
    vf_loss     = purrr::map2(vfA_vec, vfB_vec, \(a,b) names(which(a == 1L & b == 0L)))
  )

gene_changes_long <- pairs %>%
  select(Participant_id, tp_A, tp_B, vf_gain, vf_loss) %>%
  tidyr::pivot_longer(c(vf_gain, vf_loss), names_to = "change_type", values_to = "GENE") %>%
  tidyr::unnest(GENE) %>%
  mutate(change_type = dplyr::recode(change_type, vf_gain = "gained_in_B", vf_loss = "lost_in_B"))

safe_write_csv(gene_changes_long, file.path(out_dir, "gene_changes_long.csv"))

# ------------------------------ 9 · plasmid replicon Jaccard (optional) -------
rep_cols <- NULL
if (!is.null(rep_pid_tp)) {
  rep_cols <- setdiff(names(rep_pid_tp), c("Participant_id","tp_lab"))
  get_rep_vec <- function(pid, tp) {
    row <- rep_pid_tp %>% filter(Participant_id == pid, tp_lab == tp)
    if (nrow(row) != 1L) return(rep(NA_integer_, length(rep_cols)))
    as.integer(row[1, rep_cols, drop = TRUE])
  }
  pairs <- pairs %>%
    mutate(
      repA_vec = purrr::map2(Participant_id, tp_A, \(pid, tp) get_rep_vec(pid, tp)),
      repB_vec = purrr::map2(Participant_id, tp_B, \(pid, tp) get_rep_vec(pid, tp)),
      replicon_jaccard = if (length(rep_cols)) purrr::map2_dbl(repA_vec, repB_vec, jaccard_distance01) else NA_real_
    )
} else {
  pairs$replicon_jaccard <- NA_real_
}

# ------------------------------ 10 · merge nucmer/dnadiff metrics (optional) ---
if (exists("nucmer_tbl") && nrow(nucmer_tbl)) {
  nucmer_summ <- nucmer_tbl %>%
    group_by(Participant_id, tp_A, tp_B) %>%
    summarise(
      ANI  = suppressWarnings(max(AvgIdentity, na.rm = TRUE)),
      SNPs = suppressWarnings(min(TotalSNPs, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      ANI  = ifelse(is.infinite(ANI),  NA_real_, ANI),
      SNPs = ifelse(is.infinite(SNPs), NA_real_, SNPs)
    )
  pairs <- pairs %>%
    left_join(nucmer_summ, by = c("Participant_id","tp_A","tp_B"))
} else {
  pairs$ANI  <- NA_real_
  pairs$SNPs <- NA_real_
}

# ------------------------------ 11 · classify persistence vs replacement -------
pairs <- pairs %>%
  mutate(
    ST_same = !is.na(ST_A) & !is.na(ST_B) & ST_A == ST_B,
    vf_ok   = !is.na(vf_jaccard) & vf_jaccard <= VF_JAC_MAX,
    ani_ok  = is.na(ANI)  | (ANI  >= ANI_MIN),
    snp_ok  = is.na(SNPs) | (SNPs <= SNP_MAX),
    call    = dplyr::case_when(
      ST_same & vf_ok & ani_ok & snp_ok ~ "Persistent",
      !ST_same                           ~ "Replacement",
      ST_same & !vf_ok                   ~ "Replacement",
      ST_same & vf_ok & (!ani_ok | !snp_ok) ~ "Replacement",
      TRUE                               ~ "Uncertain"
    )
  )

# ------------------------------ 12 · participant-level summaries ---------------
persistence_by_participant <- pairs %>%
  group_by(Participant_id) %>%
  summarise(
    n_pairs = n(),
    n_persistent = sum(call == "Persistent", na.rm = TRUE),
    n_replacement = sum(call == "Replacement", na.rm = TRUE),
    n_uncertain = sum(call == "Uncertain", na.rm = TRUE),
    prop_persistent = ifelse(n_pairs > 0, n_persistent / n_pairs, NA_real_),
    .groups = "drop"
  )

# ------------------------------ 13 · write pairwise summary --------------------
pairwise_summary <- pairs %>%
  select(Participant_id, tp_A, tp_B, tpA_num, tpB_num,
         ST_A, ST_B, ST_same,
         vf_jaccard, replicon_jaccard, ANI, SNPs, call)

safe_write_csv(pairwise_summary, file.path(out_dir, "pairwise_summary.csv"))
safe_write_csv(persistence_by_participant, file.path(out_dir, "persistence_by_participant.csv"))

# ------------------------------ 14 · benchmark within vs between ---------------
# Build a between-participant baseline using VF Jaccard for pairs at the same
# timepoint. We compute explicit pairwise combinations to keep it stable.

vf_for_bench <- vf_mat %>% select(Participant_id, tp_lab, all_of(gene_cols))

# helper: compute all unordered pairwise Jaccard distances within a single tp_lab
.compute_between_for_tp <- function(df_tp) {
  # df_tp has columns: Participant_id, tp_lab, <genes...>
  if (nrow(df_tp) < 2) return(tibble::tibble())
  pids <- df_tp$Participant_id
  tpv  <- unique(as.character(df_tp$tp_lab))
  tpv  <- if (length(tpv)) tpv[1] else NA_character_
  
  mat  <- as.matrix(df_tp[, gene_cols, drop = FALSE])
  idx  <- t(combn(seq_len(nrow(mat)), 2))
  
  jd <- apply(idx, 1, function(ix) {
    jaccard_distance01(mat[ix[1], ], mat[ix[2], ])
  })
  
  tibble::tibble(
    kind       = "between",
    tp_lab     = tpv,
    pid_A      = pids[idx[, 1]],
    pid_B      = pids[idx[, 2]],
    vf_jaccard = as.numeric(jd)
  ) %>% filter(!is.na(vf_jaccard))
}

within_dist <- pairwise_summary %>%
  transmute(kind = "within", vf_jaccard = vf_jaccard) %>%
  filter(!is.na(vf_jaccard))

between_dist_all <- vf_for_bench %>%
  group_by(tp_lab) %>% group_split() %>%
  purrr::map_dfr(.compute_between_for_tp)

# size-match the between distribution to the number of within pairs (if needed)
set.seed(1)
n_within <- nrow(within_dist)
between_dist <- if (nrow(between_dist_all) > n_within) {
  dplyr::slice_sample(between_dist_all, n = n_within) %>% transmute(kind, vf_jaccard)
} else {
  between_dist_all %>% transmute(kind, vf_jaccard)
}

benchmark_within_vs_between <- dplyr::bind_rows(within_dist, between_dist)
safe_write_csv(benchmark_within_vs_between,
               file.path(out_dir, "benchmark_within_vs_between.csv"))

# ------------------------------ 15 · plots: VF Jaccard distributions -----------
p_vf_jacc <- ggplot(benchmark_within_vs_between,
                    aes(x = vf_jaccard, fill = kind)) +
  geom_density(alpha = 0.45, adjust = 1.2) +
  theme_bw(base_size = 12) +
  labs(title = "VF Jaccard distance: within-person vs between-person",
       x = "Jaccard distance (0 = identical, 1 = disjoint)",
       y = "Density", fill = "") +
  guides(fill = guide_legend(override.aes = list(alpha = 0.6)))

# rename existing density save (change your earlier save line)
safe_save_plot(file.path(plot_dir, "vf_jaccard_within_vs_between_density.png"),
               p_vf_jacc, width = 7, height = 4.5, dpi = 300)

## Boxplot/jitter for within vs between with Wilcoxon p-value (uses 'kind' & 'vf_jaccard')
bench_df <- readr::read_csv(file.path(out_dir, "benchmark_within_vs_between.csv"), show_col_types = FALSE)

wilcox_p <- NA_real_
if (sum(bench_df$kind == "within", na.rm = TRUE) >= 5 &&
    sum(bench_df$kind == "between", na.rm = TRUE) >= 5) {
  wilcox_p <- suppressWarnings(
    wilcox.test(bench_df$vf_jaccard[bench_df$kind == "within"],
                bench_df$vf_jaccard[bench_df$kind == "between"],
                alternative = "less")$p.value
  )
}
subtitle_txt <- if (!is.na(wilcox_p)) paste0("Wilcoxon (Within < Between) p = ", signif(wilcox_p, 3)) else "Wilcoxon p = NA"

pD <- bench_df %>%
  ggplot(aes(kind, vf_jaccard, fill = kind)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3) +
  labs(title = "Benchmark VF distances", subtitle = subtitle_txt, x = NULL, y = "VF Jaccard", fill = "")
safe_save_plot(file.path(plot_dir, "benchmark_within_vs_between.png"), pD, 6, 4.5)
# (optional) drop a text file like the initial script did
writeLines(paste0("wilcox_p (Within < Between) = ",
                  ifelse(is.na(wilcox_p), "NA", signif(wilcox_p, 6))),
           file.path(out_dir, "benchmark_wilcoxon.txt"))

# ------------------------------ 16 · plots: call summary -----------------------
p_call <- pairs %>%
  count(call, name = "n") %>%
  mutate(call = factor(call, levels = c("Persistent","Replacement","Uncertain"))) %>%
  ggplot(aes(call, n)) + geom_col() + theme_bw(base_size = 12) +
  labs(title = "Pairwise calls", x = NULL, y = "Count")

safe_save_plot(file.path(plot_dir, "pairwise_calls_bar.png"),
               p_call, width = 5, height = 4, dpi = 300)

# ------------------------------ 17 · optional: scatter VF vs ANI ---------------
if (any(!is.na(pairs$ANI))) {
  p_scatter <- pairs %>%
    ggplot(aes(x = vf_jaccard, y = ANI, shape = call)) +
    geom_point(size = 2) +
    theme_bw(base_size = 12) +
    labs(title = "Within-person pairs: VF Jaccard vs ANI",
         x = "VF Jaccard distance",
         y = "ANI (%)",
         shape = "Call")
  safe_save_plot(file.path(plot_dir, "vf_vs_ani_scatter.png"),
                 p_scatter, width = 6, height = 4.5, dpi = 300)
}

# A) Distribution of VF Jaccard by classification (matches initial script idea)
pA <- pairs %>%
  filter(!is.na(vf_jaccard)) %>%
  ggplot(aes(vf_jaccard, fill = call)) +
  geom_histogram(bins = 30, alpha = 0.8, position = "identity") +
  geom_vline(xintercept = VF_JAC_MAX, linetype = "dashed") +
  labs(title = "Within-person VF Jaccard distances",
       subtitle = paste0("Dashed line = threshold ", VF_JAC_MAX),
       x = "VF Jaccard distance (0 = identical, 1 = disjoint)",
       y = "Pairs", fill = "Call")
safe_save_plot(file.path(plot_dir, "vf_jaccard_by_classification.png"), pA, 8, 5)

# B) Trajectories of VF distance (consecutive TPs) — requires tp_num per TP
# Build a reliable PID×TP -> tp_num lookup (avoid relying on 'samples')
tp_lookup <- vf %>%
  dplyr::select(Participant_id, tp_lab, tp_num) %>%
  dplyr::distinct()

# B) Trajectories of VF distance (consecutive TPs)
traj_df <- pairs %>%
  dplyr::left_join(tp_lookup, by = c("Participant_id","tp_A" = "tp_lab")) %>%
  dplyr::rename(tp_num_A = tp_num) %>%
  dplyr::left_join(tp_lookup, by = c("Participant_id","tp_B" = "tp_lab")) %>%
  dplyr::rename(tp_num_B = tp_num) %>%
  dplyr::filter(!is.na(tp_num_A), !is.na(tp_num_B))

pC <- traj_df %>%
  ggplot2::ggplot(ggplot2::aes(tp_num_A, vf_jaccard, group = Participant_id, color = call)) +
  ggplot2::geom_line(alpha = 0.4) +
  ggplot2::geom_point() +
  ggplot2::scale_x_continuous(breaks = sort(unique(traj_df$tp_num_A))) +
  ggplot2::labs(title = "Trajectories of within-person VF distance (consecutive TPs)",
                x = "Numeric timepoint (A)", y = "VF Jaccard")
safe_save_plot(file.path(plot_dir, "vf_trajectory_consecutive.png"), pC, 8, 5)


# ------------------------------ 18 · write extra debug aids --------------------
# (These are handy when diagnosing mismatches or empty gene change tables.)
safe_write_csv(samples, file.path(out_dir, "sample_key.csv"))
safe_write_csv(pairs,    file.path(out_dir, "pairs_with_metrics.csv"))

msg("All done. Saved benchmark table and plots in results/within_person/ and results/plots/within_person/")

# Tiny README
readme_lines <- c(
  "Within-person strain comparison (Yellow RoUTIne):",
  paste0("- Persistence rule: ST_same & VF_Jaccard ≤ ", VF_JAC_MAX,
         if (exists("nucmer_tbl") && nrow(nucmer_tbl))
           paste0(" & ANI ≥ ", ANI_MIN, "% (or SNPs ≤ ", SNP_MAX, ")")
         else ""),
  "- Missing core metrics or discordant thresholds → 'Uncertain'.",
  "- See pairwise_summary.csv for per-pair metrics; debug/ for conflicts/missing data.",
  "- Optional: replicon Jaccard is computed only if Isolate_ID→Participant×Timepoint mapping is available."
)
writeLines(readme_lines, file.path(out_dir, "README_11_within_person.txt"))

# Verification table
verify_file <- function(p) file.exists(p) && isTRUE(file.size(p) > 0)
must_have <- c(file.path(out_dir, "pairwise_summary.csv"),
               file.path(out_dir, "persistence_by_participant.csv"),
               file.path(out_dir, "gene_changes_long.csv"))
ok_vec <- vapply(must_have, verify_file, logical(1))
verif_tbl <- tibble::tibble(file = must_have, exists_nonempty = ok_vec)
safe_write_csv(verif_tbl, file.path(out_dir, "verification_table.csv"))

# Session info snapshot
sess_txt <- capture.output(sessionInfo())
writeLines(sess_txt, file.path(out_dir, "sessionInfo.txt"))