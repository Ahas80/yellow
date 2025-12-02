BEGIN{ inside=0; depth=0;
  repl = "progress_tick <- function(stage, pid = NA_character_, status = \"tick\") {\n  invisible(NULL)\n}\n"
}
{
  # Start of progress_tick <- function(...
  if (!inside && $0 ~ /^[[:space:]]*progress_tick[[:space:]]*<-[[:space:]]*function[[:space:]]*\(/) {
    print "# [patched] progress_tick (no-op)"
    printf "%s", repl
    inside=1; depth=0
    nopen  = gsub(/\{/,"&"); nclose = gsub(/\}/,"&"); depth += nopen - nclose
    next
  }
  # Consume original body until braces balance
  if (inside) {
    nopen  = gsub(/\{/,"&"); nclose = gsub(/\}/,"&"); depth += nopen - nclose
    if (depth <= 0) inside=0
    next
  }
  # Lines outside the function pass through unchanged
  print
}
