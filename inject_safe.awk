BEGIN{ added=0 }
{
  print
  if (!added && $0 ~ /^options\([^)]*stringsAsFactors/) {
    print "## [patched] define safe_future_lapply early"
    print "safe_future_lapply <- function(X, FUN, ...) lapply(X, FUN, ...)"
    added=1
  }
}
