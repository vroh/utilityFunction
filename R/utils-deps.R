assert_suggested <- function(pkg, fn = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- if (is.null(fn)) {
      sprintf("Package '%s' is required.", pkg)
    } else {
      sprintf(
        "Package '%s' is required for `%s()`. Install it with install.packages('%s').",
        pkg, fn, pkg
      )
    }
    cli::cli_abort(msg)
  }
}
