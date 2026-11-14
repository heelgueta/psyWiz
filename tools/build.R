# Build the jamovi module (.jmo) and install it into the local jamovi.
#
#   Rscript tools/build.R [jamovi home] [extra R library with lavaan]
#
# Run from the repository root with jamovi's own R and with jmvtools on the
# library path; docs/DEVELOPING.md has the exact command for each platform.
#
# Normally `jmvtools::install()` is all you need. This script exists because
# psyWiz imports lavaan, which jamovi ships inside its `jmv` module, and on
# some installations (the Linux flatpak) jmvtools does not look there and
# stops with "dependency 'lavaan' is not available". The compiler has an
# `--rlibs` flag for this, which jmvtools::install() does not expose.
args <- commandArgs(trailingOnly = TRUE)
home <- if (length(args) >= 1L) args[1L] else "/app"
rlibs <- if (length(args) >= 2L) args[2L] else "/app/lib/jamovi/modules/jmv/R"

jmc <- system.file("node_modules", "jamovi-compiler", "index.js", package = "jmvtools")
if (!nzchar(jmc)) stop("jmvtools is not installed in this R library (see docs/DEVELOPING.md)")
status <- system2(node::node(), c(shQuote(jmc), "--install", ".", "--home", shQuote(home),
                                  "--rlibs", shQuote(rlibs)))
quit(status = status)
