# pam 0.2.0

* DESCRIPTION overhauled: `Authors@R` format, CRAN-quality title and
  description, `BugReports` URL, `Suggests` for testthat/knitr/rmarkdown.
* Replaced `httr`/`stringr`/`magrittr` dependencies with base R equivalents;
  added `checkmate` for input validation.
* Switched from magrittr pipe (`%>%`) to base R pipe (`|>`).
  Minimum R version raised to 4.1.0.
* Added `checkmate::assert_*` guards in `get_geo_pam()` and `url_exists()`.
* Roxygen documentation rewritten for all exported and internal functions.
* Added testthat test suite (`tests/testthat/`).
* Added introductory vignette (`vignettes/pam-intro.Rmd`).
* Added `NEWS.md`.
* Added R CMD check CI workflow (ubuntu, macOS, windows).
* Fixed `get_rodentia_mexico()` docstring (was incorrectly labelled
  "afrosoricida").

# pam 0.1.0

* Initial release with `get_geo_pam()`, `get_rodentia_mexico()`, and
  `get_world_pca()`.
