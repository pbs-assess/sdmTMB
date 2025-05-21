# # Build sdmTMB citations list
#
# These are papers/theses/preprints/reports that use sdmTM
# as a modelling platform.
#
# To build the list:
#
# 1. Add a reference to the Zotero collection.
# 2. That collection will auto-export to the .bib files shown here.
# 3. Source `build.R` in R.
# 4. That will update the README.md file.

setwd(here::here("scratch", "citations"))

rmarkdown::render("journals.Rmd")
rmarkdown::render("theses.Rmd")
rmarkdown::render("preprints.Rmd")
rmarkdown::render("reports.Rmd")

x <- readLines("theses.md")

x <- gsub("\\\\\\[", "", x)
x <- gsub("\\\\\\]", "", x)
x <- gsub("\\{\\{", "", x)
x <- gsub("\\}\\}", "", x)
x <- gsub("\\{\\{", "", x)
x <- gsub("\\}\\}", "", x)

x <- gsub("PhD thesis\\,", "PhD Thesis\\.", x)
x <- gsub("\\* PhD Thesis\\.", "\\*. PhD Thesis\\.", x)
x <- gsub("MSc Thesis\\,", "MSc Thesis\\.", x)
x <- gsub("\\* MSc Thesis\\.", "\\*. MSc Thesis\\.", x)
x <- gsub("\\*Milvus\\* \\*Milvus\\*", "Milvus milvus", x)
writeLines(x, "theses.md")

r <- c(
  "preprints.md",
  "journals.md",
  "theses.md",
  "reports.md"
)

x <- list()
for (i in seq_along(r)) {
  x[[i]] <- readLines(r[i])
  x[[i]] <- c(x[[i]], "\n\n")
}

x <- do.call("c", x)

x <- gsub("\\\\\\[", "", x)
x <- gsub("\\\\\\]", "", x)
x <- gsub("\\{\\{", "", x)
x <- gsub("\\}\\}", "", x)
x <- gsub("\\{\\{", "", x)
x <- gsub("\\}\\}", "", x)

x <- c(readLines("preamble.md"), "\n", x)

writeLines(x, "README.md")
rmarkdown::render("citations.md", output_file = "citations-preview.html")

setwd(here::here("."))
