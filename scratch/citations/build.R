# # Build sdmTMB citations list
#
# These are papers/theses/preprints/reports that use sdmTMB
# as a modelling platform.
#
# To build the list:
#
# 1. Add a reference to the Zotero collection.
# 2. That collection will auto-export to the .bib files shown here.
# 3. Source build.R in R.
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

# Add numbering to citations within each section (counting from top)
# First pass: identify sections and mark citation lines
section_starts <- which(grepl("^# ", x))
citation_lines <- which(nchar(x) > 0 & !grepl("^#", x) & !grepl("^\\s*$", x))

# For each section, number citations from top to bottom (1, 2, 3...)
for (i in seq_along(section_starts)) {
  section_start <- section_starts[i]
  section_end <- if (i < length(section_starts)) section_starts[i + 1] - 1 else length(x)

  # Find citation lines in this section
  section_citations <- citation_lines[citation_lines > section_start & citation_lines <= section_end]

  if (length(section_citations) > 0) {
    # Number from top to bottom (1, 2, 3...)
    for (j in seq_along(section_citations)) {
      line_idx <- section_citations[j]
      # Assign number counting from top
      number <- j
      x[line_idx] <- paste0(number, ". ", x[line_idx])
    }
  }
}

x <- c(readLines("preamble.md"), "\n", x)

writeLines(x, "README.md")
rmarkdown::render("citations.md", output_file = "citations-preview.html")

setwd(here::here("."))
