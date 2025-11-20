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
      # x[line_idx] <- paste0(number, ". ", x[line_idx])
      x[line_idx] <- paste0("1", ". ", x[line_idx])
    }
  }
}

x <- c(readLines("preamble.md"), "\n", x)

writeLines(x, "README.md")
rmarkdown::render("citations.md", output_file = "citations-preview.html")

# Count entries in each category
categories <- c("journals", "theses", "preprints", "reports")
category_files <- paste0(categories, ".md")

cat("\n=== Citation Statistics ===\n\n")

total_count <- 0
total_since_2021 <- 0

# Store counts for file output
counts <- list()

for (i in seq_along(category_files)) {
  if (file.exists(category_files[i])) {
    lines <- readLines(category_files[i])

    # Count citation lines (non-empty, not headers, not whitespace)
    citation_lines <- which(nchar(lines) > 0 & !grepl("^#", lines) & !grepl("^\\s*$", lines))
    count <- length(citation_lines)

    # Extract years and count since 2021 (one year per citation)
    years <- sapply(lines[citation_lines], function(line) {
      year_matches <- regmatches(line, gregexpr("\\b(19|20)\\d{2}\\b", line))[[1]]
      if (length(year_matches) > 0) {
        # Take the first year found (usually the publication year)
        return(as.numeric(year_matches[1]))
      } else {
        return(NA)
      }
    })
    since_2021 <- sum(years >= 2021, na.rm = TRUE)

    cat(sprintf("%-12s: %3d total, %3d since 2021\n",
                tools::toTitleCase(categories[i]), count, since_2021))

    counts[[categories[i]]] <- list(total = count, since_2021 = since_2021)

    total_count <- total_count + count
    total_since_2021 <- total_since_2021 + since_2021
  }
}

cat(sprintf("\n%-12s: %3d total, %3d since 2021\n",
            "TOTAL", total_count, total_since_2021))
cat("\n")

# Write to file in specified format
papers <- counts$journals$total
papers_2021 <- counts$journals$since_2021
preprints <- counts$preprints$total
preprints_2021 <- counts$preprints$since_2021
papers_preprints <- papers + preprints
papers_preprints_2021 <- papers_2021 + preprints_2021
reports <- counts$reports$total
reports_2021 <- counts$reports$since_2021
theses <- counts$theses$total
theses_2021 <- counts$theses$since_2021

output_lines <- c(
  sprintf("sdmTMB_cites_papers, %d", papers),
  sprintf("sdmTMB_cites_preprints, %d", preprints),
  sprintf("sdmTMB_cites_papers_preprints, %d", papers_preprints),
  sprintf("sdmTMB_cites_reports, %d", reports),
  sprintf("sdmTMB_cites_theses, %d", theses),
  sprintf("sdmTMB_cites_papers_2021, %d", papers_2021),
  sprintf("sdmTMB_cites_preprints_2021, %d", preprints_2021),
  sprintf("sdmTMB_cites_papers_preprints_2021, %d", papers_preprints_2021),
  sprintf("sdmTMB_cites_reports_2021, %d", reports_2021),
  sprintf("sdmTMB_cites_theses_2021, %d", theses_2021)
)

writeLines(output_lines, "citation_counts.txt")

setwd(here::here("."))
