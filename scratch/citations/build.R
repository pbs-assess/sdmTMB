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

writeLines(x, "citations.md")
rmarkdown::render("citations.md", output_file = "citations-preview.html")

setwd(here::here("."))
