## Prepare plots and tables for report

## Before:
## After:

library(icesTAF)
taf.library(spict)
taf.library(TMB)
library(dplyr)

library(rmarkdown)

mkdir("report")
cp("bootstrap/initial/report/*", "report/")

outdir <- "report/"

output_format <- NULL # "all"
quiet <- FALSE

icesTAF::msg("Report: Making catch working document")
render("report/tur.27.3a_catch_WD.Rmd", output_dir = outdir,
       #output_format = output_format,
       clean = TRUE, quiet = quiet,  encoding = 'UTF-8')

icesTAF::msg("Report: Making assessment working document")
render("report/tur.27.3a_assessment_WD.Rmd", output_dir = outdir,
       #output_format = output_format,
       clean = TRUE, quiet = quiet,  encoding = 'UTF-8')

icesTAF::msg("Report: Making catch and assessment presentation")
render("report/tur.27.3a_assessment_Presentation.Rmd",
       output_dir = outdir,
       #output_format = output_format,
       clean = TRUE, quiet = quiet,  encoding = 'UTF-8')

