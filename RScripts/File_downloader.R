# Specify URL where file is stored
url <- "ftp://ftp.ensembl.org/pub/release-96/gtf/mus_musculus/Mus_musculus.GRCm38.96.gtf.gz"
# Specify destination where file should be saved
destfile <- "G:/Mi unidad/Doctorado/Coding/Data/Mus_musculus.GRCm38.96.gtf.gz"
# Apply download.file function in R
download.file(url, destfile)
