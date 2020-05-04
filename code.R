#identify the full names of the files of interest
files <- list.files( ".", ".coverage", full.names = T)

#define a function for computing the coverage
computeCoverage <- function(file){
  cov <- read.table(file, sep = "\t", header = F, dec = ".")
  mean.cov <- sapply(split(cov, cov$V1), function(d) sum(d$V2 * d$V5))
  mean.cov
}

#apply the function to the coverage files
covs <- sapply(c("./B01.coverage", "./B02.coverage"), computeCoverage)
colnames(covs) <- gsub("^.*/(B[0-9]{2})\\.coverage$", "\\1", colnames(covs))

#compute per contig coverage 
index <- grep("genome", rownames(covs))
covs <- covs[-index,]

#plot results with the built-in heatmap function
pal <- colorRampPalette(c("black","cyan", "white", "yellow", "red"))
heatmap(covs, Rowv = NA, Colv = NA, scale = "none", col = pal(256))

#remove outliers and plot the results again
covs.mod <- covs[-grep("CP004061.1", rownames(covs)),]
covs.mod <- covs.mod[-grep("JH370354.1", rownames(covs.mod)),]
heatmap(covs.mod, Rowv = NA, Colv = NA, scale = "none", col = pal(256))

#R code for gene coverage analysis
files <- list.files("./", "genecov.bed", full.names = T)

computeGeneCoverage <- function(file){
  sample <- gsub("^.*/(B[0-9]{2})\\_genecov.bed$", "\\1", file)
  
  genecov <- read.table(file, sep = "\t", header = F, dec = ".",
                        stringsAsFactors = F)
  
  rgx <- "^.*;protein_id=([^;]+)(;.*)?$"
  data.frame(sample = sample,
             id = gsub(rgx, "\\1", genecov$V9),
             cov = genecov$V10,
             len = genecov$V12,
             stringsAsFactors = F)
}

genecov <- lapply(files, computeGeneCoverage)


# RPKM and TPM normalization
rpkm <- function(cov, len, n){
  scale <- n/1e6
  len.kb <- len/1000
  ( cov/scale ) / len.kb
}

tpm <- function(cov, len){
  len.kb <- len/1000
  rpk <- cov / len.kb
  scale <- sum(rpk)/1e6
  rpk / scale
}

#R code for rpkm and tmp calculation 
genecov <- mapply(function(x, n){
  x$rpkm <- rpkm(cov = x$cov, len = x$len, n = n)
  x$tpm <- tpm(cov = x$cov, len = x$len)
  x
}, genecov, n_sequences, SIMPLIFY = F)

mat <- do.call(cbind, lapply(genecov, function(x){
  data.frame(rpkm = x$rpkm, tpm = x$tpm, count = x$cov)
}))
colnames(mat) <- paste(rep(c("B01", "B02"), each = ncol(mat)/2), 
                       colnames(mat))


#R code of plotting gene coverage histograms
old <- par(mfrow=c(2,3))
for(i in 1:ncol(mat)){
  mu <- mean(mat[[i]])
  hist(mat[[i]], breaks = 100, 
       main = colnames(mat)[i],
       xlab = "", 
       col = "deepskyblue2")
  abline(v = mu, col = "red")
}
par(old)


#R code of plotting gene coverage heatmaps
heatmap(as.matrix(mat), Rowv = NA,
        scale = "none", col = pal(256),
        margins = c(8, 0),
        labRow = NA)

