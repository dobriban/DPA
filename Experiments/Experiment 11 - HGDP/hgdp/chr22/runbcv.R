#Select number of factors by BCV on HGDP data
hgdp <- read.csv("C:/Dropbox/Projects/Deflation/Experiments/Experiment 11 - HGDP/hgdp/chr22/hgdp.txt", header=FALSE)
library(esaBcv)
start.time <- Sys.time()
r <-EsaBcv(hgdp, r.limit = 20, niter = 3, nRepeat = 12, only.r = T) 
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

#Error in svd(m) : infinite or missing values in 'x'
#[1] "Encounter problematic data!"
#Error: $ operator is invalid for atomic vectors