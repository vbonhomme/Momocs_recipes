# to import from optimas .asc files, with covariables
library(Momocs)

# import a single .asc when given its path
import_asc1 <- function(path){
  x <- readLines(path)
  # removes all space characters beginning a line
  # and all empty lines
  x <- gsub("^[[:space:]]+", "", x[nchar(x)>0])
  # grab the first line, it will be the colnames
  fac_cn <- unlist(strsplit(x[1], "[[:space:]]"))
  # though we also need to remove "Sampled points"
  fac_cn <- fac_cn[-length(fac_cn)]
  # turns all other lines into numeric values
  xs <- strsplit(x[-1], "[[:space:]]") %>% sapply(as.numeric)
  # deduces xs lines also have fac values
  fac_ids <- which(sapply(xs, length) > 2)
  # adds the last line+1 id for the loop to come
  fac_ids <- c(fac_ids, length(xs)+1)
  #deduces the number of covariables on the first numeric line
  n <- length(xs[[1]])-2
  # on these lines extract all but the last two
  # and turn them into a data.frame
  fac <- do.call("rbind", xs[fac_ids])[, 1:n]
  colnames(fac) <- fac_cn
  # only retain the last two numeric values
  xs %<>% lapply(function(.) .[(length(.)-1) : (length(.))])
  # prepares a list to store coordinates
  coo <- vector("list", length(fac_ids)-1)
  # loop over coordinates
  # cut all coordinates lines
  # and turn them into matrices
  for (i in seq_along(coo)){
    coo[[i]] <- do.call("rbind", xs[fac_ids[i] : (fac_ids[i+1]-1)])
  }
  # return a list with the coo and the fac
  list(coo=coo, fac=fac)
}

# import many asc provided as paths to individual .asc
# eg as a vector returned by list.files
import_asc <- function(paths){
  res <- lapply(paths, import_asc1)
  coo <- lapply(res, function(.) .$coo) %>% do.call("rbind", .)
  fac <- lapply(res, function(.) .$fac) %>% do.call("rbind", .)
  list(coo=coo, fac=fac)
}

# Practically ####
# single import ----
res <- import_asc1("~/Desktop/data/um1_bavella.asc")
Out(res$coo, res$fac) %>% panel

# many .asc files ----
# let's say all .asc are in data/
lf <- list.files("~/Desktop/data/", full.names = TRUE)
res <- import_asc(lf)
Out(res$coo, res$fac)
