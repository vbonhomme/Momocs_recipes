# VB
# 8/1/18
# Extract scale from an image with some graph paper on
# (either) left/right side

# Dependencies -----
library(jpeg)

# Given a jpg image, find on the right channel,
# if graph paper gridlines can be found. Profiled columns
# actually spanning graph paper are then slightly treated
# and distance between (blue) peaks are calculated. Then
# return a median (of positive columns) of
# medians (of distance between peaks for those columns)
# Use debug mode to tune tuning points
paper_scale <- function(path,
                        quantile_white=0.95,
                        quantile_peaks=0.75,
                        smoothing_intensity=1/50,
                        channel=3,
                        debug=FALSE){

  # import jpg and pick third channel
  img <- jpeg::readJPEG(path)[,,channel]

  # normalize
  img <- img-min(img)
  img <- img/max(img)

  # select only whitish columns (where graph paper is supposed to be)
  mc <- apply(img, 2, mean)
  x <- img[, which(mc > quantile(mc, probs=quantile_white))]


  if (debug) {
    cat("\n", path)
    layout(matrix(1:4, nc=2))
    Momocs::img_plot0(img)
    Momocs::img_plot0(x)
  }


  # from these columns, profile intensity on rows
  x <- apply(x, 1, mean)

  if (debug) plot(x, type="l")

  # smoothes with moving average
  # and return a numeric with upwards peaks
  n <- round(length(x) * smoothing_intensity)
  xs <- stats::filter(x, rep(1/n,n), sides=2)
  xs <- 1 - as.numeric(xs[!is.na(xs)])
  # thresholdize peaks
  xs <- as.numeric(xs>=quantile(xs, probs=quantile_peaks))

  if (debug) plot(xs, type="l")

  # diff between their starting points
  d <- diff(which(diff(xs)==1))

  # if found, pcik the median, NA otherwise
  res <- ifelse(length(d)>0, median(d), NA)

  if (debug){
    layout(matrix(1))
    cat("\n")
    print(d)
    cat("\n")
    print(res)
    cat(rep("*", 20))
    Sys.sleep(2)
  }
  return(res)
}

#### On a single image:
img <- download.file("https://raw.githubusercontent.com/vbonhomme/Momocs_recipes/master/EXTRACT_PaperScale_example.jpg",
                     "img.jpg")
paper_scale("img.jpg")
paper_scale("img.jpg", debug=TRUE)
file.remove("img.jpg")


#### How to use it on a batch:
lf <- list.files("data/0129-wMouBo5_Initiale/", # change this
                 full=T,
                 ignore.case=TRUE, pattern="jpg$") %>% sample

lapply(lf, paper_scale, debug=T) # debug
sapply(lf, paper_scale)          # prod

