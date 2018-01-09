
# Extract scale from an image with some graph paper on
# (either) left/right side
# Tuning points are detailed with `***` below
# See https://github.com/vbonhomme/Momocs_recipes/blob/master/EXTRACT_PaperScale_example.jpg

# Dependencies -----
library(jpeg)

# Given a jpg image, find on the third channel,
# and on several column in both left/right edges
# if graph paper gridlines can be found. Profiled columns
# actually spanning graph paper are then slightly treated
# and distance between (blue) peaks are calculated. Then
# return a median (of positive columns) of
# medians (of distance between peaks for those columns)
graph_paper <- function(path){
  peak_distances <- function(x){
    # given a vector (eg a column profile of image intensity)
    # return the distances between peaks
    # if real peaks, are found in a given column,
    # ie if the graph paper spans this columns,
    # then the CV (coefficient of variation) is expected to be very low
    # below, we thus try different columns on the left/right sides
    # then select those with a very low CV

    # remove 20% on each side to remove included scale bar
    # and other artifacts
    # *** proportion of a single vector to exclude on both ends
    exclude_edges_proportion = 0.1
    r <- round(length(x)*exclude_edges_proportion)
    x <- x[r:((length(x)-r))]
    # number of smoothing iterations ~~~~~
    n <- round(length(x)/50)
    # moving average and remove resulting NAs
    x <- stats::filter(x, rep(1/n,n), sides=2)
    x <- na.omit(x)

    # retain only the upper decile ~~~~~
    # to avoid valleys noise
    # *** how high to cut on each peak
    x[x<quantile(x, probs=0.9)] <- NA

    # peak detection ~~~~~
    #inflexion points
    pos_inflexion_points <- diff(sign(diff(x))) == -2
    # distances (in pixels) between peaks
    diff(which(pos_inflexion_points))
  }


  # read jpg ~~~~~
  x <- readJPEG(path)

  # retain only the third channel ~~~~~
  # and take the complementary (just to have pos. peaks)
  # *** which channel to choose
  x <- (1- x[,, 3])

  # find columns where to profile peaks ~~~~~
  # 5 in the left 1/5th ; same in the right 1/5th
  # *** number of columns to profile, per side
  nb_per_side = 5
  # *** proportion on each side, where to profile columns
  side_proportion = 0.20
  cols_profiled <- c(#left side profiles
    round(seq(1,
            round(ncol(x)*side_proportion),
            length=nb_per_side)),
    #right side profiles
    round(seq(round(ncol(x)*(1-side_proportion)),
              ncol(x),
              length=nb_per_side)))

  # calculate peak for chosen cols ~~~~~
  ds <- apply(x[, cols_profiled], 2, peak_distances)
  CVs <- 100*sapply(ds, function(.) sd(.)/mean(.))
  #return(CVs)
  # final filterting and peak distance calculation ~~~~~
  # return for CVs below 2%
  # the median of the medians
  # *** CV tolerance
  median(sapply(ds[which(CVs < 2)], median))
}

# Example of use ------------
# change this with your folder
lf <- list.files("data/1431_cGatetN/1431_cGatetN_Initiale/",
                 full.names = TRUE, pattern="jpg$")

# single file (very bad idea)
img <- download.file("https://raw.githubusercontent.com/vbonhomme/Momocs_recipes/master/EXTRACT_PaperScale_example.jpg",
                      "img.jpg")
graph_paper("img.jpg")
file.remove("img.jpg")

# on a list of files
scale <- sapply(lf, graph_paper)
# then
data.frame(scale=scale, file=lf)

# Benchmark --------
# microbenchmark::microbenchmark(
#     lapply(lf, graph_paper),
#   times=5)
## 10 images treated per sec.