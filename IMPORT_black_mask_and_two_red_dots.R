import_black_and_red1 <-
function(path){
  # reads raw jpg
  img <- jpeg::readJPEG(path)
  # grabs the red channel
  
  R <- img[,,1]
  # grabs the green channel and prepare masks
  grey_mask <- black_mask <- G <- img[,,2]
  
  # all pixels darker than 0.6
  # for the green channel form the grain
  black_mask[black_mask > 0.6]  <- 1
  black_mask[black_mask <= 0.6] <- 0
  
  # pixels between 0.4 and 0.6
  # for the green channel form the embryo
  grey_mask[G <  0.4 | G  > 0.4]  <- 1
  grey_mask[G >= 0.4 & G <= 0.6]  <- 0
  
  # in this case, Conte will need a center which is not
  # necessarily the center of the image
  # so we calculate the barycenter of all grey and black pixels
  grey_center  <- which(grey_mask==0,  arr.ind=TRUE) %>%
    apply(2, mean) %>% round()
  black_center <- which(black_mask==0, arr.ind=TRUE) %>%
    apply(2, mean) %>% round()
  
  # ensure that we start on a black pixel
  while (black_mask[black_center[1], black_center[2]] != 0){
    black_pixels <- which(black_mask==0, arr.ind=TRUE)
    black_center <- black_pixels[sample(1:nrow(black_pixels), size = 1),]
  }
  # now we import outlines from each mask and center
  black_out <- import_Conte(black_mask, black_center)
  
  # grey_out <- import_Conte(grey_mask, grey_center)
  # pixels above 0.9 on red, but below 0.1 on green
  # are the two landmarks
  ldk <- which(R > 0.9 & G < 0.1, arr.ind=TRUE)
  # some checks
  # if (nrow(ldk)<2)
  #   stop(path, " has less than 2 red pixels")
  # if (nrow(ldk)>2)
  #   stop(path, " has more than 2 red pixels")
  
  # because the way human read jpg and the matricial
  # way to code it are different
  ldk <- ldk[, 2:1]
  ldk[, 2] <- nrow(img) - ldk[, 2]
  
  # we also actually need their position for later Out
  # that is the ids of the closest points on the grain outline
  ldk_pos <- edm_nearest(ldk, black_out, TRUE)$pos
  
  # we return this beauty
  list(out=black_out, ldk=ldk, ldk_pos=ldk_pos)
}
