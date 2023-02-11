
pre_geom_morphospace <- function(PCA, xax, yax, pos.shp, nb.shp = 24, nr.shp = 6, nc.shp = 5,
                                 amp.shp = 1, rotate.shp = 0, flipx.shp = FALSE, flipy.shp = FALSE,
                                 size.shp = 1, wdw = max(.wdw()), pts.shp = 60, col.shp = "#00000011",
                                 border.shp = "#00000055", lwd.shp = 1) {
  # stop if more than 4 partitions
  if (length(PCA$method) > 4 | is.null(PCA$method)) {
    stop("morphospacePCA can only handle up to 4 partitions of coefficients")
  }
  # grabs important components from PCA objects
  xy     <- PCA$x[, c(xax, yax)]
  rot    <- PCA$rotation[, c(xax, yax)]
  mshape <- PCA$mshape
  pos    <- morphospace_positions(xy,
                                  pos.shp = pos.shp,
                                  nb.shp = nb.shp,
                                  nr.shp = nr.shp, nc.shp = nc.shp)
  method <- PCA$method
  lm <- length(method)
  # recycle if required
  if (length(rotate.shp) != lm) {
    rotate.shp <- rep(rotate.shp, lm)
  }
  if (length(flipx.shp) != lm) {
    flipx.shp <- rep(flipx.shp, lm)
  }
  if (length(flipy.shp) != lm) {
    flipy.shp <- rep(flipy.shp, lm)
  }
  if (length(size.shp) != lm) {
    size.shp <- rep(size.shp[1], lm)
  }

  size.shp.final <- (size.shp * wdw / 14) / ifelse(lm < 2, 1, 2)
  d <- mean(size.shp.final) / 2
  if (lm == 1) {
    dx <- 0
    dy <- 0
  }
  if (lm == 2) {
    dx <- c(0, 0)
    dy <- c(d, -d)
  }
  if (lm == 3) {
    dx <- c(0, -d, d)
    dy <- c(d, -d, -d)
  }
  if (lm == 4) {
    dx <- c(-d, d, -d, d)
    dy <- c(d, d, -d, -d)
  }
  if (lm == 1) {
    col.start <- 1
    col.end <- length(mshape)
  }
  else {
    col.start <- cumsum(PCA$cuts) - PCA$cuts + 1
    col.end <- cumsum(PCA$cuts)
  }

  for (i in seq(along = method)) {
    shp <- NULL
    plot.method <- NULL
    ids <- col.start[i]:col.end[i]
    if (method[i] == "efourier") {
      shp <- Momocs:::PCA2shp_efourier(pos = pos, rot = rot[ids, ], mshape = mshape[ids], amp.shp = amp.shp, pts.shp = pts.shp)
      shp <- lapply(shp, coo_close)
      plot.method <- "poly"
    }
    if (method[i] == "rfourier") {
      shp <- PCA2shp_rfourier(pos = pos, rot = rot[ids, ], mshape = mshape[ids], amp.shp = amp.shp, pts.shp = pts.shp)
      shp <- lapply(shp, coo_close)
      plot.method <- "poly"
    }
    if (method[i] == "sfourier") {
      shp <- PCA2shp_sfourier(pos = pos, rot = rot[ids, ], mshape = mshape[ids], amp.shp = amp.shp, pts.shp = pts.shp)
      shp <- lapply(shp, coo_close)
      plot.method <- "poly"
    }
    if (method[i] == "tfourier") {
      shp <- PCA2shp_tfourier(pos = pos, rot = rot[ids, ], mshape = mshape[ids], amp.shp = amp.shp, pts.shp = pts.shp)
      shp <- lapply(shp, coo_close)
      plot.method <- "poly"
    }
    if (method[i] == "dfourier") {
      shp <- PCA2shp_dfourier(pos = pos, rot = rot[ids, ], mshape = mshape[ids], amp.shp = amp.shp, pts.shp = pts.shp)
      plot.method <- "lines"
    }
    if (method[i] == "opoly") {
      shp <- PCA2shp_polynomials(
        pos = pos, rot = rot[ids, ], mshape = mshape[ids], amp.shp = amp.shp, pts.shp = pts.shp,
        ortho = TRUE, baseline1 = PCA$baseline1[1:2 +
                                                  (i - 1) * 2], baseline2 = PCA$baseline2[1:2 +
                                                                                            (i - 1) * 2]
      )
      plot.method <- "lines"
    }
    if (method[i] == "npoly") {
      shp <- PCA2shp_polynomials(
        pos = pos, rot = rot[ids, ], mshape = mshape[ids], amp.shp = amp.shp, pts.shp = pts.shp,
        ortho = FALSE, baseline1 = PCA$baseline1[1:2 +
                                                   (i - 1) * 2], baseline2 = PCA$baseline2[1:2 +
                                                                                             (i - 1) * 2]
      )
      plot.method <- "lines"
    }
    if (method[i] == "procrustes") {
      shp <- PCA2shp_procrustes(pos = pos, rot = rot[ids, ], mshape = mshape[ids], amp.shp = amp.shp)
      plot.method <- "points"
    }
    shp <- lapply(shp, coo_template, size = size.shp.final[i])
    shp <- lapply(shp, coo_center)
    shp <- lapply(shp, coo_rotate, rotate.shp[i])
    if (flipx.shp[i]) {
      shp <- lapply(shp, coo_flipx)
    }
    if (flipy.shp[i]) {
      shp <- lapply(shp, coo_flipy)
    }

    for (s in 1:length(shp)) {
      shp[[s]] <- coo_trans(
        shp[[s]], pos[s, 1] + dx[i],
        pos[s, 2] + dy[i]
      )
    }
  }
  list(shp = shp, pos = pos)
}

x <- bot %>%
  efourier(6) %>%
  PCA()
df <- x %>% as_df(retain = 2)

shps <- pre_geom_morphospace(x, 1, 2, pos.shp = "range")
shps_ready <- map(seq_along(shps$shp),
                ~shps$shp[[.x]] %>% coo_template(size=0.01) %>%
                  coo_trans(x=shps$pos[.x, 1], y=shps$pos[.x, 2])
)

shps_df <- map_df(seq_along(shps_ready),
                  ~ tibble(i = .x,
                           x = shps_ready[[.x]][, 1],
                           y = shps_ready[[.x]][, 2]))

gg <- ggplot(df) +
  aes(PC1, PC2) +
  geom_point() +
  coord_equal()

gg + geom_path(aes(x = x, y = y, group = i), data = shps_df)


