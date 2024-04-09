
# Copied and edited from CNValidatron, ideally when the package is public
# the two versions should be merged back



plot_cnv <- function(cnv, samp, snps = NULL, adjusted_lrr = T,
                     tmp_plot = 0, min_lrr = -1.4, max_lrr = 1.3,
                     # the following parameters should not be changed by most users
                     shrink_lrr = 0.1, w = 96, z = 4, k1 = 31, k2 = 26, in_out_ratio = 6,
                     l_wind = 20000000, # top row Mbp
                     mx_lr = 2.5,  # top row LRR range
                     norm_wind_size = 6) { # top row normalisation windows size

  # w k1, k2, z and in_out_ratio are fixed for the moment


  # everything will be [0,w-1] then I will add 1 to make it [1,w]
  w <- w-1

  # load snps data, ALL chromosome is loaded now!
  dt <- load_snps_tbx(cnv, samp, snps, in_out_ratio, adjusted_lrr, min_lrr, max_lrr, shrink_lrr)
  # keep the full chromsome for the third row of the png
  dt_big <- dt[[2]]
  dt <- dt[[1]]
  if (nrow(dt) == 0) {
    warning('Empty tabix, no image generated for sample', samp$sample_ID)
    return(data.table())
  }

  # bottom and middle row, select the relevant points and
  # move position to the x coordinates in the new system
  len <- cnv$end - cnv$start + 1
  ss <- cnv$start - (in_out_ratio*len);  ee <- cnv$end + (in_out_ratio*len)
  dt <- dt[between(position, ss, ee), ]
  dt[, x := round(((position-ss)/(ee-ss)) * w)]

  # each point need to be used for both lrr and baf so dt must be duplicated
  # formally the second copy is not needed
  dt_lrr <- copy(dt)
  dt_baf <- copy(dt)

  # move lrr and baf on hte y coordinates in the new system
  dt_lrr[, y := round(((lrr-(min_lrr))/(max_lrr-(min_lrr))) * k1)]
  dt_baf[, y := round(((baf-0)/(1-0)) * k1) + k1 + z]

  # pixel coordinates must be > 0
  dt_baf[, ':=' (x = x+1, y = y+1)]
  dt_lrr[, ':=' (x = x+1, y = y+1)]

  # top row
  center <- cnv$start + len/2
  a <- center - l_wind/2
  b <- center + l_wind/2
  # if a is out of bound then the 30Mbp region must start from the first SNP
  if (a < dt_big[, min(position)]) {
    a <- dt_big[, min(position)]
    b <- a+l_wind
  }
  else
  # same thing on the other side, both cannot be true at the same time in a human genome
    if (b > dt_big[, max(position)]) {
      b <- dt_big[, max(position)]
      a <- b-l_wind
    }

  dt_big <- dt_big[between(position, a, b),]
  dt_big[, x := round(((position-a)/(b-a)) * w)]
  dt_big[, y := round(((lrr-(-mx_lr))/(mx_lr-(-mx_lr))) * k2) + (k1*2 + z*2)]
  dt_big[, ':=' (x = x+1, y = y+1)]

  w <- w+1

  # plot in the original space
  if (tmp_plot == 1) {
    a <- ggplot(dt_lrr, aes(position, lrr)) + geom_point(alpha = 0.3, colour = 'red') +
           ylim(min_lrr, max_lrr) + theme_bw() + xlim(ss, ee) +
           geom_segment(x = cnv$start, xend = cnv$end, y = 0, yend = 0, linetype = 3) +
           theme(axis.title.x = element_blank(), axis.title.y = element_blank())
    b <- ggplot(dt_baf, aes(position, baf)) + geom_point(alpha = 0.3, colour = 'blue') +
           theme_bw() + xlim(ss, ee) +
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(), axis.title.y = element_blank())
    c <- ggplot(dt_big, aes(position, lrr)) + geom_point(alpha = 0.1, colour = 'purple') +
           geom_segment(x = cnv$start, xend = cnv$end, y = 0, yend = 0) +
           ylim(min_lrr, max_lrr) + theme_bw() +
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(), axis.title.y = element_blank())
    return(cowplot::plot_grid(c, b, a, ncol = 1))
  }
#
#  # plot in the new coordinates but suitable for humans
#  if (tmp_plot == 2) {
#    a <- ggplot(dt_lrr, aes(x, y)) + geom_point(alpha = 0.1, colour = 'red') +
#           xlim(0, w) + ylim(0+1, k1+1) + theme_bw() +
#           theme(axis.title.x = element_blank(), axis.title.y = element_blank())
#    b <- ggplot(dt_baf, aes(x, y)) + geom_point(alpha = 0.1, colour = 'blue') +
#           xlim(0, w) + ylim(k1+z + 1, k1*2 + z + 1) + theme_bw() +
#           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#                 axis.title.x = element_blank(), axis.title.y = element_blank())
#    c <- ggplot(dt_big, aes(x, y)) + geom_point(alpha = 0.05, colour = 'purple') +
#           xlim(0, w) + ylim((k1*2)+(z*2) + 1, w + 1) + theme_bw() +
#           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#                 axis.title.x = element_blank(), axis.title.y = element_blank())
#    return(cowplot::plot_grid(c, b, a, ncol = 1))
#  }
#
#  # create the pixel map, and normalise values between 0 and 1
#  dt_lrr <- dt_lrr[, .N, by = c('x', 'y')][, N := N/max(N)]
#  dt_baf <- dt_baf[, .N, by = c('x', 'y')][, N := N/max(N)]
#  # the top row normalisation is performed in windows first, then on all
#  dt_big <- dt_big[, .N, by = c('x', 'y')]
#  dt_big[, wind1 := x %% norm_wind_size]
#  dt_big[, wind2 := x %% (norm_wind_size+2)]
#  dt_big[, n := N/max(N), by = wind1]
#  dt_big[, n := N/max(N), by = wind2]
#  dt_big[, N := n][, n := NULL][, N := N/max(N)]
#
#  # create the full picture
#  dt <- as.data.table(rbind(t(combn(1:w, 2)), t((combn(w:1, 2))))) # almost all combinations
#  colnames(dt) <- c('x', 'y')
#  dt <- rbind(dt, data.table(x = 1:w, y = 1:w)) # the diagonal was missing
#  # fill in the values
#  dt <- merge(dt, dt_lrr, by = c('x', 'y'), all.x = T)
#  setnames(dt, 'N', 'a')
#  dt <- merge(dt, dt_baf, by = c('x', 'y'), all.x = T)
#  setnames(dt, 'N', 'b')
#  dt <- merge(dt, dt_big, by = c('x', 'y'), all.x = T)
#  setnames(dt, 'N', 'c')
#  dt[!is.na(a), value := a]
#  dt[!is.na(b), value := b]
#  dt[!is.na(c), value := c]
#  # dt[!is.na(a) & !is.na(b), value := a+b] # this should not be a case
#  dt[is.na(a) & is.na(b) & is.na(c), value := 0]
#  dt <- dt[, .(x, y, value)]
#
#  # plot almost identical to the final PNG
#  if (tmp_plot == 3)
#    return(ggplot(dt, aes(x, y, fill = value)) + geom_tile() + theme_minimal() +
#             scale_fill_gradient(low="white", high="black"))
#
#  dt[, y := abs(y-(max(y)+1))] # to deal with how imager use the y axis
#
#  return(dt)
}



load_snps_tbx <- function(cnv, samp, snps = NULL, in_out_ratio = 1, adjusted_lrr = T,
                          min_lrr = -1.2, max_lrr = 1, shrink_lrr = NULL) {
  chr <- cnv$chr
  start <- cnv$start
  end <- cnv$end
  len <- end - start + 1
  tbx_path <- samp$file_path_tabix

  st <- start - (in_out_ratio*len)
  st <- ifelse(st < 0, 0, st)

  # load the whole chromosome now
  dt <- fread(cmd = paste0("tabix ", tbx_path, " ", chr, ":", 0,
                          "-", 1000000000), header = F)

  if (nrow(dt) == 0) {
    warning('File: ', tbx_path, ' seems empty or broken\n')
    return(data.table())
  }

  # gdk-ipsych case
  #if (ncol(dt) == 7) colnames(dt) <- c("chr", "position", "end", "LRR", "LRRadj", "BAF", "snp")
  #else {
  if (adjusted_lrr) colnames(dt) <- c("chr", "position", "end", "LRR", "BAF", "LRRadj")
  else colnames(dt) <- c("chr", "position", "end", "LRR", "BAF")
  #}

    # filter SNPs if snp object is provided
  if (!is.null(snps))
    dt <- dt[paste0(chr, position) %in% snps[, paste0(Chr, Position)], ]

  setorder(dt, position)

  if (adjusted_lrr) setnames(dt, c('BAF', 'LRRadj'), c('baf', 'lrr'))
  else setnames(dt, c('BAF', 'LRR'), c('baf', 'lrr'))

  # save the whole chromosome unprocessed for the first row
  dt <- dt[, .(chr, position, lrr, baf)]
  dt_whole <- copy(dt)
  dt <- dt[between(position, start - len*in_out_ratio, end + len*in_out_ratio), ]

  ## some preprocessing ##
  # restrict the lrr space
  dt[lrr > max_lrr, lrr := max_lrr][lrr < min_lrr, lrr := min_lrr]
  # if lrr or baf is missing exclude the point
  dt <- dt[!(is.na(lrr) | is.na(baf)),]

  # compute mean and SD in these three groups, could be simplified using dt[,,by]
  dt[position < start, group := 1][
     between(position, start, end), group := 2][position > end, group := 3]
  ms1 <- dt[group == 1, c(mean(lrr, na.rm = T), sd(lrr, na.rm = T))]
  ms2 <- dt[group == 2, c(mean(lrr, na.rm = T), sd(lrr, na.rm = T))]
  ms3 <- dt[group == 3, c(mean(lrr, na.rm = T), sd(lrr, na.rm = T))]

  if (!is.null(shrink_lrr)) {
    # snps in each group get pulled towards the group mean proportionally
    # to their distance and shrink_lrr
    dt[group == 1 & lrr > ms1[1], lrr := lrr - abs(lrr-ms1[1])*shrink_lrr][
         group == 1 & lrr < ms1[1], lrr := lrr + abs(ms1[1]-lrr)*shrink_lrr]

    dt[group == 2 & lrr > ms2[1], lrr := lrr - abs(lrr-ms2[1])*shrink_lrr][
         group == 2 & lrr < ms2[1], lrr := lrr + abs(ms2[1]-lrr)*shrink_lrr]

    dt[group == 3 & lrr > ms3[1], lrr := lrr - abs(lrr-ms3[1])*shrink_lrr][
         group == 3 & lrr < ms3[1], lrr := lrr + abs(ms3[1]-lrr)*shrink_lrr]
  }

  # ouliers removal, 3SDs. Should not do anything after the rest
  dt[group == 1 & !between(lrr, ms1[1]-3*ms1[2], ms1[1]+3*ms1[2]), lrr := NA]
  dt[group == 2 & !between(lrr, ms2[1]-3*ms2[2], ms2[1]+3*ms2[2]), lrr := NA]
  dt[group == 3 & !between(lrr, ms3[1]-3*ms3[2], ms3[1]+3*ms3[2]), lrr := NA]
  dt <- dt[!is.na(lrr), ]


  return(list(dt, dt_whole))

}
