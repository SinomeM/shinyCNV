
# Copied and edited from CNValidatron, ideally when the package is public
# the two versions should be merged back


plot_cnv <- function(cnv, samp, snps = NULL, adjusted_lrr = T,
                     min_lrr = -1.4, max_lrr = 1.3,
                     # the following parameters should not be changed by most users
                     shrink_lrr = 0.1, w = 96, z = 4, k1 = 31, k2 = 26,
                     l_wind = 20000000, # top row Mbp
                     mx_lr = 2) {  # top row LRR range

  # w k1, k2, z and in_out_ratio are fixed for the moment


  # everything will be [0,w-1] then I will add 1 to make it [1,w]
  w <- w-1

  # in_out_ratio is now a function of the length
  len <- cnv$end - cnv$start + 1
  if ('locus' %in% colnames(cnv)) len <- cnv$loc_en - cnv$loc_st + 1

  if (len <= 100000) in_out_ratio <- 9
  if (between(len, 100001, 1000000)) in_out_ratio <- 7
  if (len > 1000000) in_out_ratio <- 5

  ss <- cnv$start - (in_out_ratio*len);  ee <- cnv$end + (in_out_ratio*len)

  if ('locus' %in% colnames(cnv)) {
    ss <- cnv$loc_st - (in_out_ratio*len);  ee <- cnv$loc_en + (in_out_ratio*len)
  }


  # load snps data, ALL chromosome is loaded now!
  dt <- load_snps_tbx(cnv, samp, snps, in_out_ratio, adjusted_lrr, min_lrr, max_lrr, shrink_lrr)
  # keep the full chromsome for the third row of the png
  dt_big <- dt[[2]]
  dt <- dt[[1]]
  if (nrow(dt) == 0) {
    warning('Empty tabix, no image generated for sample', samp$sample_ID)
    return(data.table())
  }

  # each point need to be used for both lrr and baf so dt must be duplicated
  dt_lrr <- copy(dt)
  dt_baf <- dt

  # top row
  #  the top row is zoomed out by a factor of 3
  ss2 <- cnv$start - (in_out_ratio*len*3);  ee2 <- cnv$end + (in_out_ratio*len*3)
  # if it is not at least 12.5Mbp then force it to 12.5Mbp
  len_diff <- (12500000 - (ee2 - ss2 + 1)) / 2
  if (len_diff > 0) {
    ss2 <- ss2 - len_diff
    ee2 <- ee2 + len_diff
  }

  dt_big <- dt_big[between(position, ss2, ee2),]
  dt_big[, x := round(((position-ss2)/(ee2-ss2)) * w)]
  dt_big[, y := round(((lrr-(-mx_lr))/(mx_lr-(-mx_lr))) * k2) + (k1*2 + z*2)]
  dt_big[, ':=' (x = x+1, y = y+1)]

  w <- w+1

  brr <- seq(from = ss, to = ee, length.out = 20)
  brr_l <- round(brr/1000000, 1)
  brr_l <- ifelse(brr_l > 0, brr_l, 0)
  brr2 <- seq(from = ss2, to = ee2, length.out = 20)
  brr_l2 <- round(brr2/1000000, 1)
  brr_l2 <- ifelse(brr_l2 > 0, brr_l2, 0)


  dt_lrr[lrr <= min_lrr, lrr := min_lrr+0.01]
  dt_lrr[lrr >= max_lrr, lrr := max_lrr-0.01]
  dt_big[lrr <= min_lrr, lrr := min_lrr+0.01]
  dt_big[lrr >= max_lrr, lrr := max_lrr-0.01]

  # plot in the original space
    a <- ggplot(dt_lrr, aes(position, lrr)) + geom_point(alpha = 0.3, colour = 'red') +
           ylim(min_lrr, max_lrr) + theme_bw() +
           theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
           scale_x_continuous(breaks = brr, labels = brr_l, limits = c(ss, ee))
    if (cnv$GT == 1) a <- a + geom_segment(x = cnv$start, xend = cnv$end, y = 1,
                                           yend = 1, linetype = 3)
    if (cnv$GT == 2) a <- a + geom_segment(x = cnv$start, xend = cnv$end, y = -1,
                                           yend = -1, linetype = 3)

    b <- ggplot(dt_baf, aes(position, baf)) + geom_point(alpha = 0.3, colour = 'blue') +
           theme_bw() +
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(), axis.title.y = element_blank()) +
           scale_x_continuous(breaks = brr, limits = c(ss, ee))

    c <- ggplot(dt_big, aes(position, lrr)) + geom_point(alpha = 0.1, colour = 'purple') +
           ylim(min_lrr, max_lrr) + theme_bw() +
           theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
           scale_x_continuous(breaks = brr2, labels = brr_l2, limits = c(ss2, ee2))
    if (cnv$GT == 1) c <- c + geom_segment(x = cnv$start, xend = cnv$end, y = 1,
                                           yend = 1, linetype = 3)
    if (cnv$GT == 2) c <- c + geom_segment(x = cnv$start, xend = cnv$end, y = -1,
                                           yend = -1, linetype = 3)

    if ('locus' %in% colnames(cnv)) {
      a <- a + geom_segment(x = cnv$loc_st, xend = cnv$loc_st, y = -1.45, yend = -1.40) +
               geom_segment(x = cnv$loc_en, xend = cnv$loc_en, y = -1.45, yend = -1.40) +
               geom_segment(x = cnv$loc_st, xend = cnv$loc_st, y = 1.40, yend = 1.45) +
               geom_segment(x = cnv$loc_en, xend = cnv$loc_en, y = 1.40, yend = 1.45)
      c <- c + geom_segment(x = cnv$loc_st, xend = cnv$loc_st, y = -1.45, yend = -1.40) +
               geom_segment(x = cnv$loc_en, xend = cnv$loc_en, y = -1.45, yend = -1.40) +
               geom_segment(x = cnv$loc_st, xend = cnv$loc_st, y = 1.40, yend = 1.45) +
               geom_segment(x = cnv$loc_en, xend = cnv$loc_en, y = 1.40, yend = 1.45)

      b <- b + geom_segment(x = cnv$loc_st, xend = cnv$loc_st, y = -0.03, yend = -0.01) +
               geom_segment(x = cnv$loc_en, xend = cnv$loc_en, y = -0.03, yend = -0.01) +
               geom_segment(x = cnv$loc_st, xend = cnv$loc_st, y = 1.01, yend = 1.03) +
               geom_segment(x = cnv$loc_en, xend = cnv$loc_en, y = 1.01, yend = 1.03)
    }

    pl <- cowplot::plot_grid(c, b, a, ncol = 1)
    return(pl)
}



load_snps_tbx <- function(cnv, samp, snps = NULL, in_out_ratio = 1, adjusted_lrr = T,
                          min_lrr = -1.2, max_lrr = 1, shrink_lrr = NULL) {
  chr <- cnv$chr
  start <- cnv$start
  end <- cnv$end
  len <- end - start + 1
  tbx_path <- samp$file_path_tabix
  
  if ('locus' %in% colnames(cnv)) {
    start <- cnv$loc_st
    end <- cnv$loc_en
    len <- end - start + 1
  }

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
