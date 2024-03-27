cnvs <- fread('data/cnvs.txt')
cnvs[, plot_id := paste0('/home/simone/Documents/shinyCNV/simpler_app/data/plots/',
                         sample_ID, '_', chr, '_', start, '.png')]
for (i in 1:cnvs[, .N])
  ggsave(cnvs[i, plot_id], ggplot(cnvs[i], aes(start,end, label = plot_id)) + 
                             geom_point() + geom_text(size = 1),
         width = 2, height = 2)

fwrite(cnvs, './data/cnvs_with_pl.txt', sep = '\t')
