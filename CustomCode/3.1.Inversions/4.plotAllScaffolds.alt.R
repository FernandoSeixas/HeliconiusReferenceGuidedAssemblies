
## get scaffolds within specified region
scaffoldRegions = subset(minimaps, R.scaffold == as.character(sc))

## adjust start and end depending on strand
for (rw in 1:nrow(scaffoldRegions)) {
  if (scaffoldRegions$strand[rw] == "-") {
    tmp = scaffoldRegions$R.start[rw]
    scaffoldRegions$R.start[rw] = scaffoldRegions$R.end[rw]
    scaffoldRegions$R.end[rw] = tmp
  }
}

## define species order
scaffoldRegions$species = factor(scaffoldRegions$species, levels = c(mel_clade,mid_clade,era_clade,refs))

# define if scaffold maps to 1 (normal) or 2 (inversion) strands
scaffoldRegions$type = 1
scaffoldRegions$span = -1
for (spp in levels(scaffoldRegions$species) ) {
  print(spp)
  s1 = subset(scaffoldRegions, species == spp)
  for (contig in unique(s1$Q.scaffold)) {
    s2 = subset(s1, Q.scaffold == as.character(contig))
    if (length(unique(s2$strand)) == 2) {
      scaffoldRegions$span = ifelse(scaffoldRegions$species == spp & scaffoldRegions$Q.scaffold == as.character(contig), abs(min(s2$R.start,s2$R.end)-max(s2$R.start,s2$R.end)), scaffoldRegions$span)
      scaffoldRegions$type = ifelse(scaffoldRegions$species == spp & scaffoldRegions$Q.scaffold == as.character(contig), 2, scaffoldRegions$type)
      # print(s2)
    }
  }
}
scaffoldRegions$type = as.character(scaffoldRegions$type)

## plot
ggplot(subset(scaffoldRegions, type==2)) +
  # add breakpoints
  geom_vline(xintercept = bk1/1000000, col="red") +
  geom_vline(xintercept = bk2/1000000, col="red") +
  # plot DISCOVAR scaffolds
  geom_tile(aes(x=(R.start+R.end)/2000000, y=1, width=(R.start-R.end)/1000000, height=0.8, fill=strand, color=strand), alpha=0.5) +
  # geom_gene_arrow(aes(xmin = R.start/1000000, xmax = R.end/1000000, y = 1, fill=strand, color=strand)) +
  scale_fill_manual(values = c("#fdae61","#4e7dff")) +
  scale_color_manual(values = c("#fdae61","#4e7dff")) +
  # visual
  facet_grid(species ~ ., drop = TRUE, scales = "free_y") +
  ggtitle(paste0(ref," // ", sc,":",bk1,"-",bk2)) +
  theme(panel.grid = element_blank())
