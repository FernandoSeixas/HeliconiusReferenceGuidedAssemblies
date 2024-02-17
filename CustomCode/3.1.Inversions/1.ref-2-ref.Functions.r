####
maxDist = 20000


mergeAligns = function(x) {
  xMerge = data.frame()
  for (scaff in unique(x$Q.scaffold)) {
    s1 = subset(x, Q.scaffold == scaff)
    if (nrow(s1) > 1) {
      s1.order = s1[order(s1$c_Bsta),]
      s1.merge = s1.order[1,]
      # try to merge 
      for (rw in 2:nrow(s1.order)) {
        mrg.rw = nrow(s1.merge)
        # close enough & same strand
        if (s1.order$c_Bsta[rw] - s1.merge$c_Bend[mrg.rw] <= maxDist & s1.order$strand[rw] == s1.merge$strand[mrg.rw]) {
          ## update R.coordinates and C.coordinates
          s1.merge$R.end[mrg.rw] = s1.order$R.end[rw]
          s1.merge$c_Bend[mrg.rw] = s1.order$c_Bend[rw]
          ## update Q.coordinates
          if (s1.order$strand[rw] == "+") { s1.merge$Q.end[mrg.rw] = s1.order$Q.end[rw] }
          if (s1.order$strand[rw] == "-") { s1.merge$Q.start[mrg.rw] = s1.order$Q.start[rw] }
        }
        # not close enough or not in the same strand
        if (s1.order$c_Bsta[rw] - s1.merge$c_Bend[mrg.rw] > maxDist | s1.order$strand[rw] != s1.merge$strand[mrg.rw]) {
          s1.merge = rbind(s1.merge, s1.order[rw,])
        }
      }
      xMerge = rbind(xMerge, s1.merge)
    }
    if (nrow(s1) == 1) {
      xMerge = rbind(xMerge, s1)
    }
  }
  return(xMerge)
}

# x = subset(heralat_2_heradem_Minimap2, chrom == "Herato06")
# s1.mergeUP = updateRcoords(xMerge)
# ggplot(s1.mergeUP) +
#   geom_gene_arrow(aes(xmin = c_Bsta/1000000, xmax = c_Bend/1000000, y=1, color=Q.scaffold), fill=ifelse(xMerge$strand == "-", "#b2182b", "#f5f5f5"), lwd=0.6) +
#   facet_wrap(~Q.scaffold, nrow=7)


updateRcoords = function(x) {
  for (rw in 1:nrow(x)) {
    if (x$strand[rw] == "-") {
      st = x$c_Bend[rw]
      en = x$c_Bsta[rw]
      x$c_Bsta[rw] = st
      x$c_Bend[rw] = en
    }
  }
  return(x)
}
