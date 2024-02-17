#################### DETECTING INVERSIONS ####################

## function to merge regions ========================
mergeRegions <- function(x) {
  # sort by reference coordinates and separate by strand
  x.sort = x[with(x, order(R.start)),]
  x.merge = x.sort[1,]
  # try to merge/concatenate same strand alignments
  for (rw in 2:nrow(x.sort)) {
    # next alignment close enough (in terms of reference coordinates)
    if ( x.sort$R.start[rw] <= (x.merge$R.end[nrow(x.merge)] + maxGap) ) {
      ## guarantee the new alignment is not COMPLETELY within the first (Query coordinates)
      # determine coordinates
      x1 = x.merge$Q.start[nrow(x.merge)]; x2 = x.merge$Q.end[nrow(x.merge)]; lenx = x2-x1+1
      y1 = x.sort$Q.start[rw]; y2 = x.sort$Q.end[rw]; leny = y2-y1+1
      ov = 0
      if (x1 <= y2 & y1 <= x2) {
        overlap = min(x2, y2)-max(x1, y1)+1 # length of overlap
        if (lenx == overlap | leny == overlap) { ov=1 } # one segment completely within the other
      }
      ## update Q coordinates
      if (x.sort$Q.start[rw] >= x.merge$Q.start[nrow(x.merge)]) { x.merge$Q.end[nrow(x.merge)] = x.sort$Q.end[rw] }
      if (x.sort$Q.start[rw] <= x.merge$Q.start[nrow(x.merge)]) { x.merge$Q.start[nrow(x.merge)] = x.sort$Q.start[rw] }
      ## update R coordinates
      x.merge$R.end[nrow(x.merge)] = x.sort$R.end[rw]
      x.merge$M.Bases[nrow(x.merge)] = x.merge$M.Bases[nrow(x.merge)] + x.sort$M.Bases[rw]
      x.merge$Nb.Bases[nrow(x.merge)] = x.merge$Nb.Bases[nrow(x.merge)] + x.sort$Nb.Bases[rw]
      x.merge$dv[nrow(x.merge)] = min(x.merge$dv[nrow(x.merge)], x.sort$dv[rw])
    }
    # if not close enough add new segment
    if ( x.sort$R.start[rw] > (x.merge$R.end[nrow(x.merge)]+maxGap) ) {
      x.merge = rbind(x.merge, x.sort[rw,])
    }
  }
  return(x.merge)
}



## function that parses through scaffolds, filters alignments [minAlign, DV, MQ], =======
## determines if they are potentially informative for inversions ========================
ScaffoldParser <- function(x) {
  # determine scaffold categories dataframes
  dfTemp = data.frame()
  # define species name
  spp = str_replace(x, paste0(dir, subdir,"/"), "")
  spp = str_replace(spp, paste0("-2-",ref,".filter.mm2"), "")
  # spp = str_replace(x, paste0(dir, "2.minimap2/vGlued"), "")
  # spp = str_replace(spp, paste0("-2-",ref,"Glued.filter.mm2"), "")
  # print(spp)
  
  ## read data -------------------------
  inp=read.table(x)
  names(inp)=c("Q.scaffold","Q.length","Q.start","Q.end","strand","R.scaffold","R.length","R.start","R.end","M.Bases","Nb.Bases","MQ","alignType","cm","s1","s2","dv")
  
  ## filter alignments [based on length mapped to reference, divergence to reference and mapping quality] ---
  # filter alignments shorter than "minAlign"
  inp = subset(inp, Nb.Bases >= minAlign)
  # filter alignments that are too divergent
  inp$dv = as.numeric(str_replace(inp$dv, "dv:f:", ""))
  inp = subset(inp, dv < 0.25)
  # filter low MQ alignments
  inp = subset(inp, MQ >= minMQ)
  
  ## determine if alignments are potentially informative about inversions -----------------
  df1 = data.frame()
  df2 = data.frame()
  for (contig in unique(inp$Q.scaffold)) {
    # chose contig
    sub1 = droplevels(subset(inp, Q.scaffold == contig))
    # analyse mapping to one chromosome at time (if it maps to more than 1)
    for (chrom in unique(sub1$R.scaffold)) {
      sub2 = subset(sub1, R.scaffold == chrom)
      # retain if multiple mappings to same chromosome
      if (nrow(sub2) > 1) {
        # retain only if two different strands present
        plusCount =  nrow(subset(sub2, strand == "+"))
        minusCount = nrow(subset(sub2, strand == "-"))
        if (plusCount >= 1 & minusCount >= 1) {
          df1 = rbind(df1, sub2) # write to dataframe
        }
        else { df2 = rbind(df2, sub2) }
      }
      # if maps to one chromosome only but not candidate for inversion
      else { df2 = rbind(df2, sub2) }
    }
  }
  if (nrow(df1) > 0 ) {df1$class = "Pass"; df1$species = spp; dfTemp = rbind(dfTemp, df1)} # add the candidates to common dataframe
  if (nrow(df2) > 0 ) {df2$class = "Fail"; df2$species = spp; dfTemp = rbind(dfTemp, df2)} # add the candidates to common dataframe
  # print(length(unique(df1$Q.scaffold)))
  return(dfTemp)
}
