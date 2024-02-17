mergemaps <- function(x) {
        name = genes[x] # name of gene
        genemaps = subset(inp, gene == name)
        # if gene only has 1 hit
        if (nrow(genemaps) == 1) {genemerge = genemaps}
        # if more than one hit 
        else if (nrow(genemaps) > 1) {
                genemerge = genemaps[1,] # start "annotation" for that gene
                count  = 1 
                for (rw in 2:nrow(genemaps)) {
                        # different scaffold [add new annotation]
                        if (genemaps$scaffold[rw] != genemerge$scaffold[count]) {
                                genemerge = rbind(genemerge, genemaps[rw,])
                                count = count + 1
                        }
                        # same scaffold
                        else if (genemaps$scaffold[rw] == genemerge$scaffold[count]) {
                                # if NOT overlapping or nearby [add new annotation]
                                if (genemaps$sta[rw] - genemerge$end[count] > buffer) {
                                        genemerge = rbind(genemerge, genemaps[rw,])
                                        count = count + 1
                                }
                                # if overlapping or nearby
                                else if (genemaps$sta[rw] - genemerge$end[count] <= buffer) {
                                        # if different strand [add new annotation]
                                        if (genemaps$strand[rw] != genemerge$strand[count]) {
                                                genemerge = rbind(genemerge, genemaps[rw,])
                                                count = count + 1
                                        }
                                        # if same strand
                                        else if (genemaps$strand[rw] == genemerge$strand[count]) {
                                                genemerge$end[count] = genemaps$end[rw]                                        # update end of mapping region
                                                genemerge$maplen[count] = genemerge$maplen[count] + genemaps$maplen[rw]        # update length of gene hits
                                                genemerge$identity[count] = genemerge$identity[count] + genemaps$identity[rw]  # update percent identity
                                                genemerge$hits[count] = genemerge$hits[count] + 1                              # update number of hits
                                        }
                                }
                        }
                }
        }
        return(genemerge)
}


