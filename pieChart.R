#slices <- c(48.0, 48.0, 0.5, 3.4) 
#lbls <- c("Correlate1-1", "Orig_Group", "random2-1", "random3-1")
outfilename3 <- paste('AA', inFilename3, 'PieChart.png', sep="-")
png(file = outfilename3)
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Significant Sites by Metadata Categories")
dev.off()