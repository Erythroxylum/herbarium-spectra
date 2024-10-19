################################################################################
# Plot confusion matrix
################################################################################

plotConfusionMatrix = function(matrix, directory, fileName) {
  library(corrplot)
  library(naniar)
  
  cols = colorRampPalette(c('white', '#fe9929'))
  
  jpeg(filename = paste(directory, fileName, sep = "/"),
       width = 12, height = 12, units = 'in', res = 1200)
  par(mfrow = c(1,1))
  corrplot::corrplot(matrix,
                     cl.pos = 'n',
                     method = 'square',
                     tl.col = 'black',
                     na.label = 'square',
                     na.label.col = 'white',
                     addCoef.col = '#542788',
                     number.digits = 2,
                     number.cex = .7,
                     col = cols(10))
  dev.off()
}

saveRDS(plotConfusionMatrix, "functions/plotConfusionMatrix.rds")
