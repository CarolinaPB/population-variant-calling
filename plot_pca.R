library(data.table)
library("RColorBrewer")
library(ggplot2)

# create color pallete
c26 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown","navy"
)

# Read in files
sample_details <- fread("Sample_sheet_Fasil_Chicken_Blood_Samples.txt")
sample_details$no <- sub("^", "a", sample_details$no )

pca <- fread("chicken_old_samples.eigenvec", header = F)
eigenval <- scan("chicken_old_samples.eigenval")
pca <- pca[,-1]
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))



# Plot PC1-PC2 - no new samples
# pdf(file="PCA_chicken.pdf")
merged <- merge(pca[,1:5], sample_details, by.x="ind", by.y ="no", all.x = T)

## elevation (discrete color)
merged[ind=="a130AK"]$elevation <- sample_details[no=="a130"]$elevation
merged$color <- c26[as.factor(merged$elevation)]
plot(merged[,.(PC1,PC2)], pch=19, cex=0.4, xlab="", ylab="", col=merged$color)
legend("bottomleft", title = "elevation",col=unique(merged$color), legend=unique(merged$elevation), pch=19, cex = 0.6)
title(main="African chicken", xlab=paste0("PC1 (", signif(pve$pve[1], 3), "%)"), ylab=paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

## elevation color gradient
b <- ggplot(merged, aes(PC1, PC2, col = elevation)) + geom_point(size = 0.8)
# b <- b + scale_colour_manual(values = c26)
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

## village code
merged[ind=="a130AK"]$village_code <- sample_details[no=="a130"]$village_code
merged$color <- c26[as.factor(merged$village_code)]
plot(merged[,.(PC1,PC2)], pch=19, cex=0.4, xlab="", ylab="", col=merged$color)
legend("bottomleft", title = "Village",col=unique(merged$color), legend=unique(merged$village_code), pch=19, cex = 0.6)
title(main="African chicken", xlab=paste0("PC1 (", signif(pve$pve[1], 3), "%)"), ylab=paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

## district
merged[ind=="a130AK"]$district <- sample_details[no=="a130"]$district
merged$color <- c26[as.factor(merged$district)]
plot(merged[,.(PC1,PC2)], pch=19, cex=0.4, xlab="", ylab="", col=merged$color)
legend("bottomleft", title = "district",col=unique(merged$color), legend=unique(merged$district), pch=19, cex = 0.6)
title(main="African chicken", xlab=paste0("PC1 (", signif(pve$pve[1], 3), "%)"), ylab=paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

## region
merged[ind=="a130AK"]$region <- sample_details[no=="a130"]$region
merged$color <- c26[as.factor(merged$region)]
plot(merged[,.(PC1,PC2)], pch=19, cex=0.4, xlab="", ylab="", col=merged$color)
legend("bottomleft", title = "region",col=unique(merged$color), legend=unique(merged$region), pch=19, cex = 0.6)
title(main="African chicken", xlab=paste0("PC1 (", signif(pve$pve[1], 3), "%)"), ylab=paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

# dev.off()

# Plot variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# Plot PC2-PC3
# pdf(file="PCA_chicken_PC2_PC3.pdf")

# elevation (discrete color)
merged[ind=="a130AK"]$elevation <- sample_details[no=="a130"]$elevation
merged$color <- c26[as.factor(merged$elevation)]
plot(merged[,.(PC2,PC3)], pch=19, cex=0.4, xlab="", ylab="", col=merged$color)
legend("bottomleft", title = "elevation",col=unique(merged$color), legend=unique(merged$elevation), pch=19, cex = 0.6)
title(main="African chicken", xlab=paste0("PC2 (", signif(pve$pve[2], 3), "%)"), ylab=paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

# elevation (color gradient)
b <- ggplot(merged, aes(PC2, PC3, col = elevation)) + geom_point(size = 0.8)
# b <- b + scale_colour_manual(values = c26)
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

# village code
merged[ind=="a130AK"]$village_code <- sample_details[no=="a130"]$village_code
merged$color <- c26[as.factor(merged$village_code)]
plot(merged[,.(PC2,PC3)], pch=19, cex=0.4, xlab="", ylab="", col=merged$color)
legend("bottomleft", title = "Village",col=unique(merged$color), legend=unique(merged$village_code), pch=19, cex = 0.6)
title(main="African chicken", xlab=paste0("PC2 (", signif(pve$pve[2], 3), "%)"), ylab=paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

# district
merged[ind=="a130AK"]$district <- sample_details[no=="a130"]$district
merged$color <- c26[as.factor(merged$district)]
plot(merged[,.(PC2,PC3)], pch=19, cex=0.4, xlab="", ylab="", col=merged$color)
legend("bottomleft", title = "district",col=unique(merged$color), legend=unique(merged$district), pch=19, cex = 0.6)
title(main="African chicken", xlab=paste0("PC2 (", signif(pve$pve[2], 3), "%)"), ylab=paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

# region
merged[ind=="a130AK"]$region <- sample_details[no=="a130"]$region
merged$color <- c26[as.factor(merged$region)]
plot(merged[,.(PC2,PC3)], pch=19, cex=0.4, xlab="", ylab="", col=merged$color)
legend("bottomleft", title = "region",col=unique(merged$color), legend=unique(merged$region), pch=19, cex = 0.6)
title(main="African chicken", xlab=paste0("PC2 (", signif(pve$pve[2], 3), "%)"), ylab=paste0("PC3 (", signif(pve$pve[3], 3), "%)"))


# dev.off()

# GGplots 
# b <- ggplot(merged, aes(PC1, PC2, col = elevation)) + geom_point(size = 0.8)
# # b <- b + scale_colour_manual(values = c26)
# b <- b + coord_equal() + theme_light()
# b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
# 
# 
# b <- ggplot(merged, aes(PC2, PC3, col = elevation)) + geom_point(size = 0.8)
# # b <- b + scale_colour_manual(values = c26)
# b <- b + coord_equal() + theme_light()
# b + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))
# 
# b <- ggplot(merged, aes(PC3, PC4, col = elevation)) + geom_point(size = 0.8)
# # b <- b + scale_colour_manual(values = c26)
# b <- b + coord_equal() + theme_light()
# b + xlab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) + ylab(paste0("PC4 (", signif(pve$pve[4], 3), "%)"))



