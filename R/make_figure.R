##########################################################################
#         Function to Create Figures Used in the Manuscript              #
##########################################################################


figure1 <- function()
{
  RColorBrewer::palette(RcolorBrewer::brewer.pal(n = 8, name = "Set2"))
  dir = "../output/fig1_data/"
  grDevices::pdf("../output/figures/figure1.pdf", width=12, height=4)
  graphics::par(mfrow=c(1,3))
  data = utils::read.csv(paste(dir, 'harmonized_umap.csv', sep=''), row.names=1)
  train = utlis::read.csv(paste(dir, 'training.csv', sep=''))
  test = utils::read.csv(paste(dir, 'test.csv', sep=''))
  plot(data, col="grey")
  #   plot(data, col="grey")
  #plot(data[train$cell_id,], col="cornflowerblue")
  plot(data[train$cell_id,], col=as.factor(train$celltype))
  #   plot(data, col="grey")
  #plot(data[test$cell_id,], col="orange")
  plot(data[test$cell_id,], col=as.factor(test$celltype))
  utils::dev.off()
}

figure2 <- function(dir_path, pdf_file)
{
  files = c("hvg_result.csv", "pca_result.csv", "harmonized_result.csv", "harmonized_umap_result.csv")

  dir = paste("../output/harmony/ffnn/", dir_path, "ffnn_", sep="")
  pdf(paste("../output/figures/", pdf_file, sep=""), width=12, height=9)
  par(mfrow=c(3,4))

  for (file in files)
  {
    data = read.csv(paste(dir, file, sep=""), row.names=1)
    names = c(1: length(colnames(data)))
    boxplot(data, names=names, horizontal=FALSE, ylim=c(0, 1), xlab="datasets", ylab="FFNN classification accuracy", col="#1B9E77")
    expected = 1/length(colnames(data))
    abline(h=expected, col="red")
  }

  dir = paste("../output/harmony/xgb/", dir_path, "xgb_", sep="")

  for (file in files)
  {
    data = read.csv(paste(dir, file, sep=""), row.names=1)
    names = c(1: length(colnames(data)))
    boxplot(data, names=names, horizontal=FALSE, xlab="datasets", ylim=c(0, 1),  ylab="XGB classification accuracy", col="#7570B3")
    expected = 1/length(colnames(data))
    abline(h=expected, col="red")
  }

  dir = paste("../output/harmony/knn/", dir_path, "knn_", sep="")

  for (file in files)
  {
    data = read.csv(paste(dir, file, sep=""), row.names=1)
    names = c(1: length(colnames(data)))
    boxplot(data, names=names, horizontal=FALSE, xlab="datasets", ylim=c(0, 1), ylab="kNN classification accuracy", col="#E7298A")
    expected = 1/length(colnames(data))
    abline(h=expected, col="red")
  }

  dev.off()
}

figure3 <- function()
{
  pdf("../output/figures/figure3.pdf", width=20, height=8)
  par(mfrow=c(2,3))

  dir = "../output/tables/"
  files = c("pbmcsca_ffnn.txt", "pbmcsca_xgb.txt", "pbmcsca_knn.txt", "panc8_ffnn.txt", "panc8_xgb.txt", "panc8_knn.txt")

  for (file in files)
  {
    data = read.delim(paste(dir, file, sep=""), sep = "\t", row.names = 1)
    barplot(t(as.matrix(data)), beside=TRUE, ylim=c(0, 1.0), ylab = "Classification accuracy",
            col = brewer.pal(n = length(colnames(data)), name = "Paired"))
    expected = 1/length(colnames(data))
    abline(h=expected, col="red")
  }
  dev.off()
}

figure4 <- function(dir_path, pdf_file)
{

  file = paste("../output/", dir_path, sep="")
  pdf(paste("../output/figures/", pdf_file, sep=""), width=12, height=9)
  par(mfrow=c(1,3))
  data <- read.delim(file, sep="\t")
  data[,1] <- NULL
  data <- as.matrix(data)
  for (i in c(1, 22, 43))
  {
    end = i + 20
    d <- rdist::rdist(data[, i:end])
    heatmap(as.matrix(d), col = colorRampPalette(brewer.pal(8, "Set2"))(3*256), margins=c(20, 20))
  }
  dev.off()
}

figure5 <- function(pdf_file)
{
  files = c("harmony1_21.csv", "fastmnn1_21.csv", "cca1_21.csv", "sctransform1_21.csv")
  pdf(paste("../output/figures/", pdf_file, sep=""), width=12, height=6)
  par(mfrow=c(2,4))
  dir = "../output/"
  colors = c("#1B9E77", "#1B9E77", "#7570B3", "#7570B3", "#E7298A", "#E7298A")
  labels = c("HCC1395", "HCC1395BL", "HCC1395", "HCC1395BL", "HCC1395", "HCC1395BL")
  for (file in files)
  {
    data <- read.csv(paste(dir, file, sep=""), row.names=1)
    boxplot(data, ylim=c(0.4, 1),  col=colors, names=labels, ylab="Classification Accuracy", xlab = "Cell Line")
    abline(h=0.5, col="red", lty=2)
  }
  colors = brewer.pal(n = 8, name = "Set2")
  names = c(1: length(rownames(data)))
  i = 5
  for (file in files)
  {
    data <- read.csv(paste(dir, file, sep=""), row.names=1)
    boxplot(t(data), ylim=c(0.4, 1), names=names, col=colors[i], ylab="Classification Accuracy", xlab="Paired Datasets")
    abline(h=0.5, col="red", lty=2)
    i = i + 1
  }
  dev.off()
}

figure6 <- function(pdf_file)
{
  colors = brewer.pal(n = 8, name = "Set2")
  files = c("scenario5_fastmnn_size.csv", "scenario5_cca_size.csv", "scenario5_sctransform_size.csv")
  pdf(paste("../output/figures/", pdf_file, sep=""), width=12, height=6)
  par(mfrow=c(1,2))
  dir = "../output/"
  data <- read.csv(paste(dir, "scenario5_harmony_size.csv", sep=""))
  data$AVE <- c(rowMeans(data[,3:5]))
  i = 5
  plot(data$Size, data$AVE, ylim=c(0.4,1), type="o", lwd=4, col=colors[i], ylab="Classification Accuracy", xlab="Training Size")

  for (file in files)
  {
    i = i + 1
    data <- read.csv(paste(dir, file, sep=""))
    data$AVE <- c(rowMeans(data[,3:5]))
    lines(data$Size, data$AVE, type="o", lwd=4, col=colors[i])
  }
  abline(h=0.5, lwd=4, lty=2, col="red")

  i = 5
  files = c("scenario5_fastmnn_umap_size.csv", "scenario5_cca_umap_size.csv", "scenario5_sctransform_umap_size.csv")
  dir = "../output/"
  data <- read.csv(paste(dir, "scenario5_harmony_umap_size.csv", sep=""))
  data$AVE <- c(rowMeans(data[,3:5]))
  plot(data$Size, data$AVE, ylim=c(0.4,1), type="o", lwd=4, col=colors[i], ylab="Classification Accuracy", xlab="Training Size")

  for (file in files)
  {
    i = i + 1
    data <- read.csv(paste(dir, file, sep=""))
    data$AVE <- c(rowMeans(data[,3:5]))
    lines(data$Size, data$AVE, type="o", lwd=4, col=colors[i])
  }
  abline(h=0.5, lwd=4, lty=2, col="red")
  dev.off()
}

figure7 <- function(pdf_file)
{
  colors = brewer.pal(n = 8, name = "Set2")
  files = c("panc8_bias.csv", "python_umap_bias.csv")
  pdf(paste("../output/figures/", pdf_file, sep=""), width=12, height=6)
  par(mfrow=c(1,2))
  dir = "../output/"
  d1 <- read.csv(paste(dir, files[1], sep=""), row.names=1)
  d2 <- read.csv(paste(dir, files[2], sep=""), row.names=1)
  data <- rbind(d1, d2[, colnames(d2) %in% colnames(d1)])
  boxplot(t(data), ylab="Detection Accuracy", xlab="Integration Method")
  abline(h=1/length(colnames(data)), col="red")
  colors = c("#1B9E77", "#1B9E77", "#1B9E77", "#1B9E77", "#1B9E77", "#1B9E77", "#1B9E77", "#1B9E77", "#7570B3", "#7570B3", "#7570B3", "#7570B3", "#7570B3", "#7570B3", "#7570B3", "#7570B3", "#E7298A", "#E7298A", "#E7298A", "#E7298A", "#E7298A", "#E7298A", "#E7298A", "#E7298A")


  boxplot(data, col=colors, ylab="Detection Accuracy", xlab="Dataset", names=c(1:8, 1:8, 1:8))
  abline(h=1/length(colnames(data)), col="red")
  dev.off()
}

figure8 <- function(pdf_file)
{
  colors = brewer.pal(n = 8, name = "Set2")
  files = c("harmony_umap.csv", "fastmnn_umap.csv", "cca_umap.csv", "sctransform_umap.csv")
  labels = c("harmony_labels.csv", "fastmnn_labels.csv", "cca_labels.csv", "sctransform_labels.csv")
  pdf(paste("../output/figures/", pdf_file, sep=""), width=12, height=9)
  par(mfrow=c(3,4))
  dir = "../data/scenario_13/"
  for (i in c(1:length(files)))
  {
    data = read.csv(paste(dir, files[i], sep=""))
    label = read.csv(paste(dir, labels[i], sep=""))
    plot(data[,2:3])
    points(data[which(label$CELL_LINE=="HCC1395"), 2:3], col=colors[1])
    points(data[which(label$CELL_LINE=="HCC1395BL"), 2:3], col=colors[2])
  }

  dir = "../data/scenario_14/"
  for (i in c(1:length(files)))
  {
    data = read.csv(paste(dir, files[i], sep=""))
    label = read.csv(paste(dir, labels[i], sep=""))
    plot(data[,2:3])
    points(data[which(label$CELL_LINE=="HCC1395"), 2:3], col=colors[1])
    points(data[which(label$CELL_LINE=="HCC1395BL"), 2:3], col=colors[2])
  }

  dir = "../data/scenario_15/"
  for (i in c(1:length(files)))
  {
    data = read.csv(paste(dir, files[i], sep=""))
    label = read.csv(paste(dir, labels[i], sep=""))
    plot(data[,2:3])
    points(data[which(label$CELL_LINE=="HCC1395"), 2:3], col=colors[1])
    points(data[which(label$CELL_LINE=="HCC1395BL"), 2:3], col=colors[2])
  }

  dev.off()
}

figure9 <- function(pdf_file)
{
  colors = brewer.pal(n = 8, name = "Set2")
  files = c("harmony_umap.csv", "fastmnn_umap.csv", "cca_umap.csv", "sctransform_umap.csv")
  labels = c("harmony_labels.csv", "fastmnn_labels.csv", "cca_labels.csv", "sctransform_labels.csv")
  pdf(paste("../output/figures/", pdf_file, sep=""), width=12, height=9)
  par(mfrow=c(3,4))
  dir = "../data/scenario_46/"
  for (i in c(1:length(files)))
  {
    data = read.csv(paste(dir, files[i], sep=""))
    label = read.csv(paste(dir, labels[i], sep=""))
    plot(data[,2:3])
    points(data[which(label$CELL_LINE=="HCC1395"), 2:3], col=colors[1])
    points(data[which(label$CELL_LINE=="HCC1395BL"), 2:3], col=colors[2])
  }

  dir = "../data/scenario_47/"
  for (i in c(1:length(files)))
  {
    data = read.csv(paste(dir, files[i], sep=""))
    label = read.csv(paste(dir, labels[i], sep=""))
    plot(data[,2:3])
    points(data[which(label$CELL_LINE=="HCC1395"), 2:3], col=colors[1])
    points(data[which(label$CELL_LINE=="HCC1395BL"), 2:3], col=colors[2])
  }

  dir = "../data/scenario_48/"
  for (i in c(1:length(files)))
  {
    data = read.csv(paste(dir, files[i], sep=""))
    label = read.csv(paste(dir, labels[i], sep=""))
    plot(data[,2:3])
    points(data[which(label$CELL_LINE=="HCC1395"), 2:3], col=colors[1])
    points(data[which(label$CELL_LINE=="HCC1395BL"), 2:3], col=colors[2])
  }

  dev.off()
}

figure10 <- function(pdf_file)
{
  files = c("harmony_harmonized_28_48_cell_line.csv", "fastmnn_harmonized_28_48_cell_line.csv",
            "cca_pca_28_48_cell_line.csv", "sctransform_pca_28_48_cell_line.csv")
  grDevices::pdf(paste("../output/figures/", pdf_file, sep=""), width=12, height=6)
  utils::par(mfrow=c(2,4))
  dir = "../output/"
  labels = c("HCC1395", "HCC1395BL", "HCC1395", "HCC1395BL", "HCC1395", "HCC1395BL")
  colors = c("#1B9E77", "#1B9E77", "#7570B3", "#7570B3", "#E7298A", "#E7298A")
  for (file in files)
  {
    data <- read.csv(paste(dir, file, sep=""), row.names=1)
    boxplot(data, ylim=c(0, 1), col=colors, names=labels,  ylab="Classification Accuracy", xlab = "Paired Duplicated Datasets")
    abline(h=0.5, col="red", lty=2)
  }

  colors = brewer.pal(n = 8, name = "Set2")
  names = c(1: length(rownames(data)))
  i = 5
  for (file in files)
  {
    data <- read.csv(paste(dir, file, sep=""), row.names=1)
    boxplot(t(data), ylim=c(0, 1), names=names, col=colors[i], ylab="Classification Accuracy", xlab="Paired Datasets")
    abline(h=0.5, col="red", lty=2)
    i = i + 1
  }

  dev.off()
}

figure11 <- function(pdf_file)
{
  files = c("harmony_harmonized_28_48.csv", "fastmnn_harmonized_28_48.csv",
            "cca_pca_28_48.csv", "sctransform_pca_28_48.csv")
  pdf(paste("../output/figures/", pdf_file, sep=""), width=9, height=12)
  par(mfrow=c(4,3))
  dir = "../output/"
  colors = c("#1B9E77", "#1B9E77", "#7570B3", "#7570B3", "#E7298A", "#E7298A")
  for (file in files)
  {
    data <- read.csv(paste(dir, file, sep=""), row.names=1)
    labels <- rownames(data)
    boxplot(t(data[,1:4]), ylim=c(0, 1), col=colors[1],  ylab="Classification Accuracy", xlab = "Paired Duplicated Datasets")
    abline(h=0.25, col="red", lty=2)
    boxplot(t(data[,5:8]), ylim=c(0, 1), col=colors[2],  ylab="Classification Accuracy", xlab = "Paired Duplicated Datasets")
    abline(h=0.25, col="red", lty=2)
    boxplot(t(data[,9:12]), ylim=c(0, 1), col=colors[3],  ylab="Classification Accuracy", xlab = "Paired Duplicated Datasets")
    abline(h=0.25, col="red", lty=2)
  }

  dev.off()
}
figure12 <- function(pdf_file)
{
  files = c("harmony_harmonized_49_69.csv", "fastmnn_harmonized_49_69.csv",
            "cca_pca_49_69.csv", "sctransform_pca_49_69.csv")
  pdf(paste("../output/figures/", pdf_file, sep=""), width=12, height=12)
  par(mfrow=c(4,4))
  dir = "../output/"
  colors = c("#1B9E77", "#1B9E77", "#1B9E77", "#1B9E77",
             "#7570B3", "#7570B3", "#7570B3", "#7570B3",
             "#E7298A", "#E7298A", "#E7298A", "#E7298A")
  for (file in files)
  {
    data <- read.csv(paste(dir, file, sep=""))
    data$X <- NULL
    labels <- c("1", "2", "3", "4", "1", "2", "3", "4", "1", "2", "3", "4")
    boxplot(data, ylim=c(0, 1), col=colors,  names=labels, ylab="Classification Accuracy", xlab = "Duplicated Datasets")
    legend("topright", c("FFNN", "XGB", "kNN"), fill = c("#1B9E77", "#7570B3", "#E7298A"), bty='n')
    abline(h=0.25, col="red", lty=2)
  }

  colors = brewer.pal(n = 8, name = "Set2")[5:8]
  i = 1
  for (file in files)
  {
    data <- read.csv(paste(dir, file, sep=""))
    data$X <- NULL
    labels <- rownames(data)
    df <- data.frame(data[,1], data[,4])
    barplot(as.matrix(t(df)), beside=T, ylim=c(0, 0.4), ylab="Classification Accuracy", xlab = "Duplicated Datasets")
    abline(h=0.25, col="red", lty=2)
    i = i + 1
  }

  i = 1
  for (file in files)
  {
    data <- read.csv(paste(dir, file, sep=""))
    data$X <- NULL
    labels <- rownames(data)
    df <- data.frame(data[,5], data[,8])
    barplot(as.matrix(t(df)), beside=T, ylim=c(0, 0.4), ylab="Classification Accuracy", xlab = "Duplicated Datasets")
    abline(h=0.25, col="red", lty=2)
    i = i + 1
  }

  i = 1
  for (file in files)
  {
    data <- read.csv(paste(dir, file, sep=""))
    data$X <- NULL
    labels <- rownames(data)
    df <- data.frame(data[,9], data[,12])
    barplot(as.matrix(t(df)), ylim=c(0, 1), beside=T,  ylab="Classification Accuracy", xlab = "Duplicated Datasets")
    abline(h=0.25, col="red", lty=2)
    i = i + 1
  }
  dev.off()

}

figure14 <- function(pdf_file)
{
  files = c("harmony_harmonized_49_69.csv", "fastmnn_harmonized_49_69.csv",
            "cca_pca_49_69.csv", "sctransform_pca_49_69.csv")
  pdf(paste("../output/figures/", pdf_file, sep=""), width=12, height=9)
  par(mfrow=c(3,4))
  dir = "../output/"

  colors = brewer.pal(n = 8, name = "Set2")[5:8]
  i = 1
  for (file in files)
  {
    data <- read.csv(paste(dir, file, sep=""))
    data$X <- NULL
    labels <- rownames(data)
    barplot(as.matrix(t(rowMeans(data[,1:4]))), ylim=c(0, 0.4), col=colors[i],  ylab="Classification Accuracy", xlab = "Duplicated Datasets")
    abline(h=0.25, col="red", lty=2)
    i = i + 1
  }

  i = 1
  for (file in files)
  {
    data <- read.csv(paste(dir, file, sep=""))
    data$X <- NULL
    labels <- rownames(data)
    barplot(as.matrix(t(rowMeans(data[,5:8]))), ylim=c(0, 0.4), col=colors[i],  ylab="Classification Accuracy", xlab = "Duplicated Datasets")
    abline(h=0.25, col="red", lty=2)
    i = i + 1
  }

  i = 1
  for (file in files)
  {
    data <- read.csv(paste(dir, file, sep=""))
    data$X <- NULL
    labels <- rownames(data)
    barplot(as.matrix(t(rowMeans(data[,9:12]))), ylim=c(0, 0.4), col=colors[i],  ylab="Classification Accuracy", xlab = "Duplicated Datasets")
    abline(h=0.25, col="red", lty=2)
    i = i + 1
  }
  dev.off()
}

figure13 <- function(pdf_file)
{
  files = c("harmony_harmonized_49_69.csv", "fastmnn_harmonized_49_69.csv",
            "cca_pca_49_69.csv", "sctransform_pca_49_69.csv")
  pdf(paste("../output/figures/", pdf_file, sep=""), width=12, height=12)
  par(mfrow=c(4,2))
  dir = "../output/"

  colors = brewer.pal(n = 8, name = "Set2")[5:8]
  for (i in seq(1:4))
  {
    d1 <- read.csv(paste(dir, "harmony_harmonized_49_69.csv", sep=""))
    d1$X <- NULL

    d2 <- read.csv(paste(dir, "fastmnn_harmonized_49_69.csv", sep=""))
    d2$X <- NULL

    d3 <- read.csv(paste(dir, "cca_pca_49_69.csv", sep=""))
    d3$X <- NULL

    d4 <- read.csv(paste(dir, "sctransform_pca_49_69.csv", sep=""))
    d4$X <- NULL

    df1 <- data.frame(d1[,i], d2[,i], d3[,i], d4[,i])
    barplot(as.matrix(t(df1)), col=colors, ylim=c(0, 0.4), beside=T)
    abline(h=0.25, col="red", lty=2)
    df2 <- data.frame(d1[,i+4], d2[,i+4], d3[,i+4], d4[,i+4])
    barplot(as.matrix(t(df2)), col=colors, ylim=c(0, 0.4), beside=T)
    abline(h=0.25, col="red", lty=2)
  }
  dev.off()
}
#figure14("figure14.pdf")
#figure13("figure13.pdf")
#figure12("figure12.pdf")
#figure11("figure11.pdf")
#figure9("figure9.pdf")
#figure10("figure10.pdf")
#figure8("figure8.pdf")
#figure7("figure7.pdf")
#figure6("figure6.pdf")
#figure5("figure5.pdf")
#figure4("47_latex_table.csv", "figure4_1.pdf")
#figure4("48_latex_table.csv", "figure4_2.pdf")
#figure1()
#figure2("benchmark47_", "figure2_1.pdf")
#figure2("benchmark48_", "figure2_2.pdf")
#figure2("benchmark46_", "figure2_3.pdf")
#figure3()
