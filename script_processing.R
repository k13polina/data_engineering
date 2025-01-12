library(DESeq2)
library(PCAtools)
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)



#### Анализ экспрессии: контроль качества ####


# Загрузка данных
expression.data <- read.table("C:/Users/poly1/OneDrive/Рабочий стол/data_en/ctrl_vs_exp.txt", header=TRUE, row.names=1, sep="\t", na.strings="NA", dec=".")
head(expression.data) 


# Корреляция биологических повторностей
correlation <- cor(expression.data)
correlation

cc <- round(correlation**2, 3) #округление: второй аргумент определяет количество символов после запятой
cc
write.table(cc, "C:/Users/poly1/OneDrive/Рабочий стол/data_en/leaf_corr.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE, na="NA")

plot(expression.data$ctrl1, expression.data$ctrl2)
plot(expression.data$ctrl1, expression.data$exp1)


# Иерархическая кластеризация
distance <- as.dist(1-cor(expression.data)) #чтобы посчитать насколько образцы не похожи друг на друга  
tree <- hclust(distance)
plot(tree, main = "Sample Dendrogram", font.main = 2, sub = '', ylab = "1-corr")



#### Поиск дифференциально экспрессирующихся генов с помощью DESeq2 ####
head(expression.data)
expr.matrix <- as.matrix(expression.data)

is.matrix(expression.data)
is.matrix(expr.matrix)


# Дизайн эксперимента
expr.design <- data.frame(row.names = colnames(expression.data), condition = c("control", "control", "experiment", "experiment"))
expr.design


# Анализ дифференциальной экспрессии
dds <- DESeqDataSetFromMatrix(countData = expr.matrix, colData = expr.design, design = ~ condition)
dds <- DESeq(dds)


# Principal Component Analysis - PCA 
vst <- assay(vst(dds))
p.all <- pca(vst, metadata = colData(dds), removeVar = 0.1)
biplot(p.all, legendPosition = 'right')

# Результаты анализа
res <- results(dds)
head(res)  
write.table(as.data.frame(res), "C:/Users/poly1/OneDrive/Рабочий стол/data_en/DESeq2_Leaf.txt", sep="\t", col.names=TRUE, row.names=TRUE, na="NA", quote=F)


# Нормализация на размер библиотеки
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized.counts <- counts(dds, normalized=TRUE)
summary(expression.data)
summary(normalized.counts)

write.table(normalized.counts, file="C:/Users/poly1/OneDrive/Рабочий стол/data_en/Norm_Leaf.txt", sep="\t", quote=F, col.names=NA)


# Оценка дисперсии
plotDispEsts(dds)


# Результаты анализа, визуализация
plotMA(dds,ylim=c(-2,2),main="DESeq2")


# Отбор дифференциально экспрессирующихся (ДЭ) генов
sum(res$padj < 0.05, na.rm=TRUE)
resSig <- subset(res, padj < 0.05)
head(res)
head(resSig)

sum(resSig$log2FoldChange > 1, na.rm=TRUE)
sum(resSig$log2FoldChange < -1, na.rm=TRUE)

resSigUp <- subset(resSig, log2FoldChange > 1)
head(resSigUp)
write.table(as.data.frame(resSigUp), "C:/Users/poly1/OneDrive/Рабочий стол/data_en/DESeq2_Leaf_Up.xlsx", sep="\t", col.names=TRUE, row.names=TRUE, na="NA", quote=F)

resSigDown <- subset(resSig, log2FoldChange < -1)
head(resSigDown)
write.table(as.data.frame(resSigDown), "C:/Users/poly1/OneDrive/Рабочий стол/data_en/DESeq2_Leaf_Down.xlsx", sep="\t", col.names=TRUE, row.names=TRUE, na="NA", quote=F)



#### Визуализация результатов ####
# Heatmap - тепловая карта
names <- c(row.names(resSigUp), row.names(resSigDown))
de.genes <- subset(normalized.counts, rownames(normalized.counts) %in% names)
head(de.genes)

Heatmap(de.genes, show_row_names = FALSE)

log.de.genes <- log10(de.genes)
head(log.de.genes)
Heatmap(log.de.genes, show_row_names = FALSE) #ничего не получается, т.к. там были 0 и когда мы сделали log, появились неопределнные значения

de.genes <- de.genes + 0.001 #теперь вместо 0 у нас очень маленькие числа
log.de.genes <- log10(de.genes)
Heatmap(log.de.genes, show_row_names = FALSE)


# Шкалирование
scaled.row.de.genes <- t(scale(t(log.de.genes)))
Heatmap(scaled.row.de.genes, show_row_names = FALSE)


# Volcano plot
de.result <- as.data.frame(res)
head(de.result)
ggplot(de.result, aes(log2FoldChange, -log(padj,10))) +  
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))

de.result <- de.result %>% 
  mutate(Expression = case_when(log2FoldChange >= log(2) & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -log(2) & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged"))

head(de.result)

ggplot(de.result, aes(log2FoldChange, -log(padj,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  theme_bw(base_size = 12, base_family = "") +
  ylim(0,20)
