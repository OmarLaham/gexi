#https://rdrr.io/cran/RNOmni/man/RankNorm.html
library('RNOmni')


dds <- readRDS("002_dds.rds")

#pre-normalized
countsData <- as.matrix(dds@assays@data$counts)

#init iNTCounts with same values as countsData, but then iterate to apply iNT
iNTCounts <- countsData

for (i in 1:nrow(countsData)) {
	# Rank normalize
	iNTCounts[i,] <- RankNorm(countsData[i,])
}

#export tsv
write.table(iNTCounts, file='inverse_normal_transform_counts.tsv', sep='\t', row.names = T, col.names = T)

print("Exported! Done!")
