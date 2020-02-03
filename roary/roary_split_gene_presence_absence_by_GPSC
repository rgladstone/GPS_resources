library(data.table)

GPSCs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,30,31,32,33,34,37,38,39,40,41,43,47,48,50,51,52,53,54,55,56,57,58,61,62,65,67,68,70,72,76,77,78,79,80,81,90,91,93,94,97,103,105,108,115,117,131)

GPSC_lanes <- fread("GPSC_lanes.csv",data.table = FALSE)
roary <-  fread("gene_presence_absence_minimised.csv",data.table = FALSE)

emptycols <- colnames(roary)[5:10] 
for(i in emptycols){
  roary[,i] <- ""
}
std_cols <- colnames(roary)[1:14]

summaryall <- matrix(nrow = 0, ncol = 7)
colnames(summaryall) <- c("GPSC","count","genes","core>=99", "softcore>=95<99","shell>=15<95","cloud<15")

for (GPSC in GPSCs){
  lanes <- GPSC_lanes[GPSC_lanes$GPSC == GPSC,2]
  columnselect <- c(std_cols, lanes)
  subset_lanes_roary <- roary[,columnselect]
  subset_lanes_roary[is.na(subset_lanes_roary)] <-0
  subset_genes_roary <- subset_lanes_roary[rowSums(subset_lanes_roary[15:ncol(subset_lanes_roary)])!=0, ]
  subset_genes_roary <- subset_genes_roary[rev(order(rowSums(subset_genes_roary[15:ncol(subset_genes_roary)]))),]
  subset_genes_roary[,"No. isolates"] <- rowSums(subset_genes_roary[15:ncol(subset_genes_roary)])
  GPSC_count <- length(lanes)
  core99 <-  sum(subset_genes_roary$`No. isolates` >= GPSC_count*0.99)
  softcore95 <- sum(subset_genes_roary$`No. isolates` < GPSC_count*0.99 & subset_genes_roary$`No. isolates` >= GPSC_count*0.95)
  shell15 <- sum(subset_genes_roary$`No. isolates` < GPSC_count*0.95 & subset_genes_roary$`No. isolates` >= GPSC_count*0.15)
  cloud0 <- sum(subset_genes_roary$`No. isolates` < GPSC_count*0.15)
  total <- nrow(subset_genes_roary)
  summary <- c(GPSC,GPSC_count,total,core99,softcore95,shell15,cloud0)
  summaryall <- rbind(summaryall,summary)
  subset_genes_roary[subset_genes_roary==0] <- ""
  write.csv(file=paste("GPSC",GPSC,"_gene_presence_absence.csv", sep = ""),subset_genes_roary, row.names = FALSE)
}

write.csv(file = "summary_stats_by_GPSC.csv", summaryall, row.names = FALSE)
