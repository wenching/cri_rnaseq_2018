rm(list = ls())


if(Sys.info()['sysname'] == "Linux") {
  (args <- structure(
    list(
    ),
    .Names = c(
    )
  ))
} else if(Sys.info()['sysname'] == "Darwin") {
  (args <- structure(
    list(
      c.in.file.path = "/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC_full/RNAseq/DEG/deseq2/featurecounts/star/DLBC.star.featurecounts.deseq2.count.txt",
      c.out.file.path = "/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/QuantQC/featurecounts/star/DLBC.star.featurecounts.pca.pdf",
      c.v.id = "KO01,KO02,KO03,WT01,WT02,WT03",
      c.v.group = "KO,KO,KO,WT,WT,WT",
      c.log.file.path = "/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/QuantQC/featurecounts/star/run.quantQC.pca.featurecounts.star.DLBC.log"
    ),
    .Names = c(
      "c.in.file.path",
      "c.out.file.path",
      "c.v.id",
      "c.v.group",
      "c.log.file.path"
    )
  ))
}

