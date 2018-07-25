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
      c.v.in.file.path = "/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC_full/RNAseq/DEG/edger/featurecounts/star/DLBC.star.featurecounts.edger.test.DEG.txt,/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC_full/RNAseq/DEG/deseq2/featurecounts/star/DLBC.star.featurecounts.deseq2.test.DEG.txt,/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC_full/RNAseq/DEG/limma/featurecounts/star/DLBC.star.featurecounts.limma.test.DEG.txt",
      c.out.file.path = "/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/LociStat/featurecounts/star/DLBC/DLBC.star.featurecounts.overlap.txt",
      c.anchor.file.path = "/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC_full/RNAseq/DEG/deseq2/featurecounts/star/DLBC.star.featurecounts.deseq2.test.txt",
      c.log.file.path = "/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/LociStat/featurecounts/star/DLBC/run.call.featurecounts.star.DLBC.log"
    ),
    .Names = c(
      "c.v.in.file.path",
      "c.out.file.path",
      "c.anchor.file.path",
      "c.log.file.path"
    )
  ))
}

