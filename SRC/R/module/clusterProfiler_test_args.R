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
      c.v.in.file.path = "/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/LociStat/featurecounts/star/DLBC/DLBC.star.featurecounts.overlap.txt",
      c.v.comp.pair = "DLBC",
      c.out.file.path = "/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/GSEA/featurecounts/star/DLBC/DLBC.star.featurecounts.enrichGO.ALL.txt",
      c.assembly = "grch38",
      c.log.file.path = "/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/GSEA/featurecounts/star/DLBC/run.GSEA.featurecounts.star.DLBC.log"
    ),
    .Names = c(
      "c.v.in.file.path",
      "c.v.comp.pair",
      "c.out.file.path",
      "c.assembly",
      "c.log.file.path"
    )
  ))
}

