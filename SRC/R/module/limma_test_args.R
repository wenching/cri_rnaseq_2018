rm(list = ls())


if(Sys.info()['sysname'] == "Linux") {
  (args <- structure(
    list(
      b.filter.low.expr = "TRUE",
      c.appl = "RNAseq",
      c.criteria = "fc_q",
      c.gtf.file.path = "/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.gtf",
      c.in.dir.path = "/group/bioinformatics/CRI_RNAseq_2018/example/DLBC_full/RNAseq/Quantification/featurecounts/star",
      c.log.file.path = "/group/bioinformatics/CRI_RNAseq_2018/example/DLBC_full/RNAseq/DEG/limma/featurecounts/star/run.call.limma.featurecounts.star.DLBC.log",
      c.meta.data.path = "/group/bioinformatics/CRI_RNAseq_2018/example/DLBC_full.metadata.txt",
      c.out.file.path = "/group/bioinformatics/CRI_RNAseq_2018/example/DLBC_full/RNAseq/DEG/limma/featurecounts/star/DLBC.star.featurecounts.limma.count.txt",
      c.v.case.2.ctrl = "KO_vs_WT",
      d.fold.change = 1.5,
      d.pval = 0.1
    ),
    .Names = c(
      "b.filter.low.expr",
      "c.appl",
      "c.criteria",
      "c.gtf.file.path",
      "c.in.dir.path",
      "c.log.file.path",
      "c.meta.data.path",
      "c.out.file.path",
      "c.v.case.2.ctrl",
      "d.fold.change",
      "d.pval"
    )
  ))
} else if(Sys.info()['sysname'] == "Darwin") {
  (args <- structure(
    list(
      b.filter.low.expr = "TRUE",
      c.appl = "RNAseq",
      c.criteria = "fc_q",
      c.gtf.file.path = "/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.gtf",
      c.in.dir.path = "/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC_full/RNAseq/Quantification/featurecounts/star",
      c.log.file.path = "/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC_full/RNAseq/DEG/limma/featurecounts/star/run.call.limma.featurecounts.star.DLBC.log",
      c.meta.data.path = "/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC_full.metadata.txt",
      c.out.file.path = "/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC_full/RNAseq/DEG/limma/featurecounts/star/DLBC.star.featurecounts.limma.count.txt",
      c.v.case.2.ctrl = "KO_vs_WT",
      d.fold.change = 1.5,
      d.pval = 0.1
    ),
    .Names = c(
      "b.filter.low.expr",
      "c.appl",
      "c.criteria",
      "c.gtf.file.path",
      "c.in.dir.path",
      "c.log.file.path",
      "c.meta.data.path",
      "c.out.file.path",
      "c.v.case.2.ctrl",
      "d.fold.change",
      "d.pval"
    )
  ))
}

