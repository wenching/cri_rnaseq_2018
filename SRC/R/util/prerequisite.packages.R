cat("# INSTALL NECESSARY PACKAGES - START")


if(! "futile.logger" %in% rownames(installed.packages())) {
  if(! "devtools" %in% rownames(installed.packages())) {
    install.packages(
      "devtools",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  }
  
  library("devtools")
  
  devtools::install_github(
    repo = "zatonovo/futile.logger"
  )
}


if(! "DESeq2" %in% rownames(installed.packages())) {
  if(! "BiocManager" %in% rownames(installed.packages()))
    install.packages(
      "BiocManager",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  BiocManager::install("DESeq2", version = "devel")
}



if(! "vsn" %in% rownames(installed.packages())) {
  if(! "hexbin" %in% rownames(installed.packages())) {
    install.packages(
      "hexbin",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  }
  
  if(! "BiocManager" %in% rownames(installed.packages()))
    install.packages(
      "BiocManager",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  BiocManager::install("vsn", version = "devel")
}
library("vsn")


if(! "FactoMineR" %in% rownames(installed.packages())) {
  install.packages(
    "FactoMineR",
    dependencies = TRUE,
    repos = "https://cloud.r-project.org"
  )
}


if(! "factoextra" %in% rownames(installed.packages())) {
  install.packages(
    "factoextra",
    dependencies = TRUE,
    repos = "https://cloud.r-project.org"
  )
}


if(! "Biobase" %in% rownames(installed.packages())) {
  if(! "BiocManager" %in% rownames(installed.packages()))
    install.packages(
      "BiocManager",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  BiocManager::install("Biobase", version = "devel")
}


if(! "Vennerable" %in% rownames(installed.packages())) {
  if(! "RBGL" %in% rownames(installed.packages())) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("RBGL")
  }
  
  if(! "graph" %in% rownames(installed.packages())) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("graph")
  }
  
  if(! "reshape" %in% rownames(installed.packages())) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("reshape")
  }
  
  if(! "devtools" %in% rownames(installed.packages())) {
    install.packages(
      "devtools",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  }
  library("devtools")
  
  install_github(repo = 'js229/Vennerable')
  #browseURL("http://github.com/js229/Vennerable")
}


if(! "VennDiagram" %in% rownames(installed.packages())) {
  install.packages(
    "VennDiagram",
    dependencies = TRUE,
    repos = "https://cloud.r-project.org"
  )
}


c.org.db <- 'org.Hs.eg.db'
if(! c.org.db %in% rownames(installed.packages())) {
  if(! "BiocManager" %in% rownames(installed.packages()))
    install.packages(
      "BiocManager",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  BiocManager::install(c.org.db, version = "devel")
}


c.txdb <- 'TxDb.Hsapiens.UCSC.hg38.knownGene'
if(! c.txdb %in% rownames(installed.packages())) {
  if(! "BiocManager" %in% rownames(installed.packages()))
    install.packages(
      "BiocManager",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  BiocManager::install(c.txdb, version = "devel")
}


if(! "TxDb.Hsapiens.UCSC.hg19.knownGene" %in% rownames(installed.packages())) {
  if(! "BiocUpgrade" %in% rownames(installed.packages())) {
    if(! "BiocManager" %in% rownames(installed.packages()))
      install.packages(
        "BiocManager",
        dependencies = TRUE,
        repos = "https://cloud.r-project.org"
      )
    BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", version = "devel")
  }
  
  if(! "devtools" %in% rownames(installed.packages())) {
    install.packages(
      "devtools",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  }
  library(devtools)
  
  install_github(
    "GuangchuangYu/ChIPseeker",
    build_vignettes = FALSE,
    repos = BiocInstaller::biocinstallRepos(),
    dependencies = TRUE
  )
}


if(! "KEGGREST" %in% rownames(installed.packages())) {
  if(! "BiocManager" %in% rownames(installed.packages()))
    install.packages(
      "BiocManager",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  BiocManager::install("KEGGREST", version = "devel")
}


if(! "edgeR" %in% rownames(installed.packages())) {
  if(! "BiocManager" %in% rownames(installed.packages()))
    install.packages(
      "BiocManager",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  BiocManager::install("edgeR", version = "devel")
}


if(! "limma" %in% rownames(installed.packages())) {
  if(! "BiocManager" %in% rownames(installed.packages()))
    install.packages(
      "BiocManager",
      dependencies = TRUE,
      repos = "https://cloud.r-project.org"
    )
  BiocManager::install("limma", version = "devel")
}


if(! "pheatmap" %in% rownames(installed.packages())) {
  install.packages(
    "pheatmap",
    dependencies = TRUE,
    repos = "https://cloud.r-project.org"
  )
}


cat("# INSTALL NECESSARY PACKAGES - DONE")



print(sessionInfo())
q(save = "no")
