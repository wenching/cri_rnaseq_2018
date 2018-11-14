# CRI RNAseq 2018

RNA-Seq Analysis Pipeline using the [CRI](http://cri.uchicago.edu/) HPC system

# CAUTION

**THIS PACKAGE IS LARGE, PLEASE DO NOT DOWNLOAD IT TO YOUR HOME DIRECTORY**  
**USE OTHER LOCATION LIKE /gpfs/data/bioinformatics/username**

## Dataset

The RNA-seq data used in this tutorial are from a [published paper](https://www.ncbi.nlm.nih.gov/pubmed/25499759) that explores PRDM11 and lymphomagenesis. We will use the data from the PRDM11 knockdown and wildtype samples. You are welcome to explore the full dataset on GEO ([GSE56065](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56065)).

>  
> Fog et al., 2015, Loss of PRDM11 promotes MYC-driven lymphomagenesis, Blood 125:1272-1281 The Cancer Genome Atlas Research Network, 2011, Integrated genomic analyses of ovarian carcinoma, Nature, 474:609–615  
>  

## File description

This repository contains the following items:
- ```docs/README.html``` - the main tutorial documentation
- ```example/``` - folder accommodating all data (i.e., metadata file, configuration file, sequencing data folder, and references folder) for running this pipeline
- ```README.md``` - this description file
- ```SRC/``` - automatic pipelines for RNA-seq analysis

### Prerequisites

* [R](https://www.r-project.org/)
* [CRAN](https://cran.r-project.org/)
* [Bioconductor](https://www.bioconductor.org/)
* [Python 3](https://www.python.org/download/releases/3.0/)
* [BigDataScript](https://pcingola.github.io/BigDataScript/)

### Installing

```bash
# CAUTION

**THIS PACKAGE IS LARGE, PLEASE DO NOT DOWNLOAD IT TO YOUR HOME DIRECTORY**  
**USE OTHER LOCATION LIKE /gpfs/data/bioinformatics/username**

# download the package via 'git clone'
git clone git@github.com:wenching/cri_rnaseq_2018.git


# Or, download the latest package via 'wget'
wget https://github.com/wenching/cri_rnaseq_2018/archive/master.tar.gz .

# uncompress the tarball file, if using wget
# tar -zxvf master.tar.gz

# change folder name, if using wget
# mv cri_rnaseq_2018-master cri_rnaseq_2018


# change working directory to the package directory
cd cri_rnaseq_2018

# load modules
$ module purge;module load gcc udunits python/3.6.0 R/3.5.0; module update

# This step is optional but it will install all necessary R packages ahead.
# In case the pipeline was terminated due to the failure of R package installation later when running the pipeline.
$ Rscript --vanilla SRC/R/util/prerequisite.packages.R

# create directories and generate all necessary scripts
$ bash Build_RNAseq.DLBC.sh

# run the entire pipeline with just this command
$ bash Submit_cri_rnaseq_2018.sh
```

## Contributing

Please read [CONTRIBUTING.md](https://github.com/wenching/cri_rnaseq_2018/blob/master/CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

For the versions available, see the [tags on this repository](https://github.com/wenching/cri_rnaseq_2018/tags).

## Authors

* **Wen-Ching Chan** - *Initial work* - [Wen-Ching](https://github.com/wenching)

See also the list of [contributors](https://github.com/wenching/cri_rnaseq_2018/graphs/contributors) who participated in this project.

## License

This project is licensed under the [LGPLv3](https://www.gnu.org/licenses/lgpl-3.0.en.html) License - see the [LICENSE](LICENSE) file with a copy for details

## Acknowledgments

* Thank [Kyle Hernandez](https://github.com/kmhernan) for providing a private Git repository in the begining of this project
* Thank [Riyue Bao](https://github.com/riyuebao) for providing tutorial dataset, figures, and the pipeline in the Perl version
* Thank [PurpleBooth](https://gist.github.com/PurpleBooth) for providing the templates of [README.md](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2) & [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426)


