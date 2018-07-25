"""GSEA"""
import os
import sys
import logging
import datetime

import lib.clusterprofiler

import util.ddictfunc

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def gsea(args, sw_cfg, task_cfg):
    """
    GSEA
    
    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {clusterprofiler}
    :parm task_cfg: a dictionary of corresponding task configureation {main, references, clusterprofiler}
    :returns: returns corresponding code snippets, in which will be written to shell_script_path if provided
    :raises keyError: NA
    """
    # logging
    if args.log_file is None:
        log_file_path = '{}.{}.log'.format(
            SELF_FILE_PATH,
            datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
    else:
        log_file_path = args.log_file

    formatter = "%(asctime)-15s %(levelname)-8s %(message)s"
    logging.basicConfig(
        level=[logging.NOTSET, logging.DEBUG, logging.INFO,
               logging.WARNING, logging.ERROR, logging.CRITICAL][1],
        format=formatter,
        filename=log_file_path)

    for gseatool in task_cfg.keys():
        if gseatool in ["main", "references"]: continue
        logging.info("[ {} ] GSEA tool: {}".format(SELF_FILE, gseatool))
        
        for quantifier in task_cfg[gseatool].keys():
            if quantifier in ["main"]: continue
            logging.info("[ {} ] Quantifier: {}".format(SELF_FILE, quantifier))
            
            for aligner in task_cfg[gseatool][quantifier].keys():
                if aligner in ["main"]: continue
                logging.info("[ {} ] Aligner: {}\n".format(SELF_FILE, aligner))
    
                for proj_comp in task_cfg[gseatool][quantifier][aligner].keys():
                    if proj_comp in ["main"]: continue
                    logging.info("[ {} ] Comparison: {}".format(SELF_FILE, proj_comp))                                
    
                    out_dir_t_path = task_cfg[gseatool][quantifier][aligner][proj_comp]['out_dir_path']
                    if not os.path.exists(out_dir_t_path):
                        logging.debug("MKDIR: [ $out_dir_t_path ]")
                        os.makedirs(out_dir_t_path, exist_ok=True)
        
                    task_cfg_proj = {
                        'references': task_cfg['references'],
                        gseatool: task_cfg[gseatool][quantifier][aligner][proj_comp]
                    }
                    
                    task_cfg_proj[gseatool]['in_file_path_dict'] = {}
                    
                    if True:
                        task_cfg_proj[gseatool]['in_file_path_dict'].update(
                            {
                                proj_comp: task_cfg['main'][quantifier][aligner][proj_comp]['in_file_path_list'][0]
                            }
                        )

                    #util.ddictfunc.pprint_ddicts(task_cfg_proj)
                    '''
                    {'clusterprofiler': {'in_file_path_dict': {'DLBC': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/LociStat/featurecounts/star/DLBC/DLBC.star.featurecounts.overlap.txt'},
                                         'log_file_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/GSEA/featurecounts/star/DLBC/run.call.featurecounts.star.DLBC.log',
                                         'out_dir_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/GSEA/featurecounts/star/DLBC',
                                         'out_file_path_list': ['/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/GSEA/featurecounts/star/DLBC/DLBC.star.featurecounts.enrichGO.ALL.txt'],
                                         'shell_script_path': '/Users/wenching/Desktop/Sync/CRI/CRI-Pipeline/CRI_RNAseq_2018/example/DLBC/RNAseq/shell_scripts/run.lociStat.featurecounts.star.DLBC.sh'},
                     'references': {'grch38': {'anno_bed': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.bed12',
                                               'anno_bed12': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.bed12',
                                               'anno_gtf': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.gtf',
                                               'anno_refflat': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.refFlat.txt',
                                               'anno_ribosome_rna_bed': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.rRNA.bed',
                                               'anno_ribosome_rna_interval': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/gencode.v24.primary_assembly.annotation.chr11.rRNA.interval_list',
                                               'chrom_size': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/GRCh38.primary_assembly.genome.chr11.chrom.size',
                                               'chrs': 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM',
                                               'genome': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/GRCh38.primary_assembly.genome.chr11.fa',
                                               'genomedict': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11/GRCh38.primary_assembly.genome.chr11.dict',
                                               'star_index': '/group/bioinformatics/CRI_RNAseq_2018/example/reference/GRCh38.primary_Gencode24_50bp_chr11'}}}
                    
                    '''

                    if(task_cfg['main']['application'] == "RNAseq"):
                        sw_cfg_tool = util.ddictfunc.subset(sw_cfg, [gseatool])
                        lib.clusterprofiler.clusterprofiler(args, sw_cfg_tool, task_cfg_proj)
    
                
                    logging.debug("[ {} ] Comparison: {} - DONE\n".format(SELF_FILE, proj_comp))                
                logging.debug("[ {} ] Aligner: {} - DONE\n".format(SELF_FILE, aligner))
            logging.debug("[ {} ] Quantifier: {} - DONE\n".format(SELF_FILE, quantifier))

    return
