"""Picard"""
import os
import sys
import logging
import datetime

import util.ddictfunc

SELF_FILE_PATH = os.path.realpath(__file__)
SELF_DIR_PATH = os.path.dirname(SELF_FILE_PATH)
SELF_FILE = os.path.basename(SELF_FILE_PATH)


def collect_rna_seq_metrics(args, sw_cfg, task_cfg):
    """
    Picard::CollectRnaSeqMetrics

    :parm args: an argparse.Namespace object of main argumens {.log_file}
    :parm sw_cfg: a dictionary of corresponding software configureation {picard}
    :parm task_cfg: a dictionary of corresponding task configureation {module_name, references}
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

    logging.info("[ {} ] Make Shell Script\n".format(SELF_FILE))

    sw_cfg_picard = sw_cfg['picard']

    task_cfg_module = task_cfg['CollectRnaSeqMetrics']
    log_file_path_cmd = '' if 'log_file_path' not in task_cfg_module else '''&>> {}'''.format(task_cfg_module['log_file_path'])

    task_cfg_reference = task_cfg['references'][next(iter(task_cfg['references']))]

    shell_script_path = '' if 'shell_script_path' not in task_cfg_module else task_cfg_module['shell_script_path']

    cmd_list = []

    if sw_cfg_picard['module'] is not None:
        cmd_list.append(sw_cfg_picard['module'])

    cmd_list.append('''java \\
-Xmx{}G \\
-jar {} \\
CollectRnaSeqMetrics \\
INPUT={} \\
OUTPUT={} \\
REF_FLAT={} \\
RIBOSOMAL_INTERVALS={} \\
STRAND={} \\
CHART_OUTPUT={}.pdf \\
METRIC_ACCUMULATION_LEVEL={} \\
TMP_DIR={} \\
{}
\n'''.format(
         sw_cfg_picard['mem'],
         sw_cfg_picard['exe'],
         task_cfg_module['in_file_path_list'][0],
         '''{}.{}'''.format(
             task_cfg_module['out_file_base'],
             'RNA_Metrics'
         ),
         task_cfg_reference['anno_refflat'],
         task_cfg_reference['anno_ribosome_rna_interval'],
         task_cfg_module['strand_specificity'],
         '''{}.{}'''.format(
             task_cfg_module['out_file_base'],
             'RNA_Metrics'
         ),
         "SAMPLE",
         task_cfg_module['tmp_dir_path'],
         log_file_path_cmd
     ))

    if shell_script_path:
        try:
            with open(shell_script_path, "w") as outfile:
                outfile.write(
                    '''\n{}\n'''.format(
                        '\n\n'.join(cmd_list)
                    )
                )
                
                outfile.write("\nexit 0\n")
                outfile.close()
        except IOError as exc:
            print(exc)

    logging.debug("[ {} ] Make Shell Script - DONE\n".format(SELF_FILE))

    return cmd_list
