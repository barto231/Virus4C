import argparse
import os
import sys
from pathlib import Path

VERSION = '1.1.6'


def get_script_path(path=None):
    if path is None:
        cmd_args = sys.argv
        needle = "--file="
        match = [arg for arg in cmd_args if needle in arg]
        if len(match) > 0:
            return os.path.normpath(match[0].replace(needle, ""))
        else:
            ls_vars = list(sys._getframe(1).f_locals.keys())
            if "fileName" in ls_vars:
                return os.path.normpath(sys._getframe(1).f_locals["fileName"])
            else:
                return os.path.normpath(sys._getframe(1).f_locals["ofile"])
    else:
        return path


def parse_args():
    parser = argparse.ArgumentParser(description="Values stored in the configuration file (conf.yml) are used by default.")
    parser.add_argument("-v", "--vpFile", metavar="/path/to/vp_file", type=str, required=True, help="path to viewpoint file")
    parser.add_argument("-f", "--fqFolder", metavar="/path/to/FASTQ_folder/", type=str, required=True, help="path to folder containing the FASTQ files")
    parser.add_argument("-o", "--outFolder", metavar="/path/to/output_folder/", type=str, required=True, help="path to the output folder")
    parser.add_argument("-c", "--confFile", metavar="/path/to/conf.yml/", type=str, default="conf.yml", help="path to configuration file [default %(default)s]")
    parser.add_argument("-z", "--cores", metavar="number", type=int, default=None, help="Number of cores for parallelization.")
    parser.add_argument("-w", "--wig", action="store_true", default=False, help="create wig files for all samples")
    parser.add_argument("-x", "--bigwig", action="store_true", default=False, help="create big wig files for all samples")
    parser.add_argument("-p", "--plot", action="store_true", default=False, help="Create viewpoint coverage plot for all samples.")
    parser.add_argument("-g", "--genomePlot", action="store_true", default=False, help="Create genomeplot for all samples (only possible if analysis is all in vpFile).")
    parser.add_argument("-t", "--tsv", action="store_true", default=False, help="Create tab separated values file for all samples")
    parser.add_argument("-b", "--bins", action="store_true", default=False, help="Corunt reads for binned regions.")
    return parser.parse_args()


def main():
    argsL = parse_args()

    if not os.path.exists(argsL.vpFile):
        raise ValueError("vpFile does not exist: {}".format(argsL.vpFile))
    if not os.path.exists(argsL.fqFolder):
        raise ValueError("fqFolder does not exist: {}".format(argsL.fqFolder))
    if not os.path.exists(argsL.outFolder):
        raise ValueError("outFolder does not exist: {}".format(argsL.outFolder))

    sys.path.append(os.path.dirname(get_script_path()))

    from functions import *
    import ShortRead
    import GenomicRanges
    import GenomicAlignments
    import caTools
    import config

    argsL.confFile = os.path.normpath(argsL.confFile)

    if argsL.confFile == 'conf.yml':
        argsL.confFile = os.path.normpath(get_script_path().replace('pipe4C.R',

argsL.confFile = re.sub('pipe4C\\.R', 'conf.yml', script_path) if argsL.confFile == 'conf.yml' else argsL.confFile

configOpt = createConfig(confFile=argsL.confFile)
configOpt.pipeline.version = 'VERSION'

if argsL.cores:
    configOpt.cores = int(argsL.cores)
if argsL.wSize:
    configOpt.wSize = int(argsL.wSize)
if argsL.nTop:
    configOpt.nTop = int(argsL.nTop)
if argsL.wig:
    configOpt.wig = True
if argsL.bigwig:
    configOpt.bigwig = True
if argsL.plot:
    configOpt.cisplot = argsL.plot
if argsL.genomePlot:
    configOpt.genomePlot = argsL.genomePlot

from pipe4py import Run_4Cpipeline

Run_4Cpipeline(
    VPinfo_file=argsL['vpFile'],
    FASTQ_F=argsL['fqFolder'],
    OUTPUT_F=argsL['outFolder'],
    configuration=configOpt
)