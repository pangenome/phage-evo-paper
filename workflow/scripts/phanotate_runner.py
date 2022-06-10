#!/usr/bin/env python3
"""
Phanotate runner.
Usage:
    phanotate_runner.py ( --input_file_list=PATH ) ( --output_dir=PATH )
                        ( --out_format=FORMAT ) ( --threads=INT )

Options:
    -h --help                      Show this screen.
    -i --input_file_list=PATH      Path to a file containing a list of input fasta files.
    -o --output_dir=PATH           Path to the output directory.
    -f --out_format=FORMAT         Format of the output files choices=['tabular','genbank','fasta'] [default: genbank]
    -t --threads=INT               Number of threads to use for running individual anotations [default: 1].
"""
import os
import sys
import subprocess
import logging
from multiprocessing import Pool
from functools import partial
from docopt import docopt
logger = logging.getLogger(__name__)


def check_extensions(input_file):
    extensions = ('.fasta', 'fasta.gz', '.fa', '.fa.gz', '.fna', '.fna.gz')
    for ext in extensions:
        if input_file.endswith(ext):
            return ext
    raise ValueError("Input file does not have a valid extension: {}\nValid extensions are: {}".format(input_file, extensions))

def run_phanotate(input_file, output_dir, out_format, input_file_extension='', output_file_extension='',):
    """
    Run phanotate on a single input file.
    """
    out_fname = os.path.join(output_dir, os.path.basename(input_file).replace(input_file_extension, output_file_extension))
    try:
        cmd = [
            'phanotate.py',
            '-o', out_fname,
            '-f', out_format,
            input_file ]
        logger.info("Running phanotate on {}".format(input_file))
        logger.info("Running command: {}".format(' '.join(cmd)))
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        logger.error("Error running phanotate on {}".format(input_file))
        sys.exit(1)

def main(*args, **kwargs):
    logging.basicConfig(
        level = logging.INFO,
        datefmt="%Y-%m-%d %H:%M",
        format="[%(name)s][%(asctime)s][%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(),
        ]
    )
    logger.info("Arguments: {}".format(kwargs))

    assert int(kwargs['--threads']) > 0, "Number of threads must be greater than 0"

    kwargs['--output_dir'] = os.path.abspath(kwargs['--output_dir'])
    assert os.path.isdir(kwargs['--output_dir']), "Output directory does not exist: {}".format(kwargs['--output_dir'])

    file_extensions = {'tabular': '.tsv', 'genbank': '.gbk', 'fasta': '.fna'}
    assert kwargs['--out_format'] in file_extensions.keys(), "Output format must be one of: {}".format(file_extensions.keys())

    with open(kwargs['--input_file_list'], 'r') as f:
        input_files = f.read().splitlines()

    for input_file in input_files:
        assert os.path.exists(input_file), "Input file does not exist: {}".format(input_file)

    logger.info("Running phanotate on {} files".format(len(input_files)))

    logger.info("Starting phanotate annotations with {} threads".format(kwargs['--threads']))
    with Pool(int(kwargs['--threads'])) as p:
        p.map(
            partial(
                run_phanotate,
                output_dir=kwargs['--output_dir'],
                out_format=kwargs['--out_format'],
                input_file_extension=check_extensions(input_files[0]),
                output_file_extension=file_extensions[kwargs['--out_format']],
            ),
            input_files)

    logger.info("FINISHED !")

if __name__ == '__main__':
    main(**docopt(__doc__))
