#!/usr/bin/env python3
"""
prokka runner.
Usage:
    prokka_runner.py ( --input_file_list=PATH ) ( --output_dir=PATH )
                      ( --proteins=PATH ) [ --threads=INT ] [ --prokka_threads=INT ]

Options:
    -h --help                      Show this screen.
    -i --input_file_list=PATH      Path to a file containing a list of input fasta files.
    -o --output_dir=PATH           Path to the output directory.
    -d --proteins=PATH               Path to the prokka database directory.
    -t --threads=INT               Number of threads to use for running individual anotations [default: 1].
    -p --prokka_threads=INT        Number of threads to use for running prokka [default: 1].
"""
import os
import sys
import subprocess
import logging
from multiprocessing import Pool
from functools import partial
import re
from docopt import docopt
logger = logging.getLogger(__name__)

def run_prokka(input_file, output_dir, proteins, prokka_threads):

    prefix = re.findall('S\d+#1[^.]+', os.path.basename(input_file))[0]

    sample_output_dir = os.path.join(output_dir, prefix)

    logger.info("Running prokka for {} with prefix {} ...".format(input_file, prefix))

    cmd = [
        'prokka',
        '--outdir', sample_output_dir,
        '--cpus', str(prokka_threads),
        '--increment', '10',
        '--compliant',
        '--prefix', prefix,
        '--proteins', proteins,
        '--addgenes',
        input_file
    ]

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        logger.error("Error while running prokka on {}".format(input_file))
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

    logger.info("args: {}".format(kwargs))

    # assert threads are integers greater than 0
    for i in ('--threads', '--prokka_threads'):
        if i in kwargs:
            assert int(kwargs[i]) > 0, "--{} must be an integer greater than 0".format(i)
            kwargs[i] = int(kwargs[i])


    # Cheking all input files exist
    for i in ('--input_file_list', '--output_dir', '--proteins'):
        assert os.path.exists(kwargs[i]), '{} does not exist'.format(kwargs[i])
        # make absolute path
        kwargs[i] = os.path.abspath(kwargs[i])

    # Reading input file list
    with open(kwargs['--input_file_list'], 'r') as f:
        input_files = f.read().splitlines()


    parallel_processes = max(1, int(kwargs['--threads']/kwargs['--prokka_threads']))

    # Running prokka
    logger.info("Running {} parallel processes with {} threads each".format(parallel_processes, kwargs['--prokka_threads']))

    with Pool(parallel_processes) as p:
        p.map(
            partial(
                run_prokka,
                output_dir=kwargs['--output_dir'],
                proteins=kwargs['--proteins'],
                prokka_threads=kwargs['--prokka_threads']
            ),
            input_files
        )

    logger.info("FINISHED !")

if __name__ == '__main__':
    main(**docopt(__doc__))
