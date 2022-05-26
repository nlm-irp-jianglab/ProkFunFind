from os import path
import subprocess

def run_emapper(config: dict, genome_in_path: str):
    """
    """
    outpath = genome_in_path.split('/')
    genome_prefix = outpath[len(outpath)-1]
    outdir = outpath[0:len(outpath)-2].join('/')

    cmd = ['emapper.py -i ', genome_in_path, '.pff.faa -o ', genome_in_path,
           '.emapper.annotations', ' --temp_dir ', config['main']['tmp_dir'], ' --tax_scope ',
           config['emapper'].get('tax_scope', 'bacteria'), ' --override ', '--cpu ',
           config['main'].get('cpus', 1), ' --data_dir ', config['emapper']['datadir']]
    res = subprocess.run(cmd)
    if res.returncode != 0:
        raise RuntimeError("Failed to run: {}".format(" ".join(cmd)))


def run_kofamscan(config: dict, genome_in_path):
    cmd = ['exec_annotation ', '--tmp-dir ', config['main']['tmp_dir'], ' -f detail-tsv ',
           '-p ', config['kofamscan']['profile'], ' -o ', genome_in_path, '.kofam.tsv',
           ' ', genome_in_path, '.pff.faa']
    res = subprocess.run(cmd)
    if res.returncode != 0:
        raise RuntimeError("Failed to run: {}".format(" ".join(cmd)))


def run_interproscan(config: dict, genome_in_path):
    cmd = ['interproscan ', '-t p ', '--goterms ', '--pathways ', '-f tsv',
           genome_in_paths, '.pff.faa']
    res = subprocess.run(cmd)
    if res.returncode != 0:
        raise RuntimeError("Failed to run: {}".format(" ".join(cmd)))
    # file_dir = path.dirname(path.abspath(__file__))
    # file_dir += '/../../scripts/tsv2ipr_gff.py'

    # cmd2 = ['python ', filedir, genome_in_path, '.pff.gff3 ', genome_in_path,
    # '_InterProScan.tsv ', '> ', genome_in_path, '.ipr-track.gff3']
    # res = subprocess.run(cmd2)
    # if res.returncode != 0:
    #     raise RuntimeError("Failed to run: {}".format(" ".join(cmd2)))
