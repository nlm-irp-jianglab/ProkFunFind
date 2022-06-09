from os import path
import subprocess

def run_emapper(config: dict, genome_in_path: str):
    """
    """
    outpath = genome_in_path.split('/')
    genome_prefix = outpath[len(outpath)-1]
    outdir = '/'.join(outpath[0:len(outpath)-2])
    exec_path = config['emapper'].get('exec', '')+'emapper.py'
    cmd = [exec_path, '-i', genome_in_path+'.pff.faa', '-o', genome_in_path+'.emapper.annotations', 
           '--temp_dir', config['main']['tmp_dir'], '--tax_scope',
           config['emapper'].get('tax_scope', 'bacteria'), '--override', '--cpu',
           config['main'].get('cpus', 1), '--data_dir', config['emapper']['datadir']]
    res = subprocess.run(cmd)
    if res.returncode != 0:
        raise RuntimeError("Failed to run: {}".format(" ".join(cmd)))


def run_kofamscan(config: dict, genome_in_path):
    exec_path = config['kofamscan'].get('exec', '')+'exec_annotation'
    cmd = [exec_path, '--cpu', config['main'].get('cpus', 1), '--tmp-dir', config['main']['tmp_dir'], '--format', 'detail-tsv',
           '-p', config['kofamscan']['profile'], '-o', genome_in_path+'.kofam.tsv',
           genome_in_path+'.pff.faa']
    print(cmd)
    res = subprocess.run(cmd)
    if res.returncode != 0:
        raise RuntimeError("Failed to run: {}".format(" ".join(cmd)))


def run_interproscan(config: dict, genome_in_path):
    exec_path = config['interproscan'].get('exec', '')+'interproscan'
    cmd = [exec_path, '-t', 'p', '--goterms', '--pathways', '-f', 'tsv',
           genome_in_path+'.pff.faa']
    res = subprocess.run(cmd)
    if res.returncode != 0:
        raise RuntimeError("Failed to run: {}".format(" ".join(cmd)))
