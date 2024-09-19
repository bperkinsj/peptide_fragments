import argparse
import configargparse

from parser import parse_args
from modules import proc_peptides_mod

def process_input(args):
    '''
    Extract information for each peptide from the AlphaFold files and from the UniProt database.
    '''
    pipeline = proc_peptides_mod.pipeline('init_input_procs', args)
    pipeline.download_structures()
    pipeline.process_input()

def delete_files(args):
    '''
    Delete the AlphaFold2 structure files after processing.
    '''
    pipeline = proc_peptides_mod.pipeline('delete', args)
    pipeline.delete_files()

def main():
    parser = configargparse.ArgumentParser()
    config = parse_args(parser)
    args = argparse.Namespace(**config)

    process_input(args)

    if not args.keep_structures:
        print('Deleting structure files...')
        delete_files(args)

if __name__ == '__main__':
    main()