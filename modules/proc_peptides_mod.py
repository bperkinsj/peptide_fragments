import os
import pickle
import shutil
import scripts.utils as utils
import pandas as pd
import numpy as np


class pipeline:
    def __init__(self, option, args):
        self.verbose = args.verbose
        self.results_dir = args.results_dir

        # File columns
        self.uniprot_column = args.uniprot_column
        self.protein_name_column = args.protein_name_column
        self.pep_pos_column = args.pep_pos_column
        self.seq_column = args.seq_column
        self.pval_column = args.pval_column
        self.dif_column = args.dif_column
        self.id_column = args.id_column
        self.file_name_column = args.file_name_column
        self.region = args.region
        self.alpha_helix = args.alpha_helix
        self.beta_sheet = args.beta_sheet
        self.bend = args.bend
        self.domains = args.domains

        if option == 'init_input_procs':
            self._init_files(args)
            self._init_download(args)
            self._init_input_procs(args)
        elif option == 'delete':
            self._init_delete_files(args)
        else:
            raise Exception('Invalid option')


    ######
    # inits
    def _init_files(self, args):
        '''
        Initialize the files to process
        '''
        self.is_file = args.is_file
        self.input = args.input
        if self.is_file:
            self.files_to_process = [self.input]
        else:
            self.files_to_process = os.listdir(self.input)

    def _init_download(self, args):
        self.structure_dir = args.structure_dir
        self.struct_cif_dir = args.struct_cif_dir


    def _init_input_procs(self, args):

        # Pickling files
        self.pkl_fn_str = args.pkl_fn_str
        self.pikl_dir = args.pkl_dir
        self.pikled_files = []

        # Domain files
        self.csv_fn_str = args.csv_fn_str


    def _init_delete_files(self, args):
        self.structure_dir = args.structure_dir

    
    ######
    # processors
    def download_structures(self):
        '''
        Download the AlphaFold2 structures from the AlphaFold Protein Structure Database
        '''

        print('Downloading structures...')
        if self.verbose:
            print(f'Download directory: {self.structure_dir}')

        # Collect uniprot ids and fetch structures
        self.collect_uniprot_accs()

        self.fetch_structures()

    def process_input(self):
        '''
        Extract information for each peptide from the AlphaFold files and from the UniProt database.
        '''
        print('Processing input...')
        if self.verbose:
            print(f'Peptides directory: {self.input}')
            print(f'Structures directory: {self.structure_dir}')


        self.collect_uniprot_accs()
        
        # Get secondary structure information for each fragment
        self.secondary_structure()

        # Get family and domain information for each fragment
        self.domain_info()

    def delete_files(self):
        '''
        Delete the structure files once processing is complete
        '''
        if self.verbose:
            print(f'Structure directory: {self.structure_dir}')

        shutil.rmtree(self.structure_dir)

    ######
    # helpers
    def collect_uniprot_accs(self):
        '''
        Collect all unique UniProt accessions from the files to process
        '''
        uniprot_accs = []

        for f in self.files_to_process:
            fn = os.path.join(self.input, f)
            df = pd.read_csv(fn, sep='\t')
            uniprot_accs += df[self.uniprot_column].unique().tolist()
        
        self.uniprot_accs = uniprot_accs

    def collect_uniprot_fns(self):
        '''
        Collect the file names for the AlphaFold2 structures and convert to dataframe format
        '''
        uniprot_fns = {}
        files = os.listdir(self.struct_cif_dir)
        for f in files:
            uniprot = f.split('.')[0]
            uniprot_fns[uniprot] = f

        uniprot_fns = pd.DataFrame(uniprot_fns.items(), columns=[self.uniprot_column, self.file_name_column])
        self.uniprot_fns = uniprot_fns

    def fetch_structures(self):
        '''
        Download all AlphaFold2 structures for the given UniProt accessions 
        '''
        # Check if files already exist to save time
        files_exist = os.listdir(self.struct_cif_dir)
        print(self.struct_cif_dir)
        if len(files_exist) == 0:
            utils.download_alphafold_structures(self.uniprot_accs, self.structure_dir)

        self.collect_uniprot_fns()

    def secondary_structure(self):
        for f in self.files_to_process:
            fn = os.path.join(self.input, f)
            df = pd.read_csv(fn, sep='\t')
            df = pd.merge(df, self.uniprot_fns, on=self.uniprot_column, how='left')
            # Remove any rows that do not have existing files
            df = df.dropna(subset=[self.file_name_column]).reset_index(drop=True)

            # Add peptide ID and region columns
            df[self.id_column] = df[self.uniprot_column] + '-' + df[self.pep_pos_column].astype(str) + '-' + str(len(df[self.seq_column]))
            df = utils.add_region(self, df)

            results_df = utils.pep_info(self, self.struct_cif_dir, df)
            pkl_fn = f.split('.')[0] + self.pkl_fn_str
            fp = os.path.join(self.pikl_dir, pkl_fn)
            results_df.to_pickle(fp)
            self.pikled_files.append(fp)

    def domain_info(self):
        print('Matching peptides to domains...')

        for f in self.pikled_files:
            df = pd.read_pickle(f)
            df = utils.match_peps_to_domains(self, df)
            f = f.split('/')[-1]
            fn = f.split('.')[0] + self.csv_fn_str
            fp = os.path.join(self.results_dir, fn)
            df.to_csv(fp, sep=',', index=False)