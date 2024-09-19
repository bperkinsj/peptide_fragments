import os
from pathlib import Path

def add_cmd_line_args(parser):
    parser.add('-c', '--config', required=False, is_config_file=True)

    parser.add_argument('--data_dir', type=str, help='Directory containing the output data files')
    parser.add_argument('--keep_structures', action='store_true', help='Remove AlphaFold2 structure files after processing')
    parser.add_argument('--input', type=str, help='Input data, either a file or a directory')
    parser.add_argument('--dataset_name', type=str, help='Name of the dataset')
    parser.add_argument('--is_file', action='store_true', help='Denotes that the input is a file')

    parser.add_argument('--verbose', type=bool, help='Print verbose output', default=False)


    args = parser.parse_args()
    config = vars(args)
    return config

def add_hardcoded_args(config):
    
    # File columns
    config['uniprot_column'] = 'PG.ProteinAccessions'
    config['protein_name_column'] = 'PG.ProteinNames'
    config['pep_pos_column'] = 'PEP.PeptidePosition'
    config['seq_column'] = 'PEP.StrippedSequence'
    config['pval_column'] = 'pval'
    config['dif_column'] = 'dif'
    config['id_column'] = 'PEP.PeptideID'
    config['file_name_column'] = 'file_name'
    config['region'] = 'PEP.Region'
    config['alpha_helix'] = 'alpha_helix'
    config['beta_sheet'] = 'beta_sheet'
    config['bend'] = 'bend'
    config['domains'] = 'PG.Domains'

    # filename strings
    config['uniprot_acc_fn_str'] = 'uniprot_ids.npy'
    config['pkl_fn_str'] = '.pkl'
    config['csv_fn_str'] = '.csv'

    # Directory names
    config['structure_dir_str'] = 'structures'
    config['cif_dir_str'] = 'cif_files'
    config['pkl_dir_str'] = 'pkl'
    config['uniprot_dir_str'] = 'uniprots'
    config['results_dir_str'] = 'results'


def add_path(config):
    dataset_dir = os.path.join(config['data_dir'], config['dataset_name'])
    uniprot_dir = os.path.join(dataset_dir, config['uniprot_dir_str'])

    # Add directories
    config['structure_dir'] = os.path.join(dataset_dir, config['structure_dir_str'])
    config['struct_cif_dir'] = os.path.join(config['structure_dir'], config['cif_dir_str'])
    config['pkl_dir'] = os.path.join(dataset_dir, config['pkl_dir_str'])
    config['results_dir'] = os.path.join(dataset_dir, config['results_dir_str'])

    # Add filenames
    config['uniprot_acc_fn'] = os.path.join(uniprot_dir, config['uniprot_acc_fn_str'])

    #Create dir
    for dir in [dataset_dir, uniprot_dir, config['structure_dir'], config['struct_cif_dir'], config['pkl_dir'], config['results_dir']]:
        if not os.path.exists(dir):
            Path(dir).mkdir(parents=True, exist_ok=True)

def parse_args(parser):
    print('=== Parsing ===')
    config = add_cmd_line_args(parser)
    add_hardcoded_args(config)
    add_path(config)
    return config