'''
Import dependencies.
'''
import os
import requests
import pandas as pd
import numpy as np
from pymol import cmd, stored
from alphafetcher import AlphaFetcher


def download_alphafold_structures(uniprot_accs, output_dir):
    '''
    Download the AlphaFold2 structures from the AlphaFold Protein Structure Database for
    the given UniProt accessions.
    '''
    fetcher = AlphaFetcher(base_savedir=output_dir)

    # Add UniProt access codes
    fetcher.add_proteins(uniprot_accs)

    # Retrieve metadata
    fetcher.fetch_metadata(multithread=True, workers=4)

    # Download the structures
    fetcher.download_all_files(pdb=False, cif=True, multithread=True, workers=4)

def pep_info(self, struct_dir, df):

    # Dictionary to store results
    results_dict = {}

    # Iterate through each row in the dataframe
    for index, row in df.iterrows():
        # Get the fragment id, uniprot id, and region
        frag_id = row[self.id_column]
        region = row[self.region]
        fn = row[self.file_name_column]
        fp = os.path.join(struct_dir, fn)

        # Run PyMOL
        # Clear the stored data
        cmd.delete('all')
        # Load the structure
        cmd.load(fp)
        stored.secondary_structure_list_by_resnumber = []
        cmd.iterate(f'(chain A and resi {region} and name ca)', 'stored.secondary_structure_list_by_resnumber.append((resv,ss))', quiet=1)
        ss = list(stored.secondary_structure_list_by_resnumber)

        # Categorize residues into secondary structures
        # result = {'alpha_helix': [], 'beta_sheet': [], 'bend': []}
        helix_res = []
        sheet_res = []
        bend_res = []
        # Add residues to lists
        for res in ss:
            if 'H' in res:
                helix_res.append(res[0])
            elif 'S' in res:
                sheet_res.append(res[0])
            elif 'B' in res:
                bend_res.append(res[0])

        # Convert lists to strings
        helix_res_str = ','.join(str(x) for x in helix_res)
        sheet_res_str = ','.join(str(x) for x in sheet_res)
        bend_res_str = ','.join(str(x) for x in bend_res)
        # Create a dictionary to store the results
        result = {'alpha_helix': helix_res_str, 'beta_sheet': sheet_res_str, 'bend': bend_res_str}

        results_dict[frag_id] = result

    # Add results to the DataFrame
    df[self.alpha_helix] = df[self.id_column].map(lambda x: results_dict.get(x, {}).get('alpha_helix', ''))
    df[self.beta_sheet] = df[self.id_column].map(lambda x: results_dict.get(x, {}).get('beta_sheet', ''))
    df[self.bend] = df[self.id_column].map(lambda x: results_dict.get(x, {}).get('bend', ''))

    df = df.reset_index(drop=True)

    return df

def string2range(x):
    
    """
    This function takes in a `string` representing a region of interest in a
    protein. The region of interest can be a single region or multiple regions
    of a protein. Returns a range for single regions or a list of ranges for
    multiple regions.
    
    Parameters:
    
        x (string): String containing a region or several regions of interest in a 
            protein.
            Format of x: single region -> 'start-end'
                         multiple regions -> 'start1-end1,start2-end2'
                     
    Returns:
    
        range or list of ranges: For single region proteins a range is returned. For 
            multiple region proteins a list of ranges is returned

            Format: single region -> range(start, end+1)
                    multiple region -> [range(start1, end1+1), range(start2, end2+1)]
    """
    # Handle instances with more than one range
    if ',' in x:
        list_temp = x.split(sep = ',') #list_temp = ['123-456,' '789-1111']
        for y in range(len(list_temp)): 
            list_temp[y] = list_temp[y].split(sep = '-') #list_temp[y] = [['123', '456'], ['789', '1111']]
        for y in range(len(list_temp)): 
            for x in range(len(list_temp[y])):
                list_temp[y][x] = int(list_temp[y][x]) #turns each list item into an integer

        # Make a range object with the bounds of the range. Note to the 
        # end a 1 has to be added in order to include the last position in the range
        for y in range(len(list_temp)): #[1, 2] where 1=[123, 456] and 2=[789, 1111]
            for x in range(len(list_temp[y])): #[123, 456]       
                list_temp[y] = list(range(list_temp[y][x], list_temp[y][x+1]+1)) #list_temp[0][0] = [123], list_temp[0][0+1]+1 or [456] + 1 = [457]
                break

        return list(set([item for sublist in list_temp for item in sublist]))

    # Handle instances with only one range
    else:
        list_temp = x.split(sep = '-')
        for y in range(len(list_temp)):
            list_temp[y] = int(list_temp[y]) #

        # Make a range object with the bounds of the region. Note to the 
        # end a 1 has to be added in order to include the last position in the range
        return list(range(list_temp[0], list_temp[1]+1))
    
def seq_intersect(x, y):
    '''
    Check whether two ranges intersect
    '''

    # If the range is a single residue, convert it to a range
    if len(x) == 1:
        x = [x[0], x[0]]
    if len(y) == 1:
        y = [y[0], y[0]]
    
    intersect = range(max(x[0], y[0]), min(x[1], y[1]) + 1)
    return len(intersect) > 0

def get_domains(uniprot_id): # TODO: add this to a step that checks whether this uniprot id is the same as the previous one
    '''
    Get the domain information from UniProtKB
    '''
    fields = ('ft_domain,ft_coiled,ft_compbias,cc_domain,ft_motif,protein_families,ft_region,ft_repeat,ft_zn_fing')
    url = f'https://rest.uniprot.org/uniprotkb/search?query=accession:{uniprot_id}&fields={fields}'
    response = requests.get(url)
    response_dic = response.json()
    domains = {}
    try:
        features = response_dic['results'][0]['features']
        # Get the start and end of any domains
        for i in range(len(features)):
            start = str(features[i]['location']['start']['value'])
            end = str(features[i]['location']['end']['value'])
            region = f'{start}-{end}'
            name = features[i]['type']
            description = features[i]['description']

            domains[name] = {'region': region, 'description': description}

    except KeyError:
        print(f'No domains found for {uniprot_id}')
        domains = None

    return domains

def flatten_dict(data):
    '''
    Flatten a nested dictionary
    '''

    out = {}
    out['type'] = list(data.keys())[0]
    values = list(data.values())
    out['region'] = values[0]['region']
    out['description'] = values[0]['description']

    return out

def match_peps_to_domains(self, df):
    current_id = None
    domains = None

    domain_dict = {}

    for index, row in df.iterrows():
        uniprot_id = row[self.uniprot_column]
        seq_region = row[self.region]
        seq_range = string2range(seq_region)

        # If the current uniprot id is different from the previous one, get the domains
        if uniprot_id != current_id:
            domains = get_domains(uniprot_id)
            current_id = uniprot_id

        domains_to_keep = []
        # Check whether there are any domains
        if domains:
            # Check whether the peptide intersects with any of the domains
            for domain, value in domains.items():
                region = value['region']
                domain_range = string2range(region)
                if seq_intersect(seq_range, domain_range):
                    domains_to_keep.append({domain: value})
                    

        domain_dict[row[self.id_column]] = domains_to_keep
    
    # Add the domain information to the DataFrame
    df[self.domains] = df[self.id_column].map(lambda x: domain_dict.get(x, {}))

    # Drop any rows with empty dictionaries
    df[self.domains] = df[self.domains].replace('{}', np.nan)
    df = df.dropna(subset=[self.domains]).reset_index(drop=True)

    # Explode the domain dictionaries into separate rows
    df = df.explode(self.domains).dropna(subset=[self.domains]).reset_index(drop=True)

    # Flatten the dictionaries
    df[self.domains] = df[self.domains].map(flatten_dict)

    # Expand the dicts into separate columns
    df2 = pd.json_normalize(df[self.domains]).reset_index(drop=True)

    df = pd.merge(df, df2, how='left', left_index=True, right_index=True)

    # Drop the interim domains column
    df = df.drop(columns=[self.domains])

    return df

def add_region(self, df):
    '''Add region column to dataframe by taking peptide start position and adding length of peptide to it'''
    for index, row in df.iterrows():
        start_pos = row[self.pep_pos_column]
        seq_len = len(row[self.seq_column])
        end_pos = int(start_pos) + seq_len - 1
        region = f'{start_pos}-{end_pos}'
        df.loc[index, self.region] = region

    return df

def filter_domains(self, df):
    '''
    Filter the dataframe to only include rows with domains
    '''
    df = df[df[self.domains].map(lambda x: len(x) > 0)]
    df = df.reset_index(drop=True)
    df2 = pd.json_normalize(df[self.domains])

    return df