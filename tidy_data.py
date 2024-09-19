import argparse
import os
import pandas as pd

def fix_peps(df):
    '''
    Take any peptide start position with multiple values and make them into distinct rows.
    '''

    df['PEP.PeptidePosition'] = df['PEP.PeptidePosition'].astype('string')

    df['PEP.PeptidePosition'] = df['PEP.PeptidePosition'].str.split(',')

    df['PEP.PeptidePosition'] = df['PEP.PeptidePosition'].str.split(';')

    df['PG.ProteinAccessions'] = df['PG.ProteinAccessions'].str.split(';')

    df = df.explode('PEP.PeptidePosition').reset_index(drop=True)

    df = df.explode('PG.ProteinAccessions').reset_index(drop=True)

    return df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', type=str, help='Directory with files to tidy')
    args = parser.parse_args()

    files = os.listdir(args.dir)

    for f in files:
        fn = os.path.join(args.dir, f)
        df = pd.read_csv(fn, sep='\t')
        df = fix_peps(df)
        df.to_csv(fn, sep='\t', index=False)

if __name__ == '__main__':
    main()