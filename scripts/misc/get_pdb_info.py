#!/usr/bin/env python
"""
This script pulls information from the PDB about the number of results with experimental data
specified. Counts entities/unique chains, not structures/PDB entries.
"""

import argparse
import requests
import pandas as pd

def get_args():
    """ Gather input arguments from command line """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', metavar='INPUT_FILE', required=True,
                        help='Input file with protein information')
    parser.add_argument('-e', '--exp_type', metavar='EXPERIMENT_TYPE', nargs='+', 
                        required=True, help='The type of experiment to check for')
    parser.add_argument('-o', '--output_file', metavar='OUTPUT_FILE',
                        help='Output file to write protein information to. If not \
                        specified it will overwrite the input file.')
    args = parser.parse_args()

    return args

def send_request(uniprot_id, experiment_type):
    """ Submit the query request to PDB and return the number of results """
    pdb_rest_url = "http://www.rcsb.org/pdb/rest/search"
    query_text = """

<orgPdbCompositeQuery version="1.0">
 <queryRefinement>
  <queryRefinementLevel>0</queryRefinementLevel>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>
    <description>Simple query for a list of Uniprot Accession IDs: {0} </description>
    <accessionIdList>{0}</accessionIdList>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>1</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ExpTypeQuery</queryType>
    <description>Experimental Method is {1} and has Experimental Data</description>
    <mvStructure.expMethod.value>{1}</mvStructure.expMethod.value>
    <mvStructure.hasExperimentalData.value>Y</mvStructure.hasExperimentalData.value>
  </orgPdbQuery>
 </queryRefinement>
</orgPdbCompositeQuery>

"""
    query_text = query_text.format(uniprot_id, experiment_type)
    req = requests.post(url=pdb_rest_url,
                        data=query_text.encode(),
                        headers={'Content-Type': 'application/x-www-form-urlencoded'})

    if not req.status_code == 200:
        raise requests.HTTPError("Error: {} on response".format(req.status_code))

    return req.text.strip().split('\n')

def main():
    """ Main routine for downloading information per protein """

    args = get_args()

    prot_df = pd.read_csv(args.input_file, sep='\t', header=0)
    prot_df = prot_df.reindex(columns=(prot_df.columns.tolist() + args.exp_type))
    prot_ids = prot_df.loc[:,"UniProtID"]

    count = 0
    for pid in prot_ids:
        for exp in args.exp_type:
            req_list = send_request(pid, exp)
            prot_df.loc[count, exp] = len(req_list)
        count += 1

    if args.output_file:
        prot_df.to_csv(args.output_file, sep='\t')
    else:
        prot_df.to_csv(args.input_file, sep='\t')

if __name__ == '__main__':
    main()
