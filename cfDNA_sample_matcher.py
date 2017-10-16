#!/usr/bin/env python3
#Pavlos Antoniou
#10/10/17
#Canceer Analyst
#GeL

import pandas as pd
import json
import numpy as np
from pandas.io.json import json_normalize

# get list of json files: which has been collated using bash and awk - n.b. these should be latest delivery/pipeline
inputcsv = "/Users/johnambrose/Documents/test_json_list.csv"
reportcsv = "/Users/johnambrose/Documents/reports_log.20092017.TEST.csv"
cfdnacsv = "/Users/johnambrose/Documents/cfdna_samples.csv"
outfile = "/Users/johnambrose/Documents/json_gene_matching_list.csv"

# get list of genes in gene panel (supplied by Emma Walsh)
with open('/Users/johnambrose/Documents/test_gene_panel', 'r') as thermo_list:
    thermo_genes = thermo_list.read().splitlines()
    # print(thermo_genes)

def get_genes(df):
    genelist = []
    for d in df['reportedVariantCancer.reportEvents']:
        genelist.append(d[0]['genomicFeatureCancer']['geneName'])

    genes = np.asarray(genelist)
    genes_unique = np.unique(genes)
    return genes_unique

def check_sample_for_genes(json_files_df):

    sample_id = json_files_df['sample_id']
    json_filename = json_files_df['jsonpath']
    lab_sample_id = json_files_df['LAB_SAMPLE_ID']

    with open(json_filename, 'r') as json_file_to_be_read:
        json_file = json.load(json_file_to_be_read)
        json_df = pd.io.json.json_normalize(json_file)

    json_df_som = json_df.loc[json_df['somaticOrGermline'] == 'somatic']

    if any(i in thermo_genes for i in get_genes(json_df_som)):
        gene_match = True
    else:
        gene_match = False

    d = [str(sample_id), str(lab_sample_id), str(gene_match)]

    return d

# get locations of json files
datalocations = pd.read_csv(inputcsv)

# get report csv with key to lab ID and plate ID (for json matching)
report = pd.read_csv(reportcsv)

merged_json_and_report = pd.merge(datalocations, report, how='inner', left_on='sample_id', right_on='SAMPLE', indicator=False)

json_with_genes_matched = merged_json_and_report.apply(check_sample_for_genes, axis=1)

json_with_genes_matched_df = pd.DataFrame.from_items(zip(json_with_genes_matched.index, json_with_genes_matched.values)).T
json_with_genes_matched_df.columns = ["sample_id", "lab_sample_id", "gene_match"]

# read list of patients with cfDNA samples supplied on Emma Walsh spreadsheet
cfdnas = pd.read_csv(cfdnacsv)
cfdnas = cfdnas.applymap(str)

merged_json_and_cfdna = pd.merge(cfdnas, json_with_genes_matched_df, how='left', left_on='Laboratory_Sample_ID', right_on='lab_sample_id', indicator=True)

#merged_report_and_json = pd.merge(report, json_gene_matching_list, how='inner', left_on='SAMPLE', right_on='sample_id', indicator=False)
# print(merged_report_and_json.loc[merged_report_and_json['matching_genes'] == 'False'])
# print(merged_report_and_json)

header = ["Participant_ID", "Laboratory_Sample_ID", "Cancer_Type", "lab_sample_id", "gene_match", "_merge"]

merged_json_and_cfdna.to_csv(outfile, index=False, sep='\t', columns=header, na_rep='NA')

