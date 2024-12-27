import pandas as pd
import numpy as np
import pdb
import requests
import gzip
import os
import json
from collections import defaultdict


def get_ensmbl(dataset, batch_size=1000):
    gene_info = dataset['gene_info'].tolist()

    # Extract gene IDs
    filtered_gene_info = [gene.split(';')[0].split(' ')[1].split('.')[0][1:] for gene in gene_info]
    unique_gene_ids = list(set(filtered_gene_info)) 

    # Prepare batches
    batches = [unique_gene_ids[i:i+batch_size] for i in range(0, len(unique_gene_ids), batch_size)]

    # Annotation results
    annotation_dict = {}
    i = 1
    for batch in batches:
        print(f'batch {i} of {len(batches)}')
        url = "https://rest.ensembl.org/lookup/id"
        headers = {"Content-Type": "application/json"}
        response = requests.post(url, json={"ids": batch}, headers=headers)
        if response.ok:
            data = response.json()
            for gene_id in batch:
                if gene_id in data and data[gene_id] and 'description' in data[gene_id]:
                    annotation_dict[gene_id] = data[gene_id]['description']
                else:
                    annotation_dict[gene_id] = 'No annotations found'
        else:
            for gene_id in batch:
                annotation_dict[gene_id] = 'No annotations found'
        i += 1

    # Map back to dataset
    dataset['ensembl_annotation'] = dataset['gene_info'].apply(
        lambda x: annotation_dict[x.split(';')[0].split(' ')[1].split('.')[0][1:]]
    )
    dataset.to_csv('data/final_tr_annotations.tsv', sep='\t', index=False)

def create_full_dataset(dataset):
    # 1. extract gene ids
    motifs = dataset['motif'].tolist()
    starts = dataset['rep_start'].tolist()

    keys = [f"{x}_{y}" for x, y in zip(motifs, starts)]
    i = 0
    
    motif_id_map = {}
    with open('data/condensed_gene_ids.txt', 'r') as f:
        for i, (key, line) in enumerate(zip(keys, f)):
            motif_id_map[key] = line.strip()


    gene_exp_info = {}
    header = ''
    # Read the header separately to get column names
    with open('data/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct', 'r') as f:
        header = f.readline().strip().split('\t')

    # Initialize an empty dictionary for the tissue map
    tissue_map = {}

    # Read the sample attributes file
    attributes = pd.read_csv('data/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt', sep='\t')

    tissues_mapped = []
    # Process the file in chunks
    chunk_size = 10000  
    chunks = pd.read_csv(
        'data/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct',
        sep='\t',
        chunksize=chunk_size
    )
    j = 1
    for chunk in chunks:
        print(f'chunk {j} processing')
        j += 1
        cols = [col for col in chunk.columns if 'GTEX' in col]
        
        for c in cols:
            row = attributes[attributes['SAMPID'] == c]
            if not row.empty:
                name = row['SMTSD'].iloc[0]
                if 'Brain' not in name:
                    continue
                #temp = name.split('-')[0].strip() + '_' + name.split(' ')[2] + '_' + name.split(' ')[3]
                
                tissue = name.split('-')[0].strip()
                temp = name.split(' ')
                for s in temp:
                    if s not in tissue and s != '-' and '(' not in s and ')' not in s:
                        tissue += '_' + s
                
                if tissue not in tissues_mapped:
                    tissue_map[c] = tissue
                    tissues_mapped.append(tissue)
                
            
                
                
        


    
    print('extracting expression info')
    iteration = 0
    # Each gene in tpm file has expression data in each tissue/GTEX sample
    expression_data_cols = {}
    eqtl_files = os.listdir('data/GTEx_Analysis_v10_eQTL_updated')
    sqtl_files = os.listdir('data/sQTL_info')

    eqtl_dataframes = {}
    sqtl_dataframes = {}

    print('Reading tpm data')
    # columns_to_use = list(tissue_map.keys())
    # columns_to_use.append('Name')
    # tpm_data = pd.read_csv('data/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct', usecols=columns_to_use, sep='\t')
    expected_columns = header
    columns_to_use = [col for col in expected_columns if col in tissue_map.keys() or col == 'Name']

    
    tpm_data = pd.read_csv('data/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct', usecols=columns_to_use, sep='\t')
    
    # From these we need:
    # eqtl - variant_id, beta_shapes, pval_nominal, slope, qval, af, afc
    # sqtl - pip, variant_id 
    for sample in tissue_map:
        if ' ' not in tissue_map[sample]:
            tissue = tissue_map[sample][0].upper() + tissue_map[sample][1:]
        else:
            temp = tissue_map[sample].split(' ')
            tissue = ''
            for i in range(0, len(temp)):
                tissue += temp[i][0].upper() + temp[i][1:] + '_'
            tissue = tissue[:len(tissue) - 1]

        eqtl_list = [x for x in eqtl_files if tissue in x]
        if tissue not in eqtl_dataframes:
    
            for file in eqtl_list:
                if tissue not in eqtl_dataframes and 'txt' in file:
                    eqtl_dataframes[tissue] = pd.read_csv(f'data/GTEx_Analysis_v10_eQTL_updated/{file}', sep="\t", compression="gzip")
                    eqtl_dataframes[tissue].set_index("gene_id", inplace=True)
                
        sqtl_file = [x for x in sqtl_files if tissue in x]
        if tissue not in sqtl_dataframes:
            sqtl_dataframes[tissue] = pd.read_parquet(f'data/sQTL_info/{sqtl_file[0]}')
            sqtl_dataframes[tissue].set_index("gene_id", inplace=True)

    
    tissues = set(tissue_map.values())
    for tissue in tissues:
        for metric in ["eqtl_variant", "beta_shape1", "beta_shape2", "pval", "slope", "pfdr", "af", "afc", "pip", "sqtl_variant"]:
            expression_data_cols[f"{tissue}_{metric}"] = []


    gene_specific_files = {}
    
    for idx, row in dataset.iterrows():
        if iteration % 10000 == 0:
            print(f'step {iteration} out of {dataset.shape[0]}')
        curr_motif = row['motif']
        curr_start = row['rep_start']
        
        gene_id = motif_id_map[curr_motif + '_' + str(curr_start)]
        # expression_info = gene_exp_info[gene_id].pop(0)
        expression_info = tpm_data[tpm_data['Name'] == gene_id].to_dict()
        if gene_id not in gene_specific_files:
            gene_specific_files[gene_id] = []
        
        
        tissues_filled = []
        for key in tissue_map:
            if key not in expression_data_cols:
                expression_data_cols[key] = []

            # extract eqtl and sqtl info from files
            if ' ' not in tissue_map[key]:
                tissue = tissue_map[key][0].upper() + tissue_map[key][1:]
            else:
                temp = tissue_map[key].split(' ')
                tissue = ''
                for i in range(0, len(temp)):
                    tissue += temp[i][0].upper() + temp[i][1:] + '_'
                tissue = tissue[:len(tissue) - 1]
            
            
            
            if len(expression_info[key].keys()) > 1:
                print('Multiple keys')
                exit()
            
            if len(expression_info[key].keys()) == 0:
                expression_data_cols[key].append(0.0)
                if tissue not in tissues_filled:
                    expression_data_cols[f'{tissue}_eqtl_variant'].append(None)
                    expression_data_cols[f'{tissue}_beta_shape1'].append(None)
                    expression_data_cols[f'{tissue}_beta_shape2'].append(None)
                    expression_data_cols[f'{tissue}_pval'].append(None)
                    expression_data_cols[f'{tissue}_slope'].append(None)
                    expression_data_cols[f'{tissue}_pfdr'].append(None)
                    expression_data_cols[f'{tissue}_af'].append(None)
                    expression_data_cols[f'{tissue}_afc'].append(None)
                    expression_data_cols[f'{tissue}_pip'].append(None)
                    expression_data_cols[f'{tissue}_sqtl_variant'].append(None)
                    tissues_filled.append(tissue)
                continue

            index = list(expression_info[key].keys())[0]
            expression_data_cols[key].append(expression_info[key][index])
            
            if tissue in tissues_filled:
                continue
            # if it is expressed
            if expression_info[key][index] > 0:
                # extract eqtl and sqtl info from files   
                eqtl_row = eqtl_dataframes.get(tissue, pd.DataFrame())
                if gene_id not in eqtl_row.index:
                    eqtl_row = pd.DataFrame()
                else:
                    eqtl_row = eqtl_row.loc[gene_id]
                if eqtl_row.empty:
                    expression_data_cols[f'{tissue}_eqtl_variant'].append(None)
                    expression_data_cols[f'{tissue}_beta_shape1'].append(None)
                    expression_data_cols[f'{tissue}_beta_shape2'].append(None)
                    expression_data_cols[f'{tissue}_pval'].append(None)
                    expression_data_cols[f'{tissue}_slope'].append(None)
                    expression_data_cols[f'{tissue}_pfdr'].append(None)
                    expression_data_cols[f'{tissue}_af'].append(None)
                    expression_data_cols[f'{tissue}_afc'].append(None)
                else:
                    expression_data_cols[f'{tissue}_eqtl_variant'].append(eqtl_row.get("variant_id"))
                    expression_data_cols[f'{tissue}_beta_shape1'].append(eqtl_row.get("beta_shape1"))
                    expression_data_cols[f'{tissue}_beta_shape2'].append(eqtl_row.get("beta_shape2"))
                    expression_data_cols[f'{tissue}_pval'].append(eqtl_row.get("pval_nominal"))
                    expression_data_cols[f'{tissue}_slope'].append(eqtl_row.get("slope"))
                    expression_data_cols[f'{tissue}_pfdr'].append(eqtl_row.get("qval"))
                    expression_data_cols[f'{tissue}_af'].append(eqtl_row.get("af"))
                    expression_data_cols[f'{tissue}_afc'].append(eqtl_row.get("afc"))

                sqtl_row = sqtl_dataframes.get(tissue, pd.DataFrame())

                if gene_id not in sqtl_row:
                    sqtl_row = pd.DataFrame()
                else:
                    sqtl_row = sqtl_row.loc[gene_id]

                if sqtl_row.empty:
                    expression_data_cols[f'{tissue}_pip'].append(None)
                    expression_data_cols[f'{tissue}_sqtl_variant'].append(None)
                else:
                    expression_data_cols[f'{tissue}_pip'].append(sqtl_row.get("pip"))
                    expression_data_cols[f'{tissue}_sqtl_variant'].append(sqtl_row.get("variant_id"))

            else:
                expression_data_cols[f'{tissue}_eqtl_variant'].append(None)
                expression_data_cols[f'{tissue}_beta_shape1'].append(None)
                expression_data_cols[f'{tissue}_beta_shape2'].append(None)
                expression_data_cols[f'{tissue}_pval'].append(None)
                expression_data_cols[f'{tissue}_slope'].append(None)
                expression_data_cols[f'{tissue}_pfdr'].append(None)
                expression_data_cols[f'{tissue}_af'].append(None)
                expression_data_cols[f'{tissue}_afc'].append(None)
                expression_data_cols[f'{tissue}_pip'].append(None)
                expression_data_cols[f'{tissue}_sqtl_variant'].append(None)
            
            tissues_filled.append(tissue)
            

        iteration += 1
    

    # for gene in the file, we need to extract eqtl and sqtl info
    # then take that and add it to the dataset
    
    
    print('creating csv')
    for col in expression_data_cols:
        if 'GTEX' in col:
            col_name = f'{col}_({tissue_map[col]})'
            dataset[col_name] = expression_data_cols[col]
        else:
            dataset[col] = expression_data_cols[col]
    
    dataset.drop_duplicates(inplace=True)
    dataset.to_csv('data/expression_tr_annotations.txt.gz', sep='\t', index=False, compression="gzip")
    print('completed')
        

    # 5. extract tissue specific info from qtls 
    # 6. compile into one singular bed file


def parse_json(json_file, dataset):
    with open(json_file, 'r') as f:
        data = json.load(f)
        diseased_data = []
        for i in range(0, len(data)):
            chr_info = data[i]['LocationCoordinates']['hg38']
            chr = chr_info.split(':')[0]
            start_coord = chr_info.split('-')[0].split(':')[1]
            end_coord = chr_info.split('-')[1]
            
            
            motif = data[i]['Motif']
            id = list(data[i]['Diseases'].keys())[0]
            disease = data[i]['Diseases'][id]['DiseaseName']
            threshold = data[i]['Diseases'][id]['PathogenicCutoff']

            
            diseased_data.append({'chr': chr, 'start': start_coord, 'end': end_coord, 'motif': motif, 'disease': disease, 'threshold': threshold})
        

        data_found = 0
        for disease_motif in diseased_data:
            motif = disease_motif['motif']
            specific_rows_upper_bound = dataset[dataset['motif'] == motif]
            specific_rows_upper_bound = specific_rows_upper_bound[specific_rows_upper_bound['ann_start'] >= int(disease_motif['start'])]
            specific_rows_upper_bound = specific_rows_upper_bound[specific_rows_upper_bound['ann_start'] < int(disease_motif['end'])]
            if not specific_rows_upper_bound.empty:
                data_found += 1
                
            

        # df = pd.DataFrame(diseased_data)
        # df.to_csv('data/diseased_trs.bed', sep='\t')

def cross_reference(diseased_set, regular_set):
    regular_set.drop('ensembl_annotation', axis=1, inplace=True)
    condensed_list = []
    
    # regular_set['gene_info'] = regular_set['gene_info'].split(';')[0].split(' ')[1].split('"')[0].split('.')[0]
    # diseased_set['gene_info'] = diseased_set['gene_info'].split(';')[0].split(' ')[1].split('"')[0].split('.')[0]
    genes_found = []
    for idx, row in regular_set.iterrows():
        motif = row['motif']
        chr = row['chr']
        gene = row['gene_info'].split(';')[0].split(' ')[1].split('"')[0].split('.')[0]


        specific_rows = diseased_set[diseased_set['gene_info'].str.contains(gene)]
        if specific_rows.empty:
            condensed_list.append({'chr': chr, 'rep_start': row['rep_start'], 'rep_end': row['rep_end'], 'motif': motif, 'disease': 'N/A', 'threshold': 'N/A', 'ann_start': row['ann_start'], 'ann_end': row['ann_end'], 'gene_info': row['gene_info']})
        else:
            specific_rows = specific_rows[(specific_rows['motif'] == motif) & (specific_rows['chr'] == chr)]
            if specific_rows.empty:
                condensed_list.append({'chr': chr, 'rep_start': row['rep_start'], 'rep_end': row['rep_end'], 'motif': motif, 'disease': 'N/A', 'threshold': 'N/A', 'ann_start': row['ann_start'], 'ann_end': row['ann_end'], 'gene_info': row['gene_info']})
            else:
                for _, r in specific_rows.iterrows():
                    condensed_list.append({'chr': chr, 'rep_start': r['rep_start'], 'rep_end': r['rep_end'], 'motif': motif, 'disease': r['disease'], 'threshold': r['threshold'], 'ann_start': r['ann_start'], 'ann_end': r['ann_end'], 'gene_info': r['gene_info']})
                genes_found.append(gene)
    
    print(len(list(set(genes_found))))
    for _, row in diseased_set.iterrows():
        gene = row['gene_info'].split(';')[0].split(' ')[1].split('"')[0].split('.')[0]
        if gene not in genes_found:
            condensed_list.append({'chr': chr, 'rep_start': row['rep_start'], 'rep_end': row['rep_end'], 'motif': motif, 'disease': row['disease'], 'threshold': row['threshold'], 'ann_start': row['ann_start'], 'ann_end': row['ann_end'], 'gene_info': row['gene_info']})

    condensed_df = pd.DataFrame(condensed_list)
    condensed_df.to_csv('data/condensed_set.tsv', sep='\t')
    
    
    

if __name__ == "__main__":
    df = pd.read_csv('data/condensed_set.tsv', sep='\t')
    diseased_df = pd.read_csv('data/diseased_intersect.tsv', sep='\t')
    diseased_df.columns = ['chr', 'rep_start', 'rep_end', 'motif', 'disease', 'threshold', 'chr_ref', 'ann_start', 'ann_end', 'gene_info']
    create_full_dataset(df)
    #parse_json('data/diseased_str_catalog.json', df)
    #cross_reference(diseased_df, df)

    