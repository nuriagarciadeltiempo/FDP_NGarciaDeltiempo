
import pandas as pd

#import csv omim file as data frame
OMIM_file = pd.read_csv(r'OMIM_X-linked_non_adult_onset.csv', sep=';')

#import txt with the gene symbol from Genome Project database
annotated_genes = open('GP_SYMBOL.txt', 'r')

# create a list with the gene symbol from Genome Project database
gene_GP_list = []
for gene in annotated_genes:
    if gene not in gene_GP_list:
        gene_GP_list.append(gene)
    else:
        continue

# remove per each element the \n character
gene_GP_list = [x.strip() for x in gene_GP_list]

#create a file with gene_GP_list
with open('gene_GP_list.txt', 'w') as f:
    for item in gene_GP_list:
        #write them separated by a new line
        f.write("%s \n" % item)

# for each gene symbol in the list, check if it is present in the OMIM file
gene2disease = {}
for gene in gene_GP_list:
    for Gene in OMIM_file['Gene']:
        if gene == Gene:
            # add to the dictionary as a key the gene symbol and as value the associated disease
            gene2disease[gene] = OMIM_file['Associated Disease'][OMIM_file['Gene'] == Gene].values
        else:
            continue
        
#convert teh dictionary into a dataframe
dictionary_dataframe = pd.DataFrame.from_dict(gene2disease, orient='index', columns=['Associated Disease'])
#store dictionary_dataframe into a txt file
dictionary_dataframe.to_csv('gene2disease_df.txt', sep='\t')
