
import pandas as pd
import sys
import re
import csv

# read the vcf file
vcf_file = csv.reader(open (sys.argv[1], "r"), delimiter="\t") 


# create a dictionary with the information from gene2disease_df to have an easy access.
gene2disease_dict = {}
# include in the dictionary as keys the gene names and as values the recessive/dominant information
with open("OMIM_X-linked_non_adult_onset_removedcol.csv", "r") as gene2disease_df:
    for line in gene2disease_df:
        line = line.split(";")
        gene2disease_dict[line[1]] = line[2]

ID2gender = {}
with open("ID2gender.txt", "r") as id2gender_df:
    for line in id2gender_df:
        line = line.split(":")
        ID2gender[line[0]] = line[1]

# filter the vcf file to keep only the variants related to the gene present in the gene2disease_df.txt file
#MAIN 
HEADER_dict= dict()
variant_gene_dict= dict()
for line in vcf_file:
    if line[0].startswith("#"):
        # line[0] is considered to be the header
        if line[0].startswith("##INFO=<ID=CSQ"): # when in the header finds the line with the format of the CSQ field
            format = re.search(r'Format: ([^;]*)\">', line[0]).group(1)
            format_list = format.split("|") # split the format to be more readable
            #include a last element to the list called "HERITANCE"
            format_list.append("HERITANCE")
            for i in range(0, len(format_list)): # create a dictionary with the format as key and the index as value
                HEADER_dict[format_list[i]]=format_list.index(format_list[i])
        elif line[0].startswith("#CHROM"):
            header_list = []
            for i in range(9, len(line)):
                #print(line[i]) print the IDs of the samples
                for line[i] in ID2gender.keys():
                    line[i] = line[i].replace(line[i], ID2gender[line[i]])
                    header_list.append(line[i])
            print('\t'.join(line + header_list))
            
                    

    else:
        format_info = line[7].split('CSQ=')
        final_format_info = format_info[0]
        format_info_all_transcript = format_info[1].split(",")
        n=0 
        for format_info_per_transcript in format_info_all_transcript:
            format_info_list = format_info_per_transcript.split("|")
            format_info_list.insert(HEADER_dict["HERITANCE"],"")
            gene = format_info_list[HEADER_dict["SYMBOL"]]
            
            for Gene in gene2disease_dict.keys(): 
                if gene == Gene: # in case the gene present in the gene2disease.txt is the same as the gene of this variant.
                    # print(HEADER_dict) 
                    # {'Allele': 0, 'Consequence': 1, 'IMPACT': 2, 'SYMBOL': 3, 'Gene': 4, 'Feature_type': 5, 'Feature': 6, 
                    # #'BIOTYPE': 7, 'EXON': 8, 'INTRON': 9, 'HGVSc': 10, 'HGVSp': 11, 'cDNA_position': 12, 'CDS_position': 13, 
                    # 'Protein_position': 14, 'Amino_acids': 15, 'Codons': 16, 'Existing_variation': 17, 'REF_ALLELE': 18, 
                    # 'DISTANCE': 19, 'STRAND': 20, 'FLAGS': 21, 'PICK': 22, 'SYMBOL_SOURCE': 23, 'HGNC_ID': 24, 'CANONICAL': 25, 
                    # 'REFSEQ_MATCH': 26, 'REFSEQ_OFFSET': 27, 'GIVEN_REF': 28, 'USED_REF': 29, 'BAM_EDIT': 30, 'HERITANCE': 31}
                    #change the HERITANCE value from HEADER_dict to the value of gene2disease_dict from the Gene value
                    #print(format_info_list)
                    #['T', 'downstream_gene_variant', 'MODIFIER', 'ARSE', '415', 'Transcript', 'NM_000047.3', 'protein_coding',
                    #  '', '', '', '', '', '', '', '', '', '', 'C', '4991', '-1', '', '', 'EntrezGene', '', '', '', '', 'C', 'C', '']
                    
                    format_info_list.insert(HEADER_dict["HERITANCE"],gene2disease_dict[Gene])
                    n=1
                    break
                else:
                    continue #in case they are not the same, go to the next variant.

        if n == 1:
            new_line = format_info_list[HEADER_dict["SYMBOL"]] + ";" + format_info_list[HEADER_dict["HERITANCE"]]
            line[7]=new_line
            print("\t".join(line))

print(format_info_list)

