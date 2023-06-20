# Final Degree Project

During my Final Degree Project I had to create many small pieces of code to filter and merge relevant information from the data obtained from different databases such as OMIM, ClinVar, Genome Project etc...
This GitHub includes the scripts I created and used with the resulting outputs.

Furthermore, it contains the main pieces of code that I have been modifying in order to optimize and standardize the TAPES algorithm.


## GVPAT Functionalities

### Annotation
GVPAT can be used as an ANNOVAR wrapper for easy VCF/gzipped VCF/BCF annotation (this requires downloading ANNOVAR).
To annotate the VCF file with ANNOVAR, the command is:
`python3 gvpat.py annotate -i /path/to/variant_file.bcf -o /path/to/annotated_file.vcf --amcg`

### Variant Categorization
Prediction of pathogenic variants is made once the VCF has been annotated. The command used to make the variant priorization is:

`python3 gvpat.py sort -i /path/to/annotated_file.vcf -o /path/to/output/ --tab`
