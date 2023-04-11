
import pandas as pd

OMIM_file = pd.read_excel(r'/Users/nuriagarciadeltiempo/Desktop/OMIM_by_column.xls') # import the file
adult_onset = pd.read_excel(r'/Users/nuriagarciadeltiempo/Desktop/BS2_rec_dom_ad.xls', header=None) # import the file

# Create subset of the dataframe where the last column contains X-linked
df = OMIM_file[OMIM_file.iloc[:, -1].str.contains('X-linked')]

#save the new dataframe as a new excel file
df.to_excel(r'/Users/nuriagarciadeltiempo/Desktop/OMIM_X-linked.xls', index = False, header=None)

#remove the rows from adult_onset  that contains the column adult_onset == 1.
df_adult_onset = adult_onset[adult_onset.iloc[:, -1] != 1]

# save df_adult_onset as a new excel file
df_adult_onset.to_excel(r'/Users/nuriagarciadeltiempo/Desktop/BS2_rec_dom_ad_X-linked.xls', index = False, header=None)

# create a new dataframe that contains all the genes from the OMIM_X-linked.xls file that are also in the BS2_rec_dom_ad_X-linked.xls file
final_df = df[df.iloc[:, 0].isin(df_adult_onset[0])]
# save the final_df as a new excel file
final_df.to_excel(r'/Users/nuriagarciadeltiempo/Desktop/OMIM_X-linked_non_adult_onset.xls', index = False, header=None)








