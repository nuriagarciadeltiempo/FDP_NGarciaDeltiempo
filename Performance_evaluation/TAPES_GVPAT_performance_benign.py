import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.metrics import roc_curve, auc


#my results obtained with intern annotation
myresults_ClinGer = pd.read_csv("/Users/nuriagarciadeltiempo/Desktop/results/TAPES_validated_hg19.txt", sep="\t")
myresults_Zhang = pd.read_csv("/Users/nuriagarciadeltiempo/Desktop/results2/TAPES_benchmark_hg19.txt", sep="\t")

myresults_ClinGer_GVPAT = pd.read_csv("/Users/nuriagarciadeltiempo/Desktop/GVPAT/GVPAT_Benchmark_vep.txt", sep="\t")
myresults_Zhang_GVPAT = pd.read_csv("/Users/nuriagarciadeltiempo/Desktop/GVPAT/GVPAT_Validated_vep.txt", sep="\t")


#paper results
paper_results_ClinGer = pd.read_csv("/Users/nuriagarciadeltiempo/Desktop/benchmark_data/Validation/TAPES_Validation.txt", sep="\t")
paper_results_Zhang = pd.read_csv("/Users/nuriagarciadeltiempo/Desktop/benchmark_data/Initial_Benchmark/TAPES_Benchmark_try.csv", sep = ",")

paper_results_ClinGer_GVPAT = pd.read_csv("/Users/nuriagarciadeltiempo/Desktop/GVPAT/GVPAT_Benchmark_annovar.txt", sep="\t")
paper_results_Zhang_GVPAT = pd.read_csv("/Users/nuriagarciadeltiempo/Desktop/Zhang_results_GVPAT.txt", sep="\t")


# reference data used in the paper as control and reference data. This reference data is the one used as control for the entire analysis
reference_Zhang = pd.read_csv("/Users/nuriagarciadeltiempo/Desktop/S2_Table_Zhang.csv", sep="\t")
#print(reference_Zhang.iloc[:, -1])
reference_ClinGer = pd.read_csv("/Users/nuriagarciadeltiempo/Desktop/S3_Table_ClinGer.txt", sep="\t")



####### First step: Use the columns that are needed for each of the imported files #######

myresults_ClinGer = myresults_ClinGer.iloc[:, [0, 1, 3, 4, -1]]
myresults_ClinGer.columns = ["Chr", "Start", "Ref", "Alt", "Prediction_ACMG_Tapes"]
myresults_ClinGer_filtered = myresults_ClinGer.drop_duplicates(keep="last")
print(myresults_ClinGer_filtered.shape)
myresults_Zhang = myresults_Zhang.iloc[:, [0, 1, 3, 4, -1]]
myresults_Zhang.columns = ["Chr", "Start", "Ref", "Alt", "Prediction_ACMG_Tapes"]
myresults_Zhang_filtered =myresults_Zhang.drop_duplicates(keep="last")
print(myresults_Zhang_filtered.shape)

paper_results_ClinGer = paper_results_ClinGer.iloc[:, [0, 1, 3, 4, -1]]
paper_results_ClinGer.columns = ["Chr", "Start", "Ref", "Alt", "Prediction_ACMG_Tapes"]
paper_results_ClinGer_filtered = paper_results_ClinGer.drop_duplicates(keep="last")
print(paper_results_ClinGer_filtered.shape)
paper_results_Zhang = paper_results_Zhang.iloc[:, [0, 1, 3, 4, -1]]
paper_results_Zhang.columns = ["Chr", "Start", "Ref", "Alt", "Prediction_ACMG_Tapes"]
paper_results_Zhang_filtered = paper_results_Zhang.drop_duplicates(keep= "last")
print(paper_results_Zhang_filtered.shape)


reference_Zhang = reference_Zhang.iloc[:, [ 0, 1, 3, 4, -1]]
reference_Zhang.columns = ["Chr", "Start", "Ref", "Alt", "Panel_Decision"]
reference_Zhang_filtered = reference_Zhang.drop_duplicates(keep=    "last")
print(reference_Zhang_filtered.shape)
reference_ClinGer = reference_ClinGer.iloc[:, [0, 1, 3, 4, 13]]
reference_ClinGer.columns = ["Chr", "Start", "Ref", "Alt", "Panel_Decision"]
reference_ClinGer_filtered = reference_ClinGer.drop_duplicates(keep="last")
print(reference_ClinGer_filtered.shape)

myresults_ClinGer_GVPAT = myresults_ClinGer_GVPAT.iloc[:, [0, 1, 3, 4, -1]]
myresults_ClinGer_GVPAT.columns = ["Chr", "Start", "Ref", "Alt", "Prediction_ACMG_GVPAT"]
myresults_ClinGer_GVPAT_filtered = myresults_ClinGer_GVPAT.drop_duplicates(keep="last")
print(myresults_ClinGer_GVPAT_filtered.shape)
myresults_Zhang_GVPAT = myresults_Zhang_GVPAT.iloc[:, [0, 1, 3, 4, -1]]
myresults_Zhang_GVPAT.columns = ["Chr", "Start", "Ref", "Alt", "Prediction_ACMG_GVPAT"]
myresults_Zhang_GVPAT_filtered = myresults_Zhang_GVPAT.drop_duplicates(keep="last")
print(myresults_Zhang_GVPAT_filtered.shape)

paper_results_ClinGer_GVPAT = paper_results_ClinGer_GVPAT.iloc[:, [0, 1, 3, 4, -1]]
paper_results_ClinGer_GVPAT.columns = ["Chr", "Start", "Ref", "Alt", "Prediction_ACMG_GVPAT"]
paper_results_ClinGer_GVPAT_filtered = paper_results_ClinGer_GVPAT.drop_duplicates(keep="last")
print(paper_results_ClinGer_GVPAT_filtered.shape)

paper_results_Zhang_GVPAT = paper_results_Zhang_GVPAT.iloc[:, [0, 1, 3, 4, -1]]
paper_results_Zhang_GVPAT.columns = ["Chr", "Start", "Ref", "Alt", "Prediction_ACMG_GVPAT"]
paper_results_Zhang_GVPAT_filtered = paper_results_Zhang_GVPAT.drop_duplicates(keep="last")
print(paper_results_Zhang_GVPAT_filtered.shape)



## change the nomenclature of pathogenic categorization in the reference data

#remove elements that are VUS OR U or Uncertain significance
reference_ClinGer_filtered = reference_ClinGer_filtered[reference_ClinGer_filtered.iloc[:, -1] != "VUS"]
reference_ClinGer_filtered = reference_ClinGer_filtered[reference_ClinGer_filtered.iloc[:, -1] != "Uncertain Significance"]
reference_ClinGer_filtered = reference_ClinGer_filtered[reference_ClinGer_filtered.iloc[:, -1] != "U"]

reference_Zhang_filtered = reference_Zhang_filtered[reference_Zhang_filtered.iloc[:, -1] != "VUS"]
reference_Zhang_filtered = reference_Zhang_filtered[reference_Zhang_filtered.iloc[:, -1] != "Uncertain Significance"]
reference_Zhang_filtered = reference_Zhang_filtered[reference_Zhang_filtered.iloc[:, -1] != "U"]


#reference_copy_Zhang= reference_Zhang.copy()

def change_name(x):
    if (x == "P" or x== "PP" or x =="Likely Pathogenic" or x== "Pathogenic"):
        return "Pathogenic"
    if ( x=="PB" or x=="BA" or x=="B" or x == "Benign" or x =="Benign auto" or x == "Likely Benign"):
        return "Benign"
    if (x == "VUS"):
        return "VUS"

reference_Zhang_filtered['Panel_Decision'] = reference_Zhang_filtered['Panel_Decision'].apply(change_name)
reference_ClinGer_filtered['Panel_Decision'] = reference_ClinGer_filtered['Panel_Decision'].apply(change_name)
myresults_ClinGer_filtered['Prediction_ACMG_Tapes'] = myresults_ClinGer_filtered['Prediction_ACMG_Tapes'].apply(change_name)
myresults_Zhang_filtered['Prediction_ACMG_Tapes'] = myresults_Zhang_filtered['Prediction_ACMG_Tapes'].apply(change_name)
paper_results_ClinGer_filtered['Prediction_ACMG_Tapes'] = paper_results_ClinGer_filtered['Prediction_ACMG_Tapes'].apply(change_name)
paper_results_Zhang_filtered['Prediction_ACMG_Tapes'] = paper_results_Zhang_filtered['Prediction_ACMG_Tapes'].apply(change_name)

myresults_ClinGer_GVPAT_filtered['Prediction_ACMG_GVPAT'] = myresults_ClinGer_GVPAT_filtered['Prediction_ACMG_GVPAT'].apply(change_name)
myresults_Zhang_GVPAT_filtered['Prediction_ACMG_GVPAT'] = myresults_Zhang_GVPAT_filtered['Prediction_ACMG_GVPAT'].apply(change_name)
paper_results_ClinGer_GVPAT_filtered['Prediction_ACMG_GVPAT'] = paper_results_ClinGer_GVPAT_filtered['Prediction_ACMG_GVPAT'].apply(change_name)
paper_results_Zhang_GVPAT_filtered['Prediction_ACMG_GVPAT'] = paper_results_Zhang_GVPAT_filtered['Prediction_ACMG_GVPAT'].apply(change_name)

print('shape zhang',reference_Zhang_filtered.shape)
print('shape clin',reference_ClinGer_filtered.shape)

#print( paper_results_Zhang_GVPAT_filtered)
#reference_Zhang["Panel_Decision"].value_counts()
#reference_ClinGer["Panel_Decision"].value_counts()
#myresults_ClinGer["Prediction_ACMG_Tapes"].value_counts()

####### 3rd step: save the variants that have in common myresults_ClinGer with referenece_InterVar #######
myresults_ClinGer_merged = myresults_ClinGer_filtered.merge(reference_ClinGer_filtered, on=["Chr", "Start", "Ref", "Alt"], how="inner")
#print(myresults_ClinGer_merged.shape)

# save the variants that have in common myresults_Zhang with reference_Zhang
myresults_Zhang_merged  = myresults_Zhang_filtered.merge(reference_Zhang_filtered, on=["Chr", "Start", "Ref", "Alt"], how="inner")
print('Zhang annovar filtered tapes',myresults_Zhang_filtered.shape)
print('Zhang vep merged tapes',myresults_Zhang_merged.shape)

## save the variants that have in common paper_results_ClinGer with reference_ClinGer
paper_results_ClinGer_merged = paper_results_ClinGer_filtered.merge(reference_ClinGer_filtered, on=["Chr", "Start", "Ref", "Alt"], how="inner")
#print(paper_results_ClinGer_merged.shape)

# save the variants that have in common paper_results_Zhang with reference_Zhang
paper_results_Zhang_merged = paper_results_Zhang_filtered.merge(reference_Zhang_filtered, on=["Chr", "Start", "Ref", "Alt"], how="inner")
print('Zhang annovar filtered tapes',paper_results_Zhang_filtered.shape)
print('Zhang annovar merged tapes',paper_results_Zhang_merged.shape)


myresults_ClinGer_GVPAT_merged = myresults_ClinGer_GVPAT_filtered.merge(reference_ClinGer_filtered, on=["Chr", "Start", "Ref", "Alt"], how="inner")
#print(myresults_ClinGer_GVPAT_merged.shape)
myresults_Zhang_GVPAT_merged  = myresults_Zhang_GVPAT_filtered.merge(reference_Zhang_filtered, on=["Chr", "Start", "Ref", "Alt"], how="inner")
print('zhang annovar filtered gvpat',myresults_Zhang_GVPAT_filtered.shape)
print('zhang vep merged gvpat',myresults_Zhang_GVPAT_merged.shape)
paper_results_ClinGer_GVPAT_merged = paper_results_ClinGer_GVPAT_filtered.merge(reference_ClinGer_filtered, on=["Chr", "Start", "Ref", "Alt"], how="inner")
#print(paper_results_ClinGer_GVPAT_merged.shape)
paper_results_Zhang_GVPAT_merged = paper_results_Zhang_GVPAT_filtered.merge(reference_Zhang_filtered, on=["Chr", "Start", "Ref", "Alt"], how="inner")
print('zhang annovar filtered gvpat',paper_results_Zhang_GVPAT_filtered.shape)
print('zhang annovar merged gvpat',paper_results_Zhang_GVPAT_merged.shape)




####### 4th step: Create two groups, benign and pathogenic, for each of the files (Likely benign and
#  patho goes respectively on their groups) #######

#save the files 
myresults_ClinGer.to_csv("/Users/nuriagarciadeltiempo/Desktop/try/TAPES_validated_hg19_string.txt", sep="\t", index=False)
myresults_Zhang.to_csv("/Users/nuriagarciadeltiempo/Desktop/try/TAPES_benchmark_hg19_string.txt", sep="\t", index=False)
paper_results_ClinGer.to_csv("/Users/nuriagarciadeltiempo/Desktop/try/TAPES_Validation_string.txt", sep="\t", index=False)
paper_results_Zhang.to_csv("/Users/nuriagarciadeltiempo/Desktop/try/TAPES_Benchmark_try_string.txt", sep="\t", index=False)

##########################################
##### the data has to be {0, 1} for the precision-recall curve to work #############
###########################################

myresults_ClinGer_binary_b = myresults_ClinGer_merged.copy()
myresults_Zhang_binary_b = myresults_Zhang_merged.copy()
paper_results_ClinGer_binary_b = paper_results_ClinGer_merged.copy()
paper_results_Zhang_binary_b = paper_results_Zhang_merged.copy()

myresults_ClinGer_binary_b_GVPAT = myresults_ClinGer_GVPAT_merged.copy()
myresults_Zhang_binary_b_GVPAT = myresults_Zhang_GVPAT_merged.copy()
paper_results_ClinGer_binary_b_GVPAT = paper_results_ClinGer_GVPAT_merged.copy()
paper_results_Zhang_binary_b_GVPAT = paper_results_Zhang_GVPAT_merged.copy()


def categ_beningnicity(x):
    if (x == "Pathogenic" or x== "VUS"):
        return int(0)
    if (x=="Benign"):
        return int(1)

def categ_df_benign(df_binary):
    if "Prediction_ACMG_Tapes" in df_binary.columns:
        df_binary['Prediction_ACMG_Tapes'] = df_binary['Prediction_ACMG_Tapes'].apply(categ_beningnicity)
    elif "Prediction_ACMG_GVPAT" in df_binary.columns:
        df_binary['Prediction_ACMG_GVPAT'] = df_binary['Prediction_ACMG_GVPAT'].apply(categ_beningnicity)
    df_binary['Panel_Decision'] = df_binary['Panel_Decision'].apply(categ_beningnicity)
    return(df_binary)


myresults_ClinGer_binary_b = categ_df_benign(myresults_ClinGer_binary_b)
myresults_Zhang_binary_b = categ_df_benign(myresults_Zhang_binary_b)
paper_results_ClinGer_binary_b = categ_df_benign(paper_results_ClinGer_binary_b)
paper_results_Zhang_binary_b = categ_df_benign(paper_results_Zhang_binary_b)

myresults_ClinGer_binary_b_GVPAT = categ_df_benign(myresults_ClinGer_binary_b_GVPAT)
myresults_Zhang_binary_b_GVPAT = categ_df_benign(myresults_Zhang_binary_b_GVPAT)
paper_results_ClinGer_binary_b_GVPAT = categ_df_benign(paper_results_ClinGer_binary_b_GVPAT)
paper_results_Zhang_binary_b_GVPAT = categ_df_benign(paper_results_Zhang_binary_b_GVPAT)



##### 5th step: Create the confusion matrices for each of the files #######

print("These are the confusion matrices with the paper data:   ")
print('Zhang: ')
print(metrics.confusion_matrix(paper_results_Zhang_binary_b.iloc[:, -1], paper_results_Zhang_binary_b.iloc[:, -2]))
print("\n")
print('ClinGer: ')
print(metrics.confusion_matrix(paper_results_ClinGer_binary_b.iloc[:, -1], paper_results_ClinGer_binary_b.iloc[:, -2]))
print("\n")

print("These are the confusion matrices with the internal data:   ")
print('Zhang: ')
print(metrics.confusion_matrix(myresults_Zhang_binary_b.iloc[:, -1], myresults_Zhang_binary_b.iloc[:, -2]))
print("\n")
print('ClinGer: ')
print(metrics.confusion_matrix(myresults_ClinGer_binary_b.iloc[:, -1], myresults_ClinGer_binary_b.iloc[:, -2]))
print("\n")



print("These are the confusion matrices with the paper data for GVPAT:   ")
print('Zhang: ')
print(metrics.confusion_matrix(paper_results_Zhang_binary_b_GVPAT.iloc[:, -1], paper_results_Zhang_binary_b_GVPAT.iloc[:, -2]))
print("\n")
print('ClinGer: ')
print(metrics.confusion_matrix(paper_results_ClinGer_binary_b_GVPAT.iloc[:, -1], paper_results_ClinGer_binary_b_GVPAT.iloc[:, -2]))
print("\n")

print("These are the confusion matrices with the internal data for GVPAT:   ")
print('Zhang: ')
print(metrics.confusion_matrix(myresults_Zhang_binary_b_GVPAT.iloc[:, -1], myresults_Zhang_binary_b_GVPAT.iloc[:, -2]))
print("\n")
print('ClinGer: ')
print(metrics.confusion_matrix(myresults_ClinGer_binary_b_GVPAT.iloc[:, -1], myresults_ClinGer_binary_b_GVPAT.iloc[:, -2]))
print("\n")

######################################################
##### 6th step: Create the classification reports #####
####### Include precision, recall, f1-score, support #######


print("These are the classification reports for the paper data:   ")
print('Zhang: ')
print(metrics.classification_report(paper_results_Zhang_binary_b.iloc[:, -1], paper_results_Zhang_binary_b.iloc[:, -2], digits=2))
print("\n")

print('ClinGer: ')
print(metrics.classification_report(paper_results_ClinGer_binary_b.iloc[:, -1], paper_results_ClinGer_binary_b.iloc[:, -2], digits=2))
print("\n")

print("These are the classification reports for the internal data:   ")
print('Zhang: ')
print(metrics.classification_report(myresults_Zhang_binary_b.iloc[:, -1], myresults_Zhang_binary_b.iloc[:, -2], digits=2))
print("\n")

print('ClinGer: ')
print(metrics.classification_report(myresults_ClinGer_binary_b.iloc[:, -1], myresults_ClinGer_binary_b.iloc[:, -2], digits=2))
print("\n")



print("These are the classification reports for the paper data for GVPAT:   ")
print('Zhang: ')
print(metrics.classification_report(paper_results_Zhang_binary_b_GVPAT.iloc[:, -1], paper_results_Zhang_binary_b_GVPAT.iloc[:, -2], digits=2))
print("\n")

print('ClinGer: ')
print(metrics.classification_report(paper_results_ClinGer_binary_b_GVPAT.iloc[:, -1], paper_results_ClinGer_binary_b_GVPAT.iloc[:, -2], digits=2))
print("\n")

print("These are the classification reports for the internal data for GVPAT:   ")
print('Zhang: ')
print(metrics.classification_report(myresults_Zhang_binary_b_GVPAT.iloc[:, -1], myresults_Zhang_binary_b_GVPAT.iloc[:, -2], digits=2))
print("\n")

print('ClinGer: ')
print(metrics.classification_report(myresults_ClinGer_binary_b_GVPAT.iloc[:, -1], myresults_ClinGer_binary_b_GVPAT.iloc[:, -2], digits=2))
print("\n")

####### 7th step: Create Precision-Recall curves (AUC) ######


myresults_Zhang_1 = myresults_Zhang_binary_b["Panel_Decision"].astype(int).to_list()
myresults_Zhang_2 = myresults_Zhang_binary_b["Prediction_ACMG_Tapes"].astype(int).to_list()
#print(myresults_Zhang_2)
paper_results_Zhang_1 = paper_results_Zhang_binary_b["Panel_Decision"].astype(int).to_list()
paper_results_Zhang_2 = paper_results_Zhang_binary_b["Prediction_ACMG_Tapes"].astype(int).to_list()

myresults_ClinGer_1 = myresults_ClinGer_binary_b["Panel_Decision"].astype(int).to_list()
myresults_ClinGer_2 = myresults_ClinGer_binary_b["Prediction_ACMG_Tapes"].astype(int).to_list()

paper_results_ClinGer_1 = paper_results_ClinGer_binary_b["Panel_Decision"].astype(int).to_list()
paper_results_ClinGer_2 = paper_results_ClinGer_binary_b["Prediction_ACMG_Tapes"].astype(int).to_list()

myresults_ClinGer_1_GVPAT = myresults_ClinGer_binary_b_GVPAT["Panel_Decision"].astype(int).to_list()
myresults_ClinGer_2_GVPAT = myresults_ClinGer_binary_b_GVPAT["Prediction_ACMG_GVPAT"].astype(int).to_list()

myresults_Zhang_1_GVPAT = myresults_Zhang_binary_b_GVPAT["Panel_Decision"].astype(int).to_list()
myresults_Zhang_2_GVPAT = myresults_Zhang_binary_b_GVPAT["Prediction_ACMG_GVPAT"].astype(int).to_list()

paper_results_ClinGer_1_GVPAT = paper_results_ClinGer_binary_b_GVPAT["Panel_Decision"].astype(int).to_list()
paper_results_ClinGer_2_GVPAT = paper_results_ClinGer_binary_b_GVPAT["Prediction_ACMG_GVPAT"].astype(int).to_list()

paper_results_Zhang_1_GVPAT = paper_results_Zhang_binary_b_GVPAT["Panel_Decision"].astype(int).to_list()
paper_results_Zhang_2_GVPAT = paper_results_Zhang_binary_b_GVPAT["Prediction_ACMG_GVPAT"].astype(int).to_list()



##### Computing the probabilities and performin logistic regression for the precision-recall curve #######
#print(type(myresults_Zhang.loc[:, "Prediction_ACMG_tapes"]))   

#precision recall curve for my results prediction of benignicity

precision_exp2, recall_exp2, thresholds_exp2 = metrics.precision_recall_curve(myresults_Zhang_1, myresults_Zhang_2)
precision_exp1, recall_exp1, thresholds_exp1 = metrics.precision_recall_curve(paper_results_Zhang_1, paper_results_Zhang_2)
precision_exp2_gvp, recall_exp2_gvp, thresholds_exp2_gvp = metrics.precision_recall_curve(myresults_Zhang_1_GVPAT, myresults_Zhang_2_GVPAT)
precision_exp1_gvp, recall_exp1_gvp, thresholds_exp1_gvp = metrics.precision_recall_curve(paper_results_Zhang_1_GVPAT, paper_results_Zhang_2_GVPAT)
plt.plot(recall_exp1, precision_exp1, color = 'lightblue', label = ' TAPES with ANNOVAR annotation') # plot the precision-recall curve
plt.plot(recall_exp2, precision_exp2, color = 'orange',label = 'TAPES with VEP annotation') # plot the precision-recall curve
plt.plot(recall_exp1_gvp, precision_exp1_gvp, color = 'blue', label = ' GVPAT with ANNOVAR annotation') # plot the precision-recall curve
plt.plot(recall_exp2_gvp, precision_exp2_gvp, color = 'red',label = 'GVPAT with VEP annotation') # plot the precision-recall curve


plt.legend()
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve Zhang data for Benign prediction')
plt.show()

#precision recall curve for paper results prediction of benignicity

precision_inter2, recall_inter2, thresholds_inter2 = metrics.precision_recall_curve(myresults_ClinGer_2, myresults_ClinGer_1)
precision_inter1, recall_inter1, thresholds_inter1 = metrics.precision_recall_curve(paper_results_ClinGer_2, paper_results_ClinGer_1)
precision_inter2_gvp, recall_inter2_gvp, thresholds_inter2_gvp = metrics.precision_recall_curve(myresults_ClinGer_2_GVPAT, myresults_ClinGer_1_GVPAT)
precision_inter1_gvp, recall_inter1_gvp, thresholds_inter1_gvp = metrics.precision_recall_curve(paper_results_ClinGer_2_GVPAT, paper_results_ClinGer_1_GVPAT)
plt.plot(recall_inter1, precision_inter1, color = 'lightblue', label = 'TAPES with ANNOVAR annotation') # plot the precision-recall curve
plt.plot(recall_inter2, precision_inter2, color = 'orange',label = 'TAPES with VEP annotation') # plot the precision-recall curve
plt.plot(recall_inter1_gvp, precision_inter1_gvp, color = 'blue', label = 'GVPAT with ANNOVAR annotation') # plot the precision-recall curve
plt.plot(recall_inter2_gvp, precision_inter2_gvp, color = 'red',label = 'GVPAT with VEP annotation') # plot the precision-recall curve

plt.legend()
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve ClinGer data for Benign prediction')
plt.show()


####### Sixth step: Create the ROC curve for each of the performance analysis #######

fpr_ExP1, tpr_ExP1, threshold_ExP1 = roc_curve(paper_results_Zhang_2, paper_results_Zhang_1)
fpr_ExP1_gvp, tpr_ExP1_gvp, threshold_ExP1_gvp = roc_curve(paper_results_Zhang_2_GVPAT, paper_results_Zhang_1_GVPAT)
fpr_ExP2, tpr_ExP2, threshold_ExP2 = roc_curve(myresults_Zhang_2, myresults_Zhang_1)
fpr_ExP2_gvp, tpr_ExP2_gvp, threshold_ExP2_gvp = roc_curve(myresults_Zhang_2_GVPAT, myresults_Zhang_1_GVPAT)


# Calcular el área bajo la curva ROC (AUC-ROC)
auc_score1 = auc(fpr_ExP1, tpr_ExP1)
auc_score1_gvp = auc(fpr_ExP1_gvp, tpr_ExP1_gvp)
auc_score2 = auc(fpr_ExP2, tpr_ExP2)
auc_score2_gvp = auc(fpr_ExP2_gvp, tpr_ExP2_gvp)

# Graficar la curva ROC
plt.plot(fpr_ExP1, tpr_ExP1, label='ROC curve TAPES with ANNOVAR annotation(AUC = %0.2f)' % auc_score1)
plt.plot(fpr_ExP2, tpr_ExP2, label='ROC curve TAPES with VEP annotation (AUC = %0.2f)' % auc_score2)
plt.plot(fpr_ExP1_gvp, tpr_ExP1_gvp, label='ROC curve GVPAT with ANNOVAR annotation(AUC = %0.2f)' % auc_score1_gvp)
plt.plot(fpr_ExP2_gvp, tpr_ExP2_gvp, label='ROC curve GVPAT with VEP annotation (AUC = %0.2f)' % auc_score2_gvp)
plt.plot([0, 1], [0, 1], 'r--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('FPS')
plt.ylabel('TPR')
plt.title('ROC-curve for Benignicity Prediction and Zhang data')
plt.legend(loc="lower right")
plt.show()



fpr_inter1, tpr_inter1, threshold_inter1 = roc_curve(paper_results_ClinGer_2, paper_results_ClinGer_1)
fpr_inter1_gvp, tpr_inter1_gvp, threshold_inter1_gvp = roc_curve(paper_results_ClinGer_2_GVPAT, paper_results_ClinGer_1_GVPAT)
fpr_inter2, tpr_inter2, threshold_inter2 = roc_curve(myresults_ClinGer_2, myresults_ClinGer_1)
fpr_inter2_gvp, tpr_inter2_gvp, threshold_inter2_gvp = roc_curve(myresults_ClinGer_2_GVPAT, myresults_ClinGer_1_GVPAT)


# Calcular el área bajo la curva ROC (AUC-ROC)
auc_score1 = auc(fpr_inter1, tpr_inter1)
auc_score1_gvp = auc(fpr_inter1_gvp, tpr_inter1_gvp)
auc_score2 = auc(fpr_inter2, tpr_inter2)
auc_score2_gvp = auc(fpr_inter2_gvp, tpr_inter2_gvp)

# Graficar la curva ROC
plt.plot(fpr_inter1, tpr_inter1, label='ROC curve TAPES with ANNOVAR annotation (AUC = %0.2f)' % auc_score1)
plt.plot(fpr_inter2, tpr_inter2, label='ROC curve TAPES with VEP annotation(AUC = %0.2f)' % auc_score2)
plt.plot(fpr_inter1_gvp, tpr_inter1_gvp, label='ROC curve GVPAT with ANNOVAR annotation(AUC = %0.2f)' % auc_score1_gvp)
plt.plot(fpr_inter2_gvp, tpr_inter2_gvp, label='ROC curve GVPAT with VEP annotation (AUC = %0.2f)' % auc_score2_gvp)
plt.plot([0, 1], [0, 1], 'r--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('FPS')
plt.ylabel('TPR')
plt.title('ROC-curve for Benignicity Prediction and ClinGer data')
plt.legend(loc="lower right")
plt.show()


