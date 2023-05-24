import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.metrics import roc_curve, auc


#my results obtained with intern annotation
myresults_InterVar = pd.read_csv("/Users/nuriagarciadeltiempo/Desktop/results/TAPES_validated_hg19.txt", sep="\t")
myresults_ExpertPanel = pd.read_csv("/Users/nuriagarciadeltiempo/Desktop/results2/TAPES_benchmark_hg19.txt", sep="\t")

#paper results
paper_results_InterVar = pd.read_csv("/Users/nuriagarciadeltiempo/Desktop/benchmark_data/Validation/TAPES_Validation.txt", sep="\t")
paper_results_ExpertPanel = pd.read_excel("/Users/nuriagarciadeltiempo/Desktop/benchmark_data/Initial_Benchmark/TAPES_Benchmark.xlsx")

# reference data used in the paper as control and reference data. This reference data is the one used as control for the entire analysis
reference_ExpertPanel = pd.read_excel("/Users/nuriagarciadeltiempo/Desktop/S2_Table_ExpertPanel.xlsx")
#print(reference_ExpertPanel.iloc[:, -1])
reference_InterVar = pd.read_csv("/Users/nuriagarciadeltiempo/Desktop/S3_Table_InterVar.txt", sep="\t")

####### First step: Use the columns that are needed for each of the imported files #######

myresults_InterVar = myresults_InterVar.iloc[:, [0, 1, 3, 4, -1]]
myresults_ExpertPanel = myresults_ExpertPanel.iloc[:, [0, 1, 3, 4, -1]]


paper_results_InterVar = paper_results_InterVar.iloc[:, [0, 1, 3, 4, -1]]
paper_results_ExpertPanel = paper_results_ExpertPanel.iloc[:, [0, 1, 3, 4, -1]]

reference_ExpertPanel = reference_ExpertPanel.iloc[:, [ 0, 1, 3, 4, -1]]


reference_InterVar = reference_InterVar.iloc[:, [0, 1, 3, 4, 13]]
## change the nomenclature of pathogenic categorization in the reference data

#remove elements that are VUS OR U or Uncertain significance
reference_InterVar = reference_InterVar[reference_InterVar.iloc[:, -1] != "VUS"]
reference_InterVar = reference_InterVar[reference_InterVar.iloc[:, -1] != "Uncertain Significance"]
reference_InterVar = reference_InterVar[reference_InterVar.iloc[:, -1] != "U"]

reference_ExpertPanel = reference_ExpertPanel[reference_ExpertPanel.iloc[:, -1] != "VUS"]
reference_ExpertPanel = reference_ExpertPanel[reference_ExpertPanel.iloc[:, -1] != "Uncertain Significance"]
reference_ExpertPanel = reference_ExpertPanel[reference_ExpertPanel.iloc[:, -1] != "U"]

myresults_InterVar = myresults_InterVar[myresults_InterVar.iloc[:, -1] != "VUS"]
myresults_InterVar = myresults_InterVar[myresults_InterVar.iloc[:, -1] != "Uncertain Significance"]
myresults_InterVar = myresults_InterVar[myresults_InterVar.iloc[:, -1] != "U"]

myresults_ExpertPanel = myresults_ExpertPanel[myresults_ExpertPanel.iloc[:, -1] != "VUS"]
myresults_ExpertPanel = myresults_ExpertPanel[myresults_ExpertPanel.iloc[:, -1] != "Uncertain Significance"]
myresults_ExpertPanel = myresults_ExpertPanel[myresults_ExpertPanel.iloc[:, -1] != "U"]

paper_results_InterVar = paper_results_InterVar[paper_results_InterVar.iloc[:, -1] != "VUS"]
paper_results_InterVar = paper_results_InterVar[paper_results_InterVar.iloc[:, -1] != "Uncertain Significance"]
paper_results_InterVar = paper_results_InterVar[paper_results_InterVar.iloc[:, -1] != "U"]

paper_results_ExpertPanel = paper_results_ExpertPanel[paper_results_ExpertPanel.iloc[:, -1] != "VUS"]
paper_results_ExpertPanel = paper_results_ExpertPanel[paper_results_ExpertPanel.iloc[:, -1] != "Uncertain Significance"]
paper_results_ExpertPanel = paper_results_ExpertPanel[paper_results_ExpertPanel.iloc[:, -1] != "U"]

list = []
for element in reference_ExpertPanel.iloc[:, -1]:
    list.append(element)


for element in reference_InterVar.iloc[:, -1]:
    if element == "Likely benign":
        reference_InterVar.iloc[:, -1] = reference_InterVar.iloc[:, -1].replace(element, "Benign")
    elif element == "Benign auto":
        reference_InterVar.iloc[:, -1] = reference_InterVar.iloc[:, -1].replace(element, "Benign")
    elif element == "Benign":
        reference_InterVar.iloc[:, -1] = reference_InterVar.iloc[:, -1].replace(element, "Benign")
    else:
        reference_InterVar.iloc[:, -1] = reference_InterVar.iloc[:, -1].replace(element, "Pathogenic")


list_random = []
patho_count = 0
benign_count = 0
for item in list:
    if item == "PP":
        item = "Pathogenic"
        list_random.append(item)
        patho_count += 1
    elif item == "P":
        item = "Pathogenic"
        list_random.append(item)
        patho_count += 1
    else:
        item = "Benign" 
        list_random.append(item)
        benign_count += 1
#print(list_random)
reference_ExpertPanel["Panel_Decision"] = list_random
reference_ExpertPanel.drop(reference_ExpertPanel.columns[[-2]], axis=1, inplace=True)
#print(reference_ExpertPanel)
#print(patho_count)
#print(benign_count)


####### Second step: save the variants that have in common myresults_InterVar with referenece_InterVar #######
myresults_InterVar = myresults_InterVar.merge(reference_InterVar, on=["Chr", "Start", "Ref", "Alt"], how="inner")
#print(myresults_InterVar)

# save the variants that have in common myresults_ExpertPanel with reference_ExpertPanel
myresults_ExpertPanel  = myresults_ExpertPanel.merge(reference_ExpertPanel, on=["Chr", "Start", "Ref", "Alt"], how="inner")
#print(myresults_ExpertPanel)

## save the variants that have in common paper_results_InterVar with reference_InterVar
paper_results_InterVar = paper_results_InterVar.merge(reference_InterVar, on=["Chr", "Start", "Ref", "Alt"], how="inner")
#print(paper_results_InterVar)

# save the variants that have in common paper_results_ExpertPanel with reference_ExpertPanel
paper_results_ExpertPanel = paper_results_ExpertPanel.merge(reference_ExpertPanel, on=["Chr", "Start", "Ref", "Alt"], how="inner")
#print(paper_results_ExpertPanel)



####### Third step: Create two groups, benign and pathogenic, for each of the files (Likely benign and
#  patho goes respectively on their groups) #######

for element in myresults_InterVar.iloc[:, -2]:
    if element == "Likely benign":
        myresults_InterVar.iloc[:, -2] = myresults_InterVar.iloc[:, -2].replace(element, "Benign")
    elif element == "Benign auto":
        myresults_InterVar.iloc[:, -2] = myresults_InterVar.iloc[:, -2].replace(element, "Benign")
    elif element == "Benign":
        myresults_InterVar.iloc[:, -2] = myresults_InterVar.iloc[:, -2].replace(element, "Benign")
    else:
        myresults_InterVar.iloc[:, -2] = myresults_InterVar.iloc[:, -2].replace(element, "Pathogenic")
#print(myresults_InterVar) #REVISED


for element in myresults_ExpertPanel.iloc[:, -2]:
    if element == "Pathogenic":
        myresults_ExpertPanel.iloc[:, -2] = myresults_ExpertPanel.iloc[:, -2].replace(element, "Pathogenic")

    elif element == "VUS":
        myresults_ExpertPanel.iloc[:, -2] = myresults_ExpertPanel.iloc[:, -2].replace(element, "Pathogenic")

    elif element == "Likely Pathogenic":
        myresults_ExpertPanel.iloc[:, -2] = myresults_ExpertPanel.iloc[:, -2].replace(element, "Pathogenic")

    else:
        myresults_ExpertPanel.iloc[:, -2] = myresults_ExpertPanel.iloc[:, -2].replace(element, "Benign")

#print(myresults_ExpertPanel) REVISED

for element in paper_results_InterVar.iloc[:, -2]:
    if element == "Likely benign":
        paper_results_InterVar.iloc[:, -2] = paper_results_InterVar.iloc[:, -2].replace(element, "Benign")
    elif element == "Benign auto":
        paper_results_InterVar.iloc[:, -2] = paper_results_InterVar.iloc[:, -2].replace(element, "Benign")
    elif element == "Benign":
        paper_results_InterVar.iloc[:, -2] = paper_results_InterVar.iloc[:, -2].replace(element, "Benign")
    else:
        paper_results_InterVar.iloc[:, -2] = paper_results_InterVar.iloc[:, -2].replace(element, "Pathogenic")
#print(paper_results_InterVar) REVISED

benign_count = 0
patho_count = 0
for element in paper_results_ExpertPanel.iloc[:, -2]:
    if element == "Pathogenic":
        paper_results_ExpertPanel.iloc[:, -2] = paper_results_ExpertPanel.iloc[:, -2].replace(element, "Pathogenic")
        benign_count += 1
    elif element == "Likely Pathogenic":
        paper_results_ExpertPanel.iloc[:, -2] = paper_results_ExpertPanel.iloc[:, -2].replace(element, "Pathogenic")
        benign_count += 1
    elif element == "VUS":
        paper_results_ExpertPanel.iloc[:, -2] = paper_results_ExpertPanel.iloc[:, -2].replace(element, "Pathogenic")
        benign_count += 1
    else:
        paper_results_ExpertPanel.iloc[:, -2] = paper_results_ExpertPanel.iloc[:, -2].replace(element, "Benign")
        patho_count += 1
#print(paper_results_ExpertPanel) REVISED


print(paper_results_ExpertPanel)
##### Fifth step: Create the confusion matrices for each of the files #######

print("These are the confusion matrices with the paper data:   ")
print('ExpertPanel: ')
print(metrics.confusion_matrix(paper_results_ExpertPanel.iloc[:, -1], paper_results_ExpertPanel.iloc[:, -2]))
print("\n")
print('InterVar: ')
print(metrics.confusion_matrix(paper_results_InterVar.iloc[:, -1], paper_results_InterVar.iloc[:, -2]))
print("\n")

print("These are the confusion matrices with the internal data:   ")
print('ExpertPanel: ')
print(metrics.confusion_matrix(myresults_ExpertPanel.iloc[:, -1], myresults_ExpertPanel.iloc[:, -2]))
print("\n")
print('InterVar: ')
print(metrics.confusion_matrix(myresults_InterVar.iloc[:, -1], myresults_InterVar.iloc[:, -2]))
print("\n")


print("These are the classification reports for the paper data:   ")
print('ExpertPanel: ')
print(metrics.classification_report(paper_results_ExpertPanel.iloc[:, -1], paper_results_ExpertPanel.iloc[:, -2], digits=2))
print("\n")
'''''
ExpertPanel: 
              precision    recall  f1-score   support

      Benign       0.99      0.98      0.98       250
  Pathogenic       0.98      0.99      0.99       289

    accuracy                           0.99       539
   macro avg       0.99      0.98      0.99       539
weighted avg       0.99      0.99      0.99       539
'''''

print('InterVar: ')
print(metrics.classification_report(paper_results_InterVar.iloc[:, -1], paper_results_InterVar.iloc[:, -2], digits=2))
print("\n")
'''''
InterVar: 
              precision    recall  f1-score   support

      Benign       0.83      0.60      0.70        40
  Pathogenic       0.93      0.98      0.95       216

    accuracy                           0.92       256
   macro avg       0.88      0.79      0.82       256
weighted avg       0.91      0.92      0.91       256
'''''

print("These are the classification reports for the internal data:   ")
print('ExpertPanel: ')
print(metrics.classification_report(myresults_ExpertPanel.iloc[:, -1], myresults_ExpertPanel.iloc[:, -2], digits=2))
print("\n")
''''
ExpertPanel: 
             precision    recall  f1-score   support

      Benign       0.99      0.99      0.99       189
  Pathogenic       1.00      1.00      1.00       249

    accuracy                           1.00       438
   macro avg       1.00      1.00      1.00       438
weighted avg       1.00      1.00      1.00       438
'''''

print('InterVar: ')
print(metrics.classification_report(myresults_InterVar.iloc[:, -1], myresults_InterVar.iloc[:, -2], digits=2))
print("\n")
'''
InterVar: 
              precision    recall  f1-score   support

      Benign       0.83      0.28      0.42        18
  Pathogenic       0.94      1.00      0.97       217

    accuracy                           0.94       235
   macro avg       0.89      0.64      0.69       235
weighted avg       0.93      0.94      0.93       235
'''



####### Fifth step: Create Precision-Recall curves (AUC) ######

## the data has to be {0, 1} for the precision-recall curve to work #############
for element in myresults_ExpertPanel.iloc[:, -2]:
    if element == "Benign":
        myresults_ExpertPanel.iloc[:, -2] = myresults_ExpertPanel.iloc[:, -2].replace(element, 0)
    else:
        myresults_ExpertPanel.iloc[:, -2] = myresults_ExpertPanel.iloc[:, -2].replace(element, 1)

for element in myresults_ExpertPanel.iloc[:, -1]:
    if element == "Benign":
        myresults_ExpertPanel.iloc[:, -1] = myresults_ExpertPanel.iloc[:, -1].replace(element, 0)
    else:
        myresults_ExpertPanel.iloc[:, -1] = myresults_ExpertPanel.iloc[:, -1].replace(element, 1)
#print(myresults_ExpertPanel)


for element in paper_results_ExpertPanel.iloc[:, -2]:
    if element == "Benign":
        paper_results_ExpertPanel.iloc[:, -2] = paper_results_ExpertPanel.iloc[:, -2].replace(element, 0)
    else:
        paper_results_ExpertPanel.iloc[:, -2] = paper_results_ExpertPanel.iloc[:, -2].replace(element, 1)

for element in paper_results_ExpertPanel.iloc[:, -1]:
    if element == "Benign":
        paper_results_ExpertPanel.iloc[:, -1] = paper_results_ExpertPanel.iloc[:, -1].replace(element, 0)
    else:
        paper_results_ExpertPanel.iloc[:, -1] = paper_results_ExpertPanel.iloc[:, -1].replace(element, 1)



for element in myresults_InterVar.iloc[:, -2]:
    if element == "Benign":
        myresults_InterVar.iloc[:, -2] = myresults_InterVar.iloc[:, -2].replace(element, 0)
    else:
        myresults_InterVar.iloc[:, -2] = myresults_InterVar.iloc[:, -2].replace(element, 1)

for element in myresults_InterVar.iloc[:, -1]:
    if element == "Benign":
        myresults_InterVar.iloc[:, -1] = myresults_InterVar.iloc[:, -1].replace(element, 0)
    else:
        myresults_InterVar.iloc[:, -1] = myresults_InterVar.iloc[:, -1].replace(element, 1)



for element in paper_results_InterVar.iloc[:, -2]:
    if element == "Benign":
        paper_results_InterVar.iloc[:, -2] = paper_results_InterVar.iloc[:, -2].replace(element, 0)
    else:
        paper_results_InterVar.iloc[:, -2] = paper_results_InterVar.iloc[:, -2].replace(element, 1)

for element in paper_results_InterVar.iloc[:, -1]:
    if element == "Benign":
        paper_results_InterVar.iloc[:, -1] = paper_results_InterVar.iloc[:, -1].replace(element, 0)
    else:
        paper_results_InterVar.iloc[:, -1] = paper_results_InterVar.iloc[:, -1].replace(element, 1)

print(myresults_ExpertPanel)
#print(myresults_InterVar)

#precision recall curve for my results prediction of benignicity

precision_exp2, recall_exp2, thresholds_exp2 = metrics.precision_recall_curve(myresults_ExpertPanel.iloc[:, -1], myresults_ExpertPanel.iloc[:, -2])
precision_exp1, recall_exp1, thresholds_exp1 = metrics.precision_recall_curve(paper_results_ExpertPanel.iloc[:, -1], paper_results_ExpertPanel.iloc[:, -2])
plt.plot(recall_exp1, precision_exp1, color = 'lightblue', label = 'v1 TAPES annotated') # plot the precision-recall curve
plt.plot(recall_exp2, precision_exp2, color = 'orange',label = 'v1 VEP annotated') # plot the precision-recall curve

plt.legend()
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve ExpertPanel data for Pathogenic prediction')
plt.show()

#precision recall curve for paper results prediction of benignicity

precision_inter2, recall_inter2, thresholds_inter2 = metrics.precision_recall_curve(myresults_InterVar.iloc[:, -1], myresults_InterVar.iloc[:, -2])
precision_inter1, recall_inter1, thresholds_inter1 = metrics.precision_recall_curve(paper_results_InterVar.iloc[:, -1], paper_results_InterVar.iloc[:, -2])
plt.plot(recall_inter1, precision_inter1, color = 'lightblue', label = 'v1 TAPES annotated') # plot the precision-recall curve
plt.plot(recall_inter2, precision_inter2, color = 'orange',label = 'v1 VEP annotated') # plot the precision-recall curve

plt.legend()
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve InterVar data for Pathogenic prediction')
plt.show()


####### Sixth step: Create the ROC curve for each of the performance analysis #######

fpr_ExP1, tpr_ExP1, threshold_ExP1 = roc_curve(paper_results_ExpertPanel.iloc[:, -2], paper_results_ExpertPanel.iloc[:, -1])
fpr_ExP2, tpr_ExP2, threshold_ExP2 = roc_curve(myresults_ExpertPanel.iloc[:, -2], myresults_ExpertPanel.iloc[:, -1])


# Calcular el área bajo la curva ROC (AUC-ROC)
auc_score1 = auc(fpr_ExP1, tpr_ExP1)
auc_score2 = auc(fpr_ExP2, tpr_ExP2)

# Graficar la curva ROC
plt.plot(fpr_ExP1, tpr_ExP1, label='ROC curve v1 TAPES annotated (AUC = %0.2f)' % auc_score1)
plt.plot(fpr_ExP2, tpr_ExP2, label='ROC curve v1 VEP annotated (AUC = %0.2f)' % auc_score2)
plt.plot([0, 1], [0, 1], 'r--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('FPS')
plt.ylabel('TPR')
plt.title('ROC-curve for Pathogenic Prediction and ExpertPanel data')
plt.legend(loc="lower right")
plt.show()





fpr_inter1, tpr_inter1, threshold_inter1 = roc_curve(paper_results_InterVar.iloc[:, -2], paper_results_InterVar.iloc[:, -1])
fpr_inter2, tpr_inter2, threshold_inter2 = roc_curve(myresults_InterVar.iloc[:, -2], myresults_InterVar.iloc[:, -1])


# Calcular el área bajo la curva ROC (AUC-ROC)
auc_score1 = auc(fpr_inter1, tpr_inter1)
auc_score2 = auc(fpr_inter2, tpr_inter2)

# Graficar la curva ROC
plt.plot(fpr_inter1, tpr_inter1, label='ROC curve v1 TAPES annotated (AUC = %0.2f)' % auc_score1)
plt.plot(fpr_inter2, tpr_inter2, label='ROC curve v1 VEP annotated (AUC = %0.2f)' % auc_score2)
plt.plot([0, 1], [0, 1], 'r--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('FPS')
plt.ylabel('TPR')
plt.title('ROC-curve for Pathogenic Prediction and InterVar data')
plt.legend(loc="lower right")
plt.show()



###### CURVES PATHOGENICITY PREDICTION #######
# here we are changing the nomenclature of the pathogenicity. Before we were considering the Benign as 0,
# now we want to evaluate the Pathogenic prediction so we changed the Pathogenic as 0.


for index, row in myresults_ExpertPanel.iterrows():
    if myresults_ExpertPanel.iloc[index, -2] == 1:
        myresults_ExpertPanel.iloc[index, -2] = 0
    else:
        myresults_ExpertPanel.iloc[index, -2] = 1
    
for index, row in myresults_ExpertPanel.iterrows():
    if myresults_ExpertPanel.iloc[index, -1] == 1:
        myresults_ExpertPanel.iloc[index, -1] = 0
    else:
        myresults_ExpertPanel.iloc[index, -1] = 1


for index, row in paper_results_ExpertPanel.iterrows():
    if paper_results_ExpertPanel.iloc[index, -2] == 1:
        paper_results_ExpertPanel.iloc[index, -2] = 0
    else:
        paper_results_ExpertPanel.iloc[index, -2] = 1

for index, row in paper_results_ExpertPanel.iterrows():
    if paper_results_ExpertPanel.iloc[index, -1] == 1:
        paper_results_ExpertPanel.iloc[index, -1] = 0
    else:
        paper_results_ExpertPanel.iloc[index, -1] = 1

for index, row in myresults_InterVar.iterrows():
    if myresults_InterVar.iloc[index, -2] == 1:
        myresults_InterVar.iloc[index, -2] = 0
    else:
        myresults_InterVar.iloc[index, -2] = 1
for index, row in myresults_InterVar.iterrows():
    if myresults_InterVar.iloc[index, -1] == 1:
        myresults_InterVar.iloc[index, -1] = 0
    else:
        myresults_InterVar.iloc[index, -1] = 1

for index, row in paper_results_InterVar.iterrows():
    if paper_results_InterVar.iloc[index, -2] == 1:
        paper_results_InterVar.iloc[index, -2] = 0
    else:
        paper_results_InterVar.iloc[index, -2] = 1
for index, row in paper_results_InterVar.iterrows():
    if paper_results_InterVar.iloc[index, -1] == 1:
        paper_results_InterVar.iloc[index, -1] = 0
    else:
        paper_results_InterVar.iloc[index, -1] = 1



precision_exp2, recall_exp2, thresholds_exp2 = metrics.precision_recall_curve(myresults_ExpertPanel.iloc[:, -1], myresults_ExpertPanel.iloc[:, -2])
precision_exp1, recall_exp1, thresholds_exp1 = metrics.precision_recall_curve(paper_results_ExpertPanel.iloc[:, -1], paper_results_ExpertPanel.iloc[:, -2])
plt.plot(recall_exp1, precision_exp1, color = 'lightblue', label = 'v1 TAPES annotated') # plot the precision-recall curve
plt.plot(recall_exp2, precision_exp2, color = 'orange',label = 'v1 VEP annotated') # plot the precision-recall curve

plt.legend()
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve ExpertPanel data for Benign prediction')
plt.show()

#precision recall curve for paper results prediction of benignicity

precision_inter2, recall_inter2, thresholds_inter2 = metrics.precision_recall_curve(myresults_InterVar.iloc[:, -1], myresults_InterVar.iloc[:, -2])
precision_inter1, recall_inter1, thresholds_inter1 = metrics.precision_recall_curve(paper_results_InterVar.iloc[:, -1], paper_results_InterVar.iloc[:, -2])
plt.plot(recall_inter1, precision_inter1, color = 'lightblue', label = 'v1 TAPES annotated') # plot the precision-recall curve
plt.plot(recall_inter2, precision_inter2, color = 'orange',label = 'v1 VEP annotated') # plot the precision-recall curve

plt.legend()
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve InterVar data for Benign prediction')
plt.show()


fpr_ExP1, tpr_ExP1, threshold_ExP1 = roc_curve(paper_results_ExpertPanel.iloc[:, -2], paper_results_ExpertPanel.iloc[:, -1])
fpr_ExP2, tpr_ExP2, threshold_ExP2 = roc_curve(myresults_ExpertPanel.iloc[:, -2], myresults_ExpertPanel.iloc[:, -1])


# Calcular el área bajo la curva ROC (AUC-ROC)
auc_score1 = auc(fpr_ExP1, tpr_ExP1)
auc_score2 = auc(fpr_ExP2, tpr_ExP2)

# Graficar la curva ROC
plt.plot(fpr_ExP1, tpr_ExP1, label='ROC curve v1 TAPES annotated (AUC = %0.2f)' % auc_score1)
plt.plot(fpr_ExP2, tpr_ExP2, label='ROC curve v1 VEP annotated (AUC = %0.2f)' % auc_score2)
plt.plot([0, 1], [0, 1], 'r--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('FPS')
plt.ylabel('TPR')
plt.title('ROC-curve for Benign Prediction and ExpertPanel data')
plt.legend(loc="lower right")
plt.show()



fpr_inter1, tpr_inter1, threshold_inter1 = roc_curve(paper_results_InterVar.iloc[:, -2], paper_results_InterVar.iloc[:, -1])
fpr_inter2, tpr_inter2, threshold_inter2 = roc_curve(myresults_InterVar.iloc[:, -2], myresults_InterVar.iloc[:, -1])


# Calcular el área bajo la curva ROC (AUC-ROC)
auc_score1 = auc(fpr_inter1, tpr_inter1)
auc_score2 = auc(fpr_inter2, tpr_inter2)

# Graficar la curva ROC
plt.plot(fpr_inter1, tpr_inter1, label='ROC curve v1 TAPES annotated (AUC = %0.2f)' % auc_score1)
plt.plot(fpr_inter2, tpr_inter2, label='ROC curve v1 VEP annotated (AUC = %0.2f)' % auc_score2)
plt.plot([0, 1], [0, 1], 'r--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('FPS')
plt.ylabel('TPR')
plt.title('ROC-curve for Benign Prediction and InterVar data')
plt.legend(loc="lower right")
plt.show()
