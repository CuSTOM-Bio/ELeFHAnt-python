import os
from sys import argv
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as an
import bbknn
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import gseapy

from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score
from sklearn.metrics import confusion_matrix
from sklearn import preprocessing
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import VarianceThreshold
from sklearn.utils.multiclass import unique_labels

warnings.simplefilter(action='ignore', category=FutureWarning)

def CelltypeAnnotation(reference = None, query = None, downsample = False, downsample_to = 200, classification_method = "Ensemble", crossvalidationSVM = 10, validatePredictions = True, selectvarfeatures = 2000, ntree = 500, cost= = 10, classification_approach = "ClassifyCells"):
    print("Running diagnostics on reference and query")
    print("Number of cells in reference: %s" % str(reference.shape[0]))
    print("Number of cells in query: %s" % str(query.shape[0]))
    if downsample:
        print("Downsampling reference")
        reference_use = reference[reference.obs.groupby('Celltypes',group_keys=False).apply(lambda x: x.sample(n=min(downsample_to, len(x)))).index]
        print("Number of cells in reference after downsampling per celltype: %s" % str(reference_use.shape[0]))
        print("Ratio of number of cells in query vs downsampled reference: %s" % str(query.shape[0]/reference_use.shape[0]))
        if classification_approach == "ClassifyCells":
            query = ClassifyCells(reference = reference_use, query = query, downsample = downsample, downsample_to = downsample_to, classification_method = classification_method, crossvalidationSVM = crossvalidationSVM, validatePredictions = validatePredictions, selectvarfeatures = selectvarfeatures, ntree = ntree, cost = cost)
            return(query)
        elif classification_approach == "ClassifyCells_usingApproximation":
            query = ApproximationBasedCelltypeAssignment(reference = reference_use, query = query, downsample = downsample, downsample_to = downsample_to, classification_method = classification_method, crossvalidationSVM = crossvalidationSVM, validatePredictions = validatePredictions, selectvarfeatures = selectvarfeatures, ntree = ntree, cost = cost)
            return(query)
    else:
        print("Ratio of number of cells in query vs  reference: %s" % str(query.shape[0]/reference.shape[0]))
        if classification_approach == "ClassifyCells":
            query = ClassifyCells(reference = reference, query = query, downsample = downsample, downsample_to = downsample_to, classification_method = classification_method, crossvalidationSVM = crossvalidationSVM, validatePredictions = validatePredictions, selectvarfeatures = selectvarfeatures, ntree = ntree, cost = cost)
            return(query)
        elif classification_approach == "ClassifyCells_usingApproximation":
            query = ApproximationBasedCelltypeAssignment(reference = reference, query = query, downsample = downsample, downsample_to = downsample_to, classification_method = classification_method, crossvalidationSVM = crossvalidationSVM, validatePredictions = validatePredictions, selectvarfeatures = selectvarfeatures, ntree = ntree, cost = cost)
            return(query)

def ClassifyCells(reference, query, downsample, downsample_to, classification_method, crossvalidationSVM, validatePredictions, selectvarfeatures, ntree, cost):
    query_use = query
    reference_use = reference
    print("Merging reference and query")
    combined = reference_use.concatenate(query_use)
    print("Normalization, variable feature selection, and scaling")
    sc.pp.normalize_total(combined)
    sc.pp.log1p(combined)
    sc.pp.highly_variable_genes(combined,n_top_genes=selectvarfeatures)
    combined = combined[:, combined.var.highly_variable]
    sc.pp.scale(combined)
    print("Number of features selected: %s" % str(combined.shape[1]))
    x_train = combined[combined.obs["batch"] == "0",]
    x_test = combined[combined.obs["batch"] == "1",]
    if classification_method == "randomForest":
        rf_pred, rf_acc = randomForest_predictor(train = x_train.X, test = x_test.X, train_label = x_train.obs["Celltypes"], test_label = x_test.obs["Clusters"], ntree = ntree)
        query_use.obs["PredictedCelltype_UsingRF"] = rf_pred
        print("Added predicted celltypes using randomForest to query")
        if validatePredictions:
            print("Starting validation of cell type assignments using GSEA")
            reference_validation_use = reference_use[reference_use.obs['Celltypes'].isin(np.unique(query_use.obs["PredictedCelltype_UsingRF"]))]
            validation = ValidatePredictions(reference = reference_validation_use, query = query_use)
            validation.to_csv("Summary_GeneSetEnrichmentAnalysis.txt")
    elif classification_method == "SVM":
        svm_pred, svm_acc = svm_predictor(train = x_train.X, test = x_test.X, train_label = x_train.obs["Celltypes"], test_label = x_test.obs["Clusters"], crossvalidationSVM = crossvalidationSVM, cost = cost)
        query_use.obs["PredictedCelltype_UsingSVM"] = svm_pred
        print("Added predicted celltypes using SVM to query")
        if validatePredictions:
            print("Starting validation of cell type assignments using GSEA")
            reference_validation_use = reference_use[reference_use.obs['Celltypes'].isin(np.unique(query_use.obs["PredictedCelltype_UsingSVM"]))]
            validation = ValidatePredictions(reference = reference_validation_use, query = query_use)
            validation.to_csv("Summary_GeneSetEnrichmentAnalysis.txt")
    elif classification_method == "Ensemble":
        rf_pred, rf_acc = randomForest_predictor(train = x_train.X, test = x_test.X, train_label = x_train.obs["Celltypes"], test_label = x_test.obs["Clusters"], ntree = ntree)
        query_use.obs["PredictedCelltype_UsingRF"] = rf_pred
        print("Added predicted celltypes using randomForest to query")
        svm_pred, svm_acc = svm_predictor(train = x_train.X, test = x_test.X, train_label = x_train.obs["Celltypes"], test_label = x_test.obs["Clusters"], crossvalidationSVM = crossvalidationSVM, cost = cost)
        query_use.obs["PredictedCelltype_UsingSVM"] = svm_pred
        print("Added predicted celltypes using SVM to query")
        ensemble_pred = []
        for x in range(len(rf_pred)):
            if rf_pred[x] == svm_pred[x]:
                ensemble_pred.append(rf_pred[x])
            else:
                if rf_acc >= svm_acc:
                    ensemble_pred.append(rf_pred[x])
                else:
                    ensemble_pred.append(svm_pred[x])
        query_use.obs["PredictedCelltype_UsingEnsemble"] = ensemble_pred
        print("Added predicted celltypes using ensemble to query")
        if validatePredictions:
            print("Starting validation of cell type assignments using GSEA")
            reference_validation_use = reference_use[reference_use.obs['Celltypes'].isin(np.unique(query_use.obs["PredictedCelltype_UsingEnsemble"]))]
            validation = ValidatePredictions(reference = reference_validation_use, query = query_use)
            validation.to_csv("Summary_GeneSetEnrichmentAnalysis.txt")
    return(query_use)


def ApproximationBasedCelltypeAssignment(reference, query, downsample, downsample_to, classification_method, crossvalidationSVM, validatePredictions, selectvarfeatures, ntree, cost):
    query_use =  query[query.obs.groupby('Clusters',group_keys=False).apply(lambda x: x.sample(n=min(downsample_to, len(x)))).index]
    reference_use = reference
    print("Merging reference and query")
    combined = reference_use.concatenate(query_use)
    print("Normalization, variable feature selection, and scaling")
    sc.pp.normalize_total(combined)
    sc.pp.log1p(combined)
    sc.pp.highly_variable_genes(combined,n_top_genes=selectvarfeatures)
    combined = combined[:, combined.var.highly_variable]
    sc.pp.scale(combined)
    print("Number of features selected: %s" % str(combined.shape[1]))
    x_train = combined[combined.obs["batch"] == "0",]
    x_test = combined[combined.obs["batch"] == "1",]
    if classification_method == "randomForest":
        rf_pred, rf_acc = randomForest_predictor(train = x_train.X, test = x_test.X, train_label = x_train.obs["Celltypes"], test_label = x_test.obs["Clusters"], ntree = ntree)
        rf_cm = pd.crosstab(x_test.obs["Clusters"],rf_pred)
        rf_pred_approx = list(map(lambda x: rf_cm.idxmax(axis="columns").tolist()[int(rf_cm.index.values.tolist().index(x))],x_test.obs["Clusters"]))
        query_use.obs["PredictedCelltype_UsingRF"] = rf_pred_approx
        print("Added predicted celltypes using randomForest to query")
        rf_cm.to_csv("ConfusionMatrix_RF.txt")
        if validatePredictions:
            print("Starting validation of cell type assignments using GSEA")
            reference_validation_use = reference_use[reference_use.obs['Celltypes'].isin(np.unique(query_use.obs["PredictedCelltype_UsingRF"]))]
            validation = ValidatePredictions(reference = reference_validation_use, query = query_use)
            validation.to_csv("Summary_GeneSetEnrichmentAnalysis.txt")
    elif classification_method == "SVM":
        svm_pred, svm_acc = svm_predictor(train = x_train.X, test = x_test.X, train_label = x_train.obs["Celltypes"], test_label = x_test.obs["Clusters"], crossvalidationSVM = crossvalidationSVM, cost = cost)
        svm_cm = pd.crosstab(x_test.obs["Clusters"],svm_pred)
        svm_pred_approx = list(map(lambda x: svm_cm.idxmax(axis="columns").tolist()[int(svm_cm.index.values.tolist().index(x))],x_test.obs["Clusters"]))
        query_use.obs["PredictedCelltype_UsingSVM"] = svm_pred_approx
        print("Added predicted celltypes using SVM to query")
        svm_cm.to_csv("ConfusionMatrix_SVM.txt")
        if validatePredictions:
            print("Starting validation of cell type assignments using GSEA")
            reference_validation_use = reference_use[reference_use.obs['Celltypes'].isin(np.unique(query_use.obs["PredictedCelltype_UsingSVM"]))]
            validation = ValidatePredictions(reference = reference_validation_use, query = query_use)
            validation.to_csv("Summary_GeneSetEnrichmentAnalysis.txt")
    elif classification_method == "Ensemble":
        rf_pred, rf_acc = randomForest_predictor(train = x_train.X, test = x_test.X, train_label = x_train.obs["Celltypes"], test_label = x_test.obs["Clusters"], ntree = ntree)
        rf_cm = pd.crosstab(x_test.obs["Clusters"],rf_pred)
        rf_pred_approx = list(map(lambda x: rf_cm.idxmax(axis="columns").tolist()[int(rf_cm.index.values.tolist().index(x))],x_test.obs["Clusters"]))
        query_use.obs["PredictedCelltype_UsingRF"] = rf_pred_approx
        print("Added predicted celltypes using randomForest to query")
        rf_cm.to_csv("ConfusionMatrix_RF.txt")
        svm_pred, svm_acc = svm_predictor(train = x_train.X, test = x_test.X, train_label = x_train.obs["Celltypes"], test_label = x_test.obs["Clusters"], crossvalidationSVM = crossvalidationSVM, cost = cost)
        svm_cm = pd.crosstab(x_test.obs["Clusters"],svm_pred)
        svm_pred_approx = list(map(lambda x: svm_cm.idxmax(axis="columns").tolist()[int(svm_cm.index.values.tolist().index(x))],x_test.obs["Clusters"]))
        query_use.obs["PredictedCelltype_UsingSVM"] = svm_pred_approx
        print("Added predicted celltypes using SVM to query")
        svm_cm.to_csv("ConfusionMatrix_SVM.txt")
        print("Combining into final ensemble predictions")
        rf_cm = rf_cm * rf_acc
        svm_cm = svm_cm * svm_acc
        consensus_cm = (rf_cm/rf_cm.max().max()).add(svm_cm/svm_cm.max().max(),fill_value=0)
        ensemble_pred_approx = list(map(lambda x: consensus_cm.idxmax(axis="columns").tolist()[int(consensus_cm.index.values.tolist().index(x))],x_test.obs["Clusters"]))
        query_use.obs["PredictedCelltype_UsingEnsemble"] = ensemble_pred_approx
        print("Added predicted celltypes using ensemble to query")
        consensus_cm.to_csv("ConfusionMatrix_Ensemble.txt")
        if validatePredictions:
            print("Starting validation of cell type assignments using GSEA")
            reference_validation_use = reference_use[reference_use.obs['Celltypes'].isin(np.unique(query_use.obs["PredictedCelltype_UsingEnsemble"]))]
            validation = ValidatePredictions(reference = reference_validation_use, query = query_use)
            validation.to_csv("Summary_GeneSetEnrichmentAnalysis.txt")

    return(query_use)


def LabelHarmonization (data_objects = None, perform_integration = True, integrated_atlas = None, downsample = True, downsample_to = 100, npcs = 30, resolution = 0.5, classification_method = "Ensemble", crossvalidationSVM = 10, validatePredictions = True, ntree = 500, cost = 10, selectvarfeatures=2000):
    if perform_integration:
    	print("Integrating datasets")
        if downsample:
            for x in range(len(data_objects)):
                data_objects[x] = data_objects[x][data_objects[x].obs.groupby('Celltypes',group_keys=False).apply(lambda x: x.sample(n=min(downsample_to, len(x)))).index]
        integrated = an.AnnData.concatenate(*data_objects)
        sc.pp.normalize_total(integrated)
        sc.pp.log1p(integrated)
        sc.pp.highly_variable_genes(integrated,n_top_genes=selectvarfeatures)
        integrated = integrated[:, integrated.var.highly_variable]
        sc.pp.scale(integrated)
        sc.tl.pca(integrated)
        sc.external.pp.bbknn(integrated, batch_key='batch')
        sc.tl.umap(integrated)
        sc.tl.leiden(integrated, resolution=resolution)
        integrated.obs["Clusters"] = integrated.obs["leiden"]
    else:
        if downsample:
            integrated = integrated_atlas[integrated_atlas.obs.groupby('Celltypes',group_keys=False).apply(lambda x: x.sample(n=min(downsample_to, len(x)))).index]
        else:
            integrated = integrated_atlas
    print("Generating train and test datasets using stratification -- 60% for training & 40% for testing")
    integrated_data = integrated.X
    integrated_clusters = integrated.obs["Clusters"]
    integrated_celltypes = integrated.obs["Celltypes"]
    X_train,X_test,clusters_train,clusters_test,cells_train,cells_test = train_test_split(integrated_data,integrated_clusters,integrated_celltypes,test_size=0.4, stratify=integrated_clusters,random_state=0)
    if classification_method == "randomForest":
        rf_pred, rf_acc = randomForest_predictor(train = X_train, test = X_test, train_label = cells_train, test_label = clusters_test, ntree = ntree)
        rf_cm = pd.crosstab(clusters_test,rf_pred)
        rf_pred_approx = list(map(lambda x: rf_cm.idxmax(axis="columns").tolist()[int(rf_cm.index.values.tolist().index(x))],integrated_clusters))
        integrated.obs["HarmonizedLabels_UsingRF"] = rf_pred_approx
        print("Added harmonized labels using randomForest to integrated object")
        rf_cm.to_csv("ConfusionMatrix_RF.txt")
        if validatePredictions:
            print("Starting validation of cell type assignments using GSEA")
            integrated_validation_use = integrated[integrated.obs['Celltypes'].isin(np.unique(integrated.obs["HarmonizedLabels_UsingRF"]))]
            validation = ValidatePredictions(reference = integrated_validation_use, query = integrated)
            validation.to_csv("Summary_GeneSetEnrichmentAnalysis.txt")
    elif classification_method == "SVM":
        svm_pred, svm_acc = svm_predictor(train = X_train, test = X_test, train_label = cells_train, test_label = clusters_test, crossvalidationSVM = crossvalidationSVM, cost = cost)
        svm_cm = pd.crosstab(clusters_test,svm_pred)
        svm_pred_approx = list(map(lambda x: svm_cm.idxmax(axis="columns").tolist()[int(svm_cm.index.values.tolist().index(x))],integrated_clusters))
        integrated.obs["HarmonizedLabels_UsingSVM"] = svm_pred_approx
        print("Added harmonized labels using SVM to integrated object")
        svm_cm.to_csv("ConfusionMatrix_SVM.txt")
        if validatePredictions:
            print("Starting validation of cell type assignments using GSEA")
            integrated_validation_use = integrated[integrated.obs['Celltypes'].isin(np.unique(integrated.obs["HarmonizedLabels_UsingSVM"]))]
            validation = ValidatePredictions(reference = integrated_validation_use, query = integrated)
            validation.to_csv("Summary_GeneSetEnrichmentAnalysis.txt")
    elif classification_method == "Ensemble":
        rf_pred, rf_acc = randomForest_predictor(train = X_train, test = X_test, train_label = cells_train, test_label = clusters_test, ntree = ntree)
        rf_cm = pd.crosstab(clusters_test,rf_pred)
        rf_pred_approx = list(map(lambda x: rf_cm.idxmax(axis="columns").tolist()[int(rf_cm.index.values.tolist().index(x))],integrated_clusters))
        integrated.obs["HarmonizedLabels_UsingRF"] = rf_pred_approx
        print("Added harmonized labels using randomForest to integrated object")
        rf_cm.to_csv("ConfusionMatrix_RF.txt")
        svm_pred, svm_acc = svm_predictor(train = X_train, test = X_test, train_label = cells_train, test_label = clusters_test, crossvalidationSVM = crossvalidationSVM, cost = cost)
        svm_cm = pd.crosstab(clusters_test,svm_pred)
        svm_pred_approx = list(map(lambda x: svm_cm.idxmax(axis="columns").tolist()[int(svm_cm.index.values.tolist().index(x))],integrated_clusters))
        integrated.obs["HarmonizedLabels_UsingSVM"] = svm_pred_approx
        print("Added harmonized labels using SVM to integrated object")
        svm_cm.to_csv("ConfusionMatrix_SVM.txt")
        print("Combining into final ensemble predictions")
        rf_cm = rf_cm * rf_acc
        svm_cm = svm_cm * svm_acc
        consensus_cm = (rf_cm/rf_cm.max().max()).add(svm_cm/svm_cm.max().max(),fill_value=0)
        ensemble_pred_approx = list(map(lambda x: consensus_cm.idxmax(axis="columns").tolist()[int(consensus_cm.index.values.tolist().index(x))],integrated_clusters))
        integrated.obs["HarmonizedLabels_UsingEnsemble"] = ensemble_pred_approx
        print("Added harmonized labels using ensemble to integrated object")
        consensus_cm.to_csv("ConfusionMatrix_Ensemble.txt")
        if validatePredictions:
            print("Starting validation of cell type assignments using GSEA")
            integrated_validation_use = integrated[integrated.obs['Celltypes'].isin(np.unique(integrated.obs["HarmonizedLabels_UsingEnsemble"]))]
            validation = ValidatePredictions(reference = integrated_validation_use, query = integrated)
            validation.to_csv("Summary_GeneSetEnrichmentAnalysis.txt")
    return(integrated)

def DeduceRelationship(reference1 = None, reference2 = None, downsample = False, downsample_to = 100, classification_method = "Ensemble", crossvalidationSVM = 5, selectvarfeatures = 2000, ntree = 500, cost = 10):
    if downsample:
        print("Downsampling reference")
        reference1_use = reference1[reference1.obs.groupby('Celltypes',group_keys=False).apply(lambda x: x.sample(n=min(downsample_to, len(x)))).index]
        reference2_use = reference2[reference2.obs.groupby('Celltypes',group_keys=False).apply(lambda x: x.sample(n=min(downsample_to, len(x)))).index]
    else:
        reference1_use = reference1
        reference2_use = reference2
    print("Merging references")
    combined = reference1_use.concatenate(reference2_use)
    print("Normalization, variable feature selection, and scaling")
    sc.pp.normalize_total(combined)
    sc.pp.log1p(combined)
    sc.pp.highly_variable_genes(combined,n_top_genes=selectvarfeatures)
    combined = combined[:, combined.var.highly_variable]
    sc.pp.scale(combined)
    print("Number of features selected: %s" % str(combined.shape[1]))
    x_train = combined[combined.obs["batch"] == "0",]
    x_test = combined[combined.obs["batch"] == "1",]
    if classification_method == "randomForest":
        rf_pred, rf_acc = randomForest_predictor(train = x_train.X, test = x_test.X, train_label = x_train.obs["Celltypes"], test_label = x_test.obs["Celltypes"], ntree = ntree)
        rf_cm = pd.crosstab(x_test.obs["Celltypes"],rf_pred)
        rf_cm.to_csv("ConfusionMatrix_RF.txt")
        rf_cm_norm = rf_cm.divide(rf_cm.max(axis=1),axis=0).round(3)
        heatmap = sns.clustermap(rf_cm_norm.transpose(),xticklabels=True,yticklabels=True,cmap="vlag",cbar_pos=(0, .4, .03, .4))
        heatmap.ax_row_dendrogram.set_visible(False)
        heatmap.ax_col_dendrogram.set_visible(False)
        heatmap.ax_heatmap.set_ylabel("reference 1")
        heatmap.ax_heatmap.set_xlabel("reference 2")
        plt.savefig('randomForest_Heatmap.png')
    elif classification_method == "SVM":
        svm_pred, svm_acc = svm_predictor(train = x_train.X, test = x_test.X, train_label = x_train.obs["Celltypes"], test_label = x_test.obs["Celltypes"], crossvalidationSVM = crossvalidationSVM, cost = cost)
        svm_cm = pd.crosstab(x_test.obs["Celltypes"],svm_pred)
        svm_cm.to_csv("ConfusionMatrix_SVM.txt")
        svm_cm_norm = svm_cm.divide(svm_cm.max(axis=1),axis=0).round(3)
        heatmap = sns.clustermap(svm_cm_norm.transpose(),xticklabels=True,yticklabels=True,cmap="vlag",cbar_pos=(0, .4, .03, .4))
        heatmap.ax_row_dendrogram.set_visible(False)
        heatmap.ax_col_dendrogram.set_visible(False)
        heatmap.ax_heatmap.set_ylabel("reference 1")
        heatmap.ax_heatmap.set_xlabel("reference 2")
        plt.savefig('SVM_Heatmap.png')
    elif classification_method == "Ensemble":
        rf_pred, rf_acc = randomForest_predictor(train = x_train.X, test = x_test.X, train_label = x_train.obs["Celltypes"], test_label = x_test.obs["Celltypes"], ntree = ntree)
        rf_cm = pd.crosstab(x_test.obs["Celltypes"],rf_pred)
        rf_cm.to_csv("ConfusionMatrix_RF.txt")
        svm_pred, svm_acc = svm_predictor(train = x_train.X, test = x_test.X, train_label = x_train.obs["Celltypes"], test_label = x_test.obs["Celltypes"], crossvalidationSVM = crossvalidationSVM, cost = cost)
        svm_cm = pd.crosstab(x_test.obs["Celltypes"],svm_pred)
        svm_cm.to_csv("ConfusionMatrix_SVM.txt")
        print("Combining into final ensemble predictions")
        rf_cm = rf_cm * rf_acc
        svm_cm = svm_cm * svm_acc
        consensus_cm = (rf_cm/rf_cm.max().max()).add(svm_cm/svm_cm.max().max(),fill_value=0)
        consensus_cm.to_csv("ConfusionMatrix_Ensemble.txt")
        consensus_cm_norm = consensus_cm.divide(consensus_cm.max(axis=1),axis=0).round(3)
        heatmap = sns.clustermap(consensus_cm_norm.transpose(),xticklabels=True,yticklabels=True,cmap="vlag",cbar_pos=(0, .4, .03, .4))
        heatmap.ax_row_dendrogram.set_visible(False)
        heatmap.ax_col_dendrogram.set_visible(False)
        heatmap.ax_heatmap.set_ylabel("reference 1")
        heatmap.ax_heatmap.set_xlabel("reference 2")
        plt.savefig('Ensemble_Heatmap.png')
    return(heatmap)

def ValidatePredictions (reference = None, query = None):
    sc.tl.rank_genes_groups(reference, 'Celltypes', method='wilcoxon')
    celltype_markers = pd.DataFrame(reference.uns['rank_genes_groups']['names']).head(100)
    gene_list = celltype_markers.to_dict("list")
    sc.tl.rank_genes_groups(query, 'Clusters', method='wilcoxon')
    cluster_markers = pd.DataFrame(query.uns['rank_genes_groups']['names']).head(100)
    result = query.uns['rank_genes_groups']
    groups = list(result['names'].dtype.names)
    cluster_data = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'logfoldchanges']})
    gsea_combined = []
    for x in range(0, len(cluster_data.columns), 2):
        subset = cluster_data.iloc[:,x:x+2]
        subset.columns = [0,1]
        subset.sort_values(by=subset.columns[1],ascending=False, inplace=True)
        gsea_res = gseapy.prerank(rnk=subset, gene_sets=gene_list, permutation_num=10, min_size=10, max_size=1000)
        gsea_res = gsea_res.res2d.sort_values(by="nes",key=abs,ascending=False).head(10)
        gsea_res["Cluster"] = groups.pop(0)
        gsea_combined.append(gsea_res)
    return(pd.concat(gsea_combined))

def randomForest_predictor(train = None, test = None, train_label = None, test_label = None, ntree = None):
    clf = RandomForestClassifier(n_estimators=ntree,oob_score=True)
    print("Training randomForest classifier")
    clf.fit(train,train_label)
    print("OOB score: %s" % str(clf.oob_score_))
    print("Predicting using randomForest classifier")
    predicted = clf.predict(test)
    return predicted,clf.oob_score_

def svm_predictor(train = None, test = None, train_label = None, test_label = None, crossvalidationSVM = None, cost = None):
    Classifier = LinearSVC(max_iter=10000, C=cost, dual=True)
    clf = CalibratedClassifierCV(Classifier)
    accuracy = 0
    if crossvalidationSVM > 0:
        print("Performing cross-validation for SVM")
        skf = StratifiedKFold(n_splits=crossvalidationSVM)
        true_label = []
        pred_label = []
        for train_index, test_index in skf.split(train, train_label):
            x_train, x_test = train[train_index], train[test_index]
            y_train, y_test = train_label[train_index], train_label[test_index]
            clf.fit(x_train, y_train)
            cv_predict = clf.predict(x_test)
            true_label.extend(y_test.values)
            pred_label.extend(cv_predict)
        accuracy = accuracy_score(true_label, pred_label)
        print("Accuracy: %s" % str(accuracy))
    print("Training SVM classifier")
    clf.fit(train, train_label)
    print("Predicting using SVM classifier")
    predicted = clf.predict(test)
    return predicted,accuracy