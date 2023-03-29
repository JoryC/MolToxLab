# Below is your R command history : 
ReadTabExpressData("raw.txt");
PerformDataAnnot("dre", "count", "emblgene", "sum");
GetSummaryData()
PlotDataBox("qc_boxplot_0", "72", "png");
PlotDataPCA("qc_pca_0","72", "png", "NA");
PerformExpressNormalization("logcount", 30, 10,"dre");
PlotDataBox("qc_norm_boxplot_0", "72", "png");
PlotDataPCA("qc_norm_pca_0","72", "png", "NA");
qc.density("qc_norm_density_0","72", "png", "NA");
SetupDesignMatrix("limma");
PerformDEAnal("reference", "DOSE_0");
GetSigGenes("Result_1", 0.05, 1.0,"true",1,"dre");
PlotMAPlot("ma_plot_0","72","png","0.05","1.0","1");
PlotVolcanoPlot("volcano_plot_0","72", "png", "0.05", "1.0" ,"1");
GetMetaColLength();
nrow(dataSet$sig.mat)
PrepareDataForDoseResponse();
GetSigDRItems( 0.05,1.0,"true","false", 0.05);
GetNumberDoses();
PerformAPIDRFit();
PerformDRFit();
PerformBMDCalc();
PerformExpressNormalization("none", 30, 10,"dre");
PlotDataBox("qc_norm_boxplot_1", "72", "png");
PlotDataPCA("qc_norm_pca_1","72", "png", "NA");
qc.density("qc_norm_density_1","72", "png", "NA");
SetupDesignMatrix("limma");
PerformDEAnal("reference", "DOSE_0");
GetSigGenes("Result_2", 0.05, 1.0,"true",1,"dre");
PlotMAPlot("ma_plot_1","72","png","0.05","1.0","1");
PlotVolcanoPlot("volcano_plot_1","72", "png", "0.05", "1.0" ,"1");
GetMetaColLength();
nrow(dataSet$sig.mat)
PerformExpressNormalization("TMM", 30, 10,"dre");
PlotDataBox("qc_norm_boxplot_2", "72", "png");
PlotDataPCA("qc_norm_pca_2","72", "png", "NA");
qc.density("qc_norm_density_2","72", "png", "NA");
SetupDesignMatrix("limma");
PerformDEAnal("reference", "DOSE_0");
GetSigGenes("Result_3", 0.05, 1.0,"true",1,"dre");
PlotMAPlot("ma_plot_2","72","png","0.05","1.0","1");
PlotVolcanoPlot("volcano_plot_2","72", "png", "0.05", "1.0" ,"1");
GetMetaColLength();
nrow(dataSet$sig.mat)
PerformExpressNormalization("logcount", 30, 10,"dre");
PlotDataBox("qc_norm_boxplot_3", "72", "png");
PlotDataPCA("qc_norm_pca_3","72", "png", "NA");
qc.density("qc_norm_density_3","72", "png", "NA");
SetupDesignMatrix("limma");
PerformDEAnal("reference", "DOSE_0");
GetSigGenes("Result_4", 0.05, 0.0,"true",1,"dre");
PlotMAPlot("ma_plot_3","72","png","0.05","0.0","1");
PlotVolcanoPlot("volcano_plot_3","72", "png", "0.05", "0.0" ,"1");
GetMetaColLength();
nrow(dataSet$sig.mat)
SetupDesignMatrix("limma");
PerformDEAnal("reference", "DOSE_0");
GetSigGenes("Result_5", 0.05, 1.0,"true",1,"dre");
PlotMAPlot("ma_plot_4","72","png","0.05","1.0","1");
PlotVolcanoPlot("volcano_plot_4","72", "png", "0.05", "1.0" ,"1");
GetMetaColLength();
nrow(dataSet$sig.mat)
PlotMAPlot("ma_plot_5","72","png","0.05","1.0","5");
PlotVolcanoPlot("volcano_plot_5","72", "png", "0.05", "1.0" ,"5");
nrow(dataSet$sig.mat)
SetupDesignMatrix("limma");
PerformDEAnal("reference", "DOSE_0");
GetSigGenes("Result_6", 0.05, 1.0,"true",5,"dre");
PlotMAPlot("ma_plot_6","72","png","0.05","1.0","5");
PlotVolcanoPlot("volcano_plot_6","72", "png", "0.05", "1.0" ,"5");
GetMetaColLength();
nrow(dataSet$sig.mat)
PlotSelectedGene("560232")
PlotSelectedGene("641564")
PlotSelectedGene("447910")
PlotSelectedGene("447910")
PlotSelectedGene("100003790")
PlotSelectedGene("557394")
SetupDesignMatrix("limma");
PerformDEAnal("reference", "DOSE_0");
GetSigGenes("Result_7", 0.05, 1.0,"false",5,"dre");
PlotMAPlot("ma_plot_7","72","png","0.05","1.0","5");
PlotVolcanoPlot("volcano_plot_7","72", "png", "0.05", "1.0" ,"5");
GetMetaColLength();
nrow(dataSet$sig.mat)
SetupDesignMatrix("limma");
PerformDEAnal("reference", "DOSE_0");
GetSigGenes("Result_8", 0.05, 1.0,"false",5,"dre");
PlotMAPlot("ma_plot_8","72","png","0.05","1.0","5");
PlotVolcanoPlot("volcano_plot_8","72", "png", "0.05", "1.0" ,"5");
GetMetaColLength();
nrow(dataSet$sig.mat)
PlotMAPlot("ma_plot_9","72","png","0.05","1.0","1");
PlotVolcanoPlot("volcano_plot_9","72", "png", "0.05", "1.0" ,"1");
nrow(dataSet$sig.mat)
PrepareDataForDoseResponse();
GetSigDRItems( 0.05,1.0,"false","false", 0.05);
PlotDRHistogram("dr_histogram_0", 72, "png", "mg/kg", "log10");
GetNumberDoses();
PerformAPIDRFit();
FilterDRFit();
PerformBMDCalc();
PlotDRModelBars("dr_barplot_0", 72, "png");
PlotDRHistogram("dr_histogram_1", 72, "png", "mg/kg", "log10");
PlotDRHistogram("dr_histogram_2", 72, "png", "mg/L", "log10");
PlotGeneDRCurve("407680","gpa33a","Exp4",0.00265188513066143,0.730753037479413,NaN,11.8039594952377,62.3,165.11,288.89,"log10")
PlotGeneDRCurve("492707","ncapg","Exp4",0.00249984045784109,0.513718039992646,NaN,9.40939387531697,57.93,234.09,438.69,"log10")
PlotGeneDRCurve("606595","dhx8","Exp4",0.00144686892002367,1.43612175144579,NaN,6.28505497614171,143.82,245.98,378.41,"log10")
PlotGeneDRCurve("100149909","si:ch211-181d7.1","Hill",1.16426226179537,2.13925859197042,6.42122100082147,3355.42423169274,326.46,1698.59,2885.5,"log10")
PlotGeneDRCurve("394015","def6a","Hill",1.49157457052608,2.63143882611913,5.81580227139929,2112.51065772545,408.12,1727.94,2640.73,"log10")
PreparePODJSON(fileNm = "pod_1.json", scale = "log10", geneDB = "kegg", org = "dre");
GetNumberDoses();
PerformAPIDRFit();
FilterDRFit();
PerformBMDCalc();
PlotDRModelBars("dr_barplot_1", 72, "png");
PlotDRHistogram("dr_histogram_3", 72, "png", "mg/L", "log10");
PreparePODJSON(fileNm = "pod_2.json", scale = "log10", geneDB = "kegg", org = "dre");
