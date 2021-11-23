# Below is your R command history : 
ReadTabExpressData("raw.txt");
PerformDataAnnot("noAnn", "NA", "custom", "sum");
GetSummaryData()
PlotDataBox("qc_boxplot_0", "72", "png");
PlotDataPCA("qc_pca_0","72", "png", "NA");
PerformExpressNormalization("none", 0, 0,"noAnn");
PlotDataBox("qc_norm_boxplot_0", "72", "png");
PlotDataPCA("qc_norm_pca_0","72", "png", "NA");
qc.density("qc_norm_density_0","72", "png", "NA");
PerformExpressNormalization("none", 0, 0,"noAnn");
PlotDataBox("qc_norm_boxplot_1", "72", "png");
PlotDataPCA("qc_norm_pca_1","72", "png", "NA");
qc.density("qc_norm_density_1","72", "png", "NA");
SetupDesignMatrix("limma");
PerformDEAnal("reference", "dose_0.0003");
GetSigGenes("Result_1", 1.0, 0.0,"false",1,"noAnn");
PlotMAPlot("ma_plot_0","72","png","1.0","0.0","1");
PlotVolcanoPlot("volcano_plot_0","72", "png", "1.0", "0.0" ,"1");
GetMetaColLength();
nrow(dataSet$sig.mat)
PrepareDataForDoseResponse();
GetSigDRItems( 1.0,0.0,"false","false", 0.05);
InitDrcFitObj();
FilterDRFit();
PerformBMDCalc();
PlotDRModelBars("dr_barplot_0", 72, "png");
PlotDRHistogram("dr_histogram_0", 72, "png", "mg/kg", "natural");
