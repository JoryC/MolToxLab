# Below is your R command history : 
ReadTabExpressData("raw.txt");
PerformDataAnnot("dre", "NA", "emblgene", "sum");
GetSummaryData()
PlotDataBox("qc_boxplot_0", "72", "png");
PlotDataPCA("qc_pca_0","72", "png", "NA");
PerformExpressNormalization("none", 30, 15,"dre");
PlotDataBox("qc_norm_boxplot_0", "72", "png");
PlotDataPCA("qc_norm_pca_0","72", "png", "NA");
qc.density("qc_norm_density_0","72", "png", "NA");
SetupDesignMatrix("limma");
PerformDEAnal("reference", "dose_0");
GetSigGenes("Result_1", 0.05, 1.0,"true",1,"dre");
PlotMAPlot("ma_plot_0","72","png","0.05","1.0","1");
PlotVolcanoPlot("volcano_plot_0","72", "png", "0.05", "1.0" ,"1");
GetMetaColLength();
nrow(dataSet$sig.mat)
PrepareDataForDoseResponse();
GetSigDRItems( 0.05,1.0,"false","false", 0.05);
InitDrcFitObj();
FilterDRFit();
PerformBMDCalc();
PlotDRModelBars("dr_barplot_0", 72, "png");
PlotDRHistogram("dr_histogram_0", 72, "png", "mg/kg", "natural");
PlotDRHistogram("dr_histogram_1", 72, "png", "mg/kg", "log2");
