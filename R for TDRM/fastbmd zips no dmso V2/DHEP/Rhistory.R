# Below is your R command history : 
ReadTabExpressData("raw.txt");
ReadTabExpressData("raw.txt");
PerformDataAnnot("noAnn", "NA", "custom", "sum");
GetSummaryData()
PlotDataBox("qc_boxplot_0", "72", "png");
PlotDataPCA("qc_pca_0","72", "png", "NA");
PerformExpressNormalization("none", 0, 0,"noAnn");
PlotDataBox("qc_norm_boxplot_0", "72", "png");
PlotDataPCA("qc_norm_pca_0","72", "png", "NA");
qc.density("qc_norm_density_0","72", "png", "NA");
SetupDesignMatrix("limma");
PerformDEAnal("reference", "dose_0.0003");
GetSigGenes("Result_1", 0.05, 1.0,"false",1,"noAnn");
PlotMAPlot("ma_plot_0","72","png","0.05","1.0","1");
PlotVolcanoPlot("volcano_plot_0","72", "png", "0.05", "1.0" ,"1");
GetMetaColLength();
nrow(dataSet$sig.mat)
PlotMAPlot("ma_plot_1","72","png","0.05","1.0","4");
PlotVolcanoPlot("volcano_plot_1","72", "png", "0.05", "1.0" ,"4");
nrow(dataSet$sig.mat)
PlotMAPlot("ma_plot_2","72","png","0.05","1.0","3");
PlotVolcanoPlot("volcano_plot_2","72", "png", "0.05", "1.0" ,"3");
nrow(dataSet$sig.mat)
PlotMAPlot("ma_plot_3","72","png","0.05","1.0","2");
PlotVolcanoPlot("volcano_plot_3","72", "png", "0.05", "1.0" ,"2");
nrow(dataSet$sig.mat)
PrepareDataForDoseResponse();
GetSigDRItems( 0.05,1.0,"false","false", 0.05);
InitDrcFitObj();
FilterDRFit();
PerformBMDCalc();
PlotDRModelBars("dr_barplot_0", 72, "png");
PlotDRHistogram("dr_histogram_0", 72, "png", "mg/kg", "natural");
PreparePDFReport("doseresponse", "guest17443610682709879313", "Fit Curves")

ReadTabExpressData("raw.txt");
PerformDataAnnot("noAnn", "NA", "custom", "sum");
GetSummaryData()
PlotDataBox("qc_boxplot_1", "72", "png");
PlotDataPCA("qc_pca_1","72", "png", "NA");
PerformExpressNormalization("none", 0, 0,"noAnn");
PlotDataBox("qc_norm_boxplot_1", "72", "png");
PlotDataPCA("qc_norm_pca_1","72", "png", "NA");
qc.density("qc_norm_density_1","72", "png", "NA");
SetupDesignMatrix("limma");
PerformDEAnal("reference", "dose_0.0003");
GetSigGenes("Result_2", 0.05, 1.0,"false",2,"noAnn");
PlotMAPlot("ma_plot_4","72","png","0.05","1.0","2");
PlotVolcanoPlot("volcano_plot_4","72", "png", "0.05", "1.0" ,"2");
GetMetaColLength();
nrow(dataSet$sig.mat)
PrepareDataForDoseResponse();
GetSigDRItems( 0.05,1.0,"false","false", 0.05);
InitDrcFitObj();
FilterDRFit();
PerformBMDCalc();
PlotDRModelBars("dr_barplot_1", 72, "png");
PlotDRHistogram("dr_histogram_1", 72, "png", "mg/kg", "natural");
ReadTabExpressData("raw.txt");
PerformDataAnnot("noAnn", "NA", "custom", "sum");
GetSummaryData()
ReadTabExpressData("raw.txt");
PerformDataAnnot("noAnn", "NA", "custom", "sum");
GetSummaryData()
PlotDataBox("qc_boxplot_2", "72", "png");
PlotDataPCA("qc_pca_2","72", "png", "NA");
PerformExpressNormalization("none", 0, 0,"noAnn");
PlotDataBox("qc_norm_boxplot_2", "72", "png");
PlotDataPCA("qc_norm_pca_2","72", "png", "NA");
qc.density("qc_norm_density_2","72", "png", "NA");
SetupDesignMatrix("limma");
PerformDEAnal("reference", "dose_0.0003");
GetSigGenes("Result_3", 0.05, 1.0,"false",2,"noAnn");
PlotMAPlot("ma_plot_5","72","png","0.05","1.0","2");
PlotVolcanoPlot("volcano_plot_5","72", "png", "0.05", "1.0" ,"2");
GetMetaColLength();
nrow(dataSet$sig.mat)
PrepareDataForDoseResponse();
GetSigDRItems( 0.05,1.0,"false","false", 0.05);
InitDrcFitObj();
FilterDRFit();
PerformBMDCalc();
PlotDRModelBars("dr_barplot_2", 72, "png");
PlotDRHistogram("dr_histogram_2", 72, "png", "mg/kg", "natural");
PlotGeneDRCurve("ENSDARG00000012275","ENSDARG00000012275","Exp4",257.204725142756,0.763683316121233,0.0,4.76290295021803,0.0,0.0,0.0,"natural")
PlotDRHistogram("dr_histogram_3", 72, "png", "mg/kg", "log10");
PlotGeneDRCurve("ENSDARG00000012275","ENSDARG00000012275","Exp4",257.204725142756,0.763683316121233,0.0,4.76290295021803,0.0,0.0,0.0,"log10")
PlotGeneDRCurve("ENSDARG00000097890","ENSDARG00000097890","Exp3",-0.0056385330139892,0.0,0.200379262434537,5.31298300130767,0.0,0.01,0.02,"log10")
ReadTabExpressData("raw.txt");
PerformDataAnnot("noAnn", "NA", "custom", "sum");
GetSummaryData()
PlotDataBox("qc_boxplot_3", "72", "png");
PlotDataPCA("qc_pca_3","72", "png", "NA");
PerformExpressNormalization("none", 0, 0,"noAnn");
PlotDataBox("qc_norm_boxplot_3", "72", "png");
PlotDataPCA("qc_norm_pca_3","72", "png", "NA");
qc.density("qc_norm_density_3","72", "png", "NA");
SetupDesignMatrix("limma");
PerformDEAnal("reference", "dose_0.0003");
GetSigGenes("Result_4", 0.05, 1.0,"false",2,"noAnn");
PlotMAPlot("ma_plot_6","72","png","0.05","1.0","2");
PlotVolcanoPlot("volcano_plot_6","72", "png", "0.05", "1.0" ,"2");
GetMetaColLength();
nrow(dataSet$sig.mat)
PrepareDataForDoseResponse();
GetSigDRItems( 0.05,1.0,"false","false", 0.05);
InitDrcFitObj();
FilterDRFit();
PerformBMDCalc();
PlotDRModelBars("dr_barplot_3", 72, "png");
PlotDRHistogram("dr_histogram_4", 72, "png", "mg/kg", "log10");
