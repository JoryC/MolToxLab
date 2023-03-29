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
SetupDesignMatrix("limma");
PerformDEAnal("reference", "DOSE_0");
GetSigGenes("Result_2", 0.05, 1.0,"false",1,"dre");
PlotMAPlot("ma_plot_1","72","png","0.05","1.0","1");
PlotVolcanoPlot("volcano_plot_1","72", "png", "0.05", "1.0" ,"1");
GetMetaColLength();
nrow(dataSet$sig.mat)
SetupDesignMatrix("limma");
PerformDEAnal("reference", "DOSE_0");
GetSigGenes("Result_3", 0.05, 0.5,"false",1,"dre");
PlotMAPlot("ma_plot_2","72","png","0.05","0.5","1");
PlotVolcanoPlot("volcano_plot_2","72", "png", "0.05", "0.5" ,"1");
GetMetaColLength();
nrow(dataSet$sig.mat)
PlotSelectedGene("30518")
PlotSelectedGene("64260")
PlotSelectedGene("559229")
PlotSelectedGene("559931")
PlotSelectedGene("567675")
SetupDesignMatrix("limma");
PerformDEAnal("reference", "DOSE_0");
GetSigGenes("Result_4", 0.05, 1.0,"false",1,"dre");
PlotMAPlot("ma_plot_3","72","png","0.05","1.0","1");
PlotVolcanoPlot("volcano_plot_3","72", "png", "0.05", "1.0" ,"1");
GetMetaColLength();
nrow(dataSet$sig.mat)
PrepareDataForDoseResponse();
GetSigDRItems( 0.05,1.0,"false","false", 0.05);
PlotDRHistogram("dr_histogram_0", 72, "png", "mg/L", "natural");
PlotDRHistogram("dr_histogram_1", 72, "png", "mg/L", "log10");
GetNumberDoses();
PerformAPIDRFit();
FilterDRFit();
PerformBMDCalc();
PlotDRModelBars("dr_barplot_0", 72, "png");
PlotDRHistogram("dr_histogram_2", 72, "png", "mg/L", "log10");
PreparePODJSON(fileNm = "pod_1.json", scale = "log10", geneDB = "kegg", org = "dre");
