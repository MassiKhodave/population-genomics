#### GF14h ####
# Ensure these match the correct data file names: 
input="GF14h_logFoldChange.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="geneUniverse.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
source("gomwu.functions.R")


# ------------- Calculating stats
#### Molecular Function ####
goDivision="MF" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Strawberry/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.

# ----------- Plotting results

results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.1,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
results[[1]]

# ----------- Simply repeat for other GO categories
#### Biological Process ####
goDivision="BP"

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Strawberry/perl/bin/perl.exe", 
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25,
)

results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.1,10),
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5,
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
)
results[[1]]

#### Cellular Component ####
goDivision="CC"

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Strawberry/perl/bin/perl.exe", 
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25,
)

results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.1,10),
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5,
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
)
results[[1]]

#### HDAC19 ####
# Ensure these match the correct data file names: 
input="HDAC19_logFoldChange.csv"

#### Molecular Function ####
goDivision="MF"

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Strawberry/perl/bin/perl.exe", 
           largest=0.1,  
           smallest=5,
           clusterCutHeight=0.25, 
)

results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.1,10),
                  level1=0.1, 
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5,
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
)
results[[1]]

#### Biological Process ####
goDivision="BP"

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Strawberry/perl/bin/perl.exe", 
           largest=0.01, # Reduced fraction size such that less GO terms would be on dendrogram - It was too large with fractional cutoff of 0.1
           smallest=5,
           clusterCutHeight=0.25,
)

results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.1,10),
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1,
                  treeHeight=0.5,
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
)
results[[1]]

#### Cellular Component ####
goDivision="CC"

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Strawberry/perl/bin/perl.exe", 
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25,
)

results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=-log(0.1,10),
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5,
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
)
results[[1]]