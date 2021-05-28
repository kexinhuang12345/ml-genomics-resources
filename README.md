# Machine Learning for Genomics and Therapeutics Resources

This repo accompanies our survey paper: 

[Machine Learning Applications for Therapeutic Tasks with Genomics Data](https://arxiv.org/abs/2105.01171).
*Kexin Huang, Cao Xiao, Lucas M. Glass, Cathy W. Critchlow, Greg Gibson, Jimeng Sun*

We list tools, algorithms, data for this area. Feel free to make a pull request for new resources.

---- 


## Machine Learning for Genomics in Target Discovery

### Theme 1: Facilitating Understanding of Human Biology

#### Task 1: DNA-protein and RNA-protein binding prediction

**Task Description** Given a set of DNA/RNA sequences predict their binding scores. After training,use feature importance attribution methods to identify the motifs. 

67, 68, 66, 69, 70, 41, 71


#### Task 2: Methylation state prediction

**Task Description** For a DNA/RNA position with missing methylation status, given its availableneighboring methylation states and the DNA/RNA sequence, predict the methylation status on the positionof interest. 

72-77

#### Task 3: RNA splicing prediction

**Task Description** Given an RNA sequence and its cell type, if available, for each nucleotide,predicts the probability of being a spliced breakpoint and the splicing level. 

78-81

#### Task 4: Spatial gene expression inference

**Task Description** Given the histopathology image of the tissue, predict the gene expression forevery gene at each spatial transcriptomics spot.

82-83, 62

#### Task 5: Cell composition analysis

**Task Description** Given the gene expressions of a set of cells (in bulk RNA-seq or a spot in spatialtranscriptomics), infer proportion estimates of each cell type for this set.

84-89

#### Task 6: Gene network construction

**Task Description** Given a set of gene expression profiles of a gene set, identify the gene regulatorynetwork by predicting all pairs of interacting genes. 

90-92

### Theme 2: Identifying Druggable Biomarkers

#### Task 1: Variant calling

**Task Description** Given the aligned sequencing data ((1) read pileup image, which is a matrix ofdimension M and N, with M the number of reads and N the length of reads; or (2) the raw reads, which are aset of sequences strings) for each locus, classify the multi-class variant status.

94-98

#### Task 2: Variant pathogenicity prioritization

**Task Description** Given features about a variant, predict its corresponding disease risk and thenrank all variants based on the disease risk. Alternatively, given the DNA sequence or other related genomicsfeatures, predict the likelihood of disease risk for this sequence and retrieve the variant in the sequence thatcontributes highly to the risk prediction. 

99-104

#### Task 3: Rare disease detection
**Task Description** Given the gene expression data and other auxiliary data of a patient predictwhether this patient has a rare disease. Also, identify genetic variants for this rare disease

105-108

#### Task 4: Gene-disease association prediction
**Task Description** Given the known gene-disease association network and auxiliary information,predict the association likelihood for every unknown gene-disease pair.

109-114

#### Task 5: Pathway analysis and prediction
**Task Description** 

115-121

## Machine Learning for Genomics in Therapeutics Discovery
### Theme 1: Towards Precision Medicine

#### Task 1: Drug Response Prediction
**Task Description** 

123-130

#### Task 2: Drug Combination Therapy Prediction

**Task Description** 

131-136

### Theme 2: Improving Efficacy and Delivery of Gene Therapy
#### Task 1: CRISPR on-target outcome prediction

**Task Description** 

137-144

#### Task 2: CRISPR off-target prediction
**Task Description** 

145-152

#### Task 3: Virus vector design
**Task Description** 

153-155

## Machine Learning for Genomics in Clinical Study
**Task Description** 

### Theme 1: Translating Preclinical Animal Models to Humans

#### Task 1: Cross-species genotype-phenotype translation
**Task Description** 

157-162

### Theme 2: Curating High-quality Cohort

#### Task 1: Patient stratification/disease sub-typing
**Task Description** 

165-173

#### Task 2: Matching patients for genome-driven trials
**Task Description** 

174-179

### Theme 3: Inferring Causal Effects

#### Task 1: Mendelian randomization
**Task Description**

180-184

## Machine Learning for Genomics in Post-Market Study
### Theme 1: Mining Real-World Evidence

#### Task 1: Mining genomics-related markers from clinical free texts
**Task Description** 

185-189

#### Task 2: Discovering drug-gene/disease-gene interactions from scientific literature
**Task Description**
190-199
