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

[Jian Zhou and Olga G Troyanskaya. Predicting effects of noncoding variants with deep learning–based
sequence model. Nature Methods, 12(10):931–934, 2015.](https://pubmed.ncbi.nlm.nih.gov/26301843/)

[Babak Alipanahi, Andrew Delong, Matthew T Weirauch, and Brendan J Frey. Predicting the sequence specificities of dna-and rna-binding proteins by deep learning. Nature Biotechnology, 33(8):831–838, 2015.](https://www.nature.com/articles/nbt.3300)

[Martin Kircher, Daniela M Witten, Preti Jain, Brian J O’Roak, Gregory M Cooper, and Jay Shendure. A
general framework for estimating the relative pathogenicity of human genetic variants. Nature Genetics,
46(3):310–315, 2014.](https://pubmed.ncbi.nlm.nih.gov/24487276/)

[69] David R Kelley, Jasper Snoek, and John L Rinn. Basset: learning the regulatory code of the accessible
genome with deep convolutional neural networks. Genome Research, 26(7):990–999, 2016.

[70] Qinhu Zhang, Lin Zhu, and De-Shuang Huang. High-order convolutional neural network architecture
for predicting dna-protein binding sites. IEEE/ACM transactions on Computational Biology and
Bioinformatics, 16(4):1184–1192, 2018.

[71] Zhen Cao and Shihua Zhang. Simple tricks of convolutional neural network architectures improve
dna–protein binding prediction. Bioinformatics, 35(11):1837–1843, 2019.


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
**Task Description** Given the gene expression data for a phenotype and known gene relations, identify a set of genes corresponding to disease pathways.

115-121

## Machine Learning for Genomics in Therapeutics Discovery

### Theme 1: Improving Context-specific Drug Response

#### Task 1: Drug Response Prediction
**Task Description** Given a pair of drug compound molecular structure and gene expression profile of the cell line, predict the drug response in this context

123-130

#### Task 2: Drug Combination Therapy Prediction

**Task Description** Given a combination of drug compound structures and a cell line’s genomics profile, predict the combination response.

131-136

### Theme 2: Improving Efficacy and Delivery of Gene Therapy
#### Task 1: CRISPR on-target outcome prediction

**Task Description** With a fixed target, given the gRNA sequence and other auxiliary information such as target gene expression and epigenetic profile, predict its on-target repair outcome.

137-144

#### Task 2: CRISPR off-target prediction
**Task Description** Given the gRNA sequence and the off-target DNA sequence, predict its off-target effect.

145-152

#### Task 3: Virus vector design
**Task Description** Given a set of virus sequences and their labels for a property X, obtain an accurate predictor oracle and conduct various generation modeling to generate de novo virus variants with a high score in X and high diversity.

153-155

## Machine Learning for Genomics in Clinical Study
**Task Description** 

### Theme 1: Translating Preclinical Animal Models to Humans

#### Task 1: Cross-species genotype-phenotype translation
**Task Description** : Given genotype-phenotype data of animals and only the genotype data of humans, train the model to fit phenotype from the genotype and transfer this model to human. 

157-162

### Theme 2: Curating High-quality Cohort

#### Task 1: Patient stratification/disease sub-typing
**Task Description** Given the gene expression and other auxiliary information for a set of patients produce criteria for patient stratification.

165-173

#### Task 2: Matching patients for genome-driven trials
**Task Description** Given a pair of patient data (genomics, EHR, etc.) and trial eligibility criteria (text description), predict the matching likelihood.

174-179

### Theme 3: Inferring Causal Effects

#### Task 1: Mendelian randomization
**Task Description** Given observation data of the genomic factor, exposure, outcome, and other auxiliary information formulate or identify the causal relations among them and compute the effect of the exposure to the outcome.

180-184

## Machine Learning for Genomics in Post-Market Study
### Theme 1: Mining Real-World Evidence

#### Task 1: Mining genomics-related markers from clinical free texts
**Task Description**  Given a clinical note document, predict the genomic biomarker variable of interest.

185-189

#### Task 2: Discovering drug-gene/disease-gene interactions from scientific literature
**Task Description** Given a document from literature, extract the drug-gene, drug-disease terms, and predict the interaction types from the text.

190-199
