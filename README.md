# MUTClass package
Java implementation of MUTClass algorithm for the identification of targeted panels of genes able to classify cancer patients.
The Java package includes the following 3 software:
- **DMGSFinder**: Extract from a mutation matrix a panel of mutated genes which maximizes (or minimizes) the differential coverage between two classes of samples.
- **MUTClassCV**: Run a cross-validation test of MUTClass algorithm on a mutation matrix 
- **MUTClass**: Train MUTClass on a mutation matrix A and test the algorithm on a mutation matrix B in order to classify samples of B.
<br />

## DMGSFinder
DMGSFinder is a greedy algorithm for maximizing (or minimizing) the differential coverage with respect to a set of positive samples and a set of negative samples.
DMGSFinder is the algorithm used by MUTClass to solve the k-MaxDiffCov (or k-MinDiffCov) problem

<br />

**Usage**:

java -cp ./out DMGSFinder -m \<mutationsFile\> -d \<listGenes\> -k \<panelSize\> -min -o \<resultsFile\>
		
REQUIRED PARAMETERS:

&emsp;-m &emsp; Mutation matrix file

OPTIONAL PARAMETERS:

&emsp;-d &emsp;&emsp;List of driver genes

&emsp;-k &emsp;&emsp;Size of gene panel (default=10)

&emsp;-min &emsp;Solve the k-MinDiffCov problem (default=solve k-MaxDiffCov problem)

&emsp;-o &emsp;&emsp;Results file (default=print results to standard output)

<br />

**Files format**:

<br />

MUTATION MATRIX FILE

The first line contains the names of the samples.
If no list of driver genes is provided, the second line must contain sample classes.
The following lines contain the name of the gene followed by a list of numbers denoting the mutation frequency of the gene in each sample (0 means that a gene is mutated in the sample, any other value means that the gene is mutated in the sample).
Values, sample names (and sample classes, if provided) are separated by tabs (\t).

Example without sample classes:

Sample1&emsp;Sample2&emsp;Sample3&emsp;Sample4

Gene1&emsp;0&emsp;0&emsp;2&emsp;1

Gene2&emsp;1&emsp;0&emsp;3&emsp;0

Example with sample classes:

Sample1&emsp;Sample2&emsp;Sample3&emsp;Sample4

P&emsp;P&emsp;N&emsp;P

Gene1&emsp;0&emsp;0&emsp;2&emsp;1

Gene2&emsp;1&emsp;0&emsp;3&emsp;0

<br />

RESULTS FILE

The results file contains the panel found and the following statistics separated by tab character (\t): 
1. Average coverage in the set of positive samples;	
2. Average coverage in the set of negative samples;
3. Differential coverage.

Example:

DMGS&emsp;Average positive coverage&emsp;Average negative coverage&emsp;Differential coverage

[DOCK11, LOC100506083, TMEM138, IL1R2, CASP4, IL31RA, PLEKHA6, SPATA31C1, FER1L6, TANC2]&emsp;74.39%&emsp;19.0%&emsp;55.39%

<br />

**Examples**:

1. Find a panel with 10 genes which maximizes the differential coverage, starting from a set of driver genes D = {BRCA1,BRCA2}. The positive set will be the set of samples where at least one gene among BRCA1 and BRCA2 is mutated, the negative set will be the one in which neither BRCA1 nor BRCA2 are mutated.


java -cp ./out DMGSFinder -m Data/BRCA_snp_gene_matrix.txt -d BRCA1,BRCA2


2. Find a panel with 20 genes which minimizes the differential coverage, starting from a mutation matrix file with sample classes and save final results to "results.txt"


java -cp ./out DMGSFinder -m Data/BRCA_snp_gene_matrix_with_classes.txt -min -k 20 -o results.txt

<br />

## MUTClassCV

<br />

## MUTClass
