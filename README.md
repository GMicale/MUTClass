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

java -cp ./out DMGSFinder -m \<mutationsFile\> -d \<listDriverGenes\> -k \<panelSize\> -min -o \<resultsFile\>
		
REQUIRED PARAMETERS:

&emsp;-m &emsp; Mutation matrix file

OPTIONAL PARAMETERS:

&emsp;-d &emsp;&emsp;List of driver genes (default=matrix file with classes)

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

MUTClassCV performs a cross-validation test of MUTClass algorithm on a mutation matrix

<br />

**Usage**:

java -cp ./out MUTClassCV -m \<mutationsFile\> -d \<listDriverGenes\> -kmax \<positivePanelSize\> -kmin \<negativePanelSize\> -cv \<crossValidationIterations\> -f \<crossValidationFolds\>
		
REQUIRED PARAMETERS:

&emsp;-m &emsp; Mutation matrix file

OPTIONAL PARAMETERS:

&emsp;-d &emsp;&emsp;List of driver genes (default=matrix file with classes)

&emsp;-kmax &emsp;&emsp;Size of positive gene panel (default=10)

&emsp;-kmin &emsp;&emsp;Size of negative gene panel (default=10)

&emsp;-cv &emsp;Number of iterations of cross validation (default=10)

&emsp;-f &emsp;&emsp;Number of folds for cross validation (default=5)

<br />

**Files format**:

<br />

MUTATION MATRIX FILE

Same mutation matrix file format described for DMGSFinder. If no list of driver genes is provided, sample classes must be provided in the mutation matrix file. 

<br />

**Output**:

A list of performance statistics returned by test cross-validation test, including true positives, true negatives, false positives, false negatives, unclassified samples, precision, recall, FPR, FNR, specificity, accuracy and F1 score.

<br />

**Example**:

Run 5-fold 10 times repeated cross-validation on BRCA mutation matrix with sample classes, with kmax=5 and kmin=50.

java -cp ./out MUTClassCV -m Data/BRCA_snp_gene_matrix_with_classes.txt -kmax 5 -kmin 50 -cv 10 -f 5

<br />


## MUTClass

Train MUTClass on a mutation matrix A and test the algorithm on a mutation matrix B in order to classify samples of B.

<br />

**Usage**:

java -cp ./out MUTClass -train \<trainingSetFile\> -test \<testSetFile\> -d \<listDriverGenes\> -kmax \<positivePanelSize\> -kmin \<negativePanelSize\> -o \<resultsFile\>
		
REQUIRED PARAMETERS:

&emsp;-train &emsp; Training mutation matrix file

&emsp;-test &emsp; Test mutation matrix file

OPTIONAL PARAMETERS:

&emsp;-d &emsp;&emsp;List of driver genes (default=training matrix file with classes)

&emsp;-kmax &emsp;&emsp;Size of positive gene panel (default=5)

&emsp;-kmin &emsp;&emsp;Size of negative gene panel (default=50)

&emsp;-o &emsp;Output result file (default='results.txt')

<br />

**Files format**:

<br />

TRAINING MUTATION MATRIX FILE

Same mutation matrix file format described for DMGSFinder. If no list of driver genes is provided, sample classes must be provided in the mutation matrix file. 

<br />

TEST MUTATION MATRIX FILE

A mutation matrix file with no sample classes in the same format described for DMGSFinder.

<br />

OUTPUT RESULT FILE

A text file listing the positive panel, the negative panel and the predicted classes for each sample in the test mutation matrix.

Example:

Positive panel:

[ACAB1, BRAF, DOCK11, LOC100506083, TMEM138]

Negative panel:

[IL1R2, CASP4, IL31RA, PLEKHA6, SPATA31C1]

Predicted classes:

Sample1&emsp;P

Sample2&emsp;N

Sample3&emsp;N

Sample4&emsp;P

Sample5&emsp;P

<br />

**Example**:

Train MUTClass with kmax=10 and kmin=20 on BRCA mutation matrix with sample classes and test it on PCAWG-BRCA mutation matrix.

java -cp ./out MUTClass -train Data/BRCA_snp_gene_matrix_with_classes.txt -test Data/PCAWG-BRCA_gene_matrix.txt -kmax 10 -kmin 20

<br />
