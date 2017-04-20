# SCOUP

SCOUP : a probabilistic model based on the Ornstein-Uhlenbeck process to analyze single-cell expression data during differentiation.

## Reference

https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1109-3

## Requirements

The following two libraries are necessary for pseudo-time estimation based on the shortest path on the PCA space.
** This pseudo-time is only used for initialing SCOUP, and hence, pseudo-time estimates from other methods or experimental time can be substituted for initialization.**

* LAPACK
* BLAS

## How to build

```
git clone https://github.com/hmatsu1226/SCOUP
cd SCOUP
make
```
Or download from "Download ZIP" button and unzip it.

## Running SP
Estimate pseudo-time based on shortest path on the PCA space.

##### Usage
```
./sp <Input_file1> <Input_file2> <Output_file1> <Output_file2> <G> <C> <D>
```

* Input_file1 : G x C matrix of expression data
* Input_file2 : Initial distribution data
* Output_file1 : Pseudo-time estimates
* Output_file2 : Coordinates of PCA
* G : The number of genes
* C : The number of cells
* D : The number of PCA dimensions

##### Format of Input_file1
The Input_file1 is the G x C matrix of expression data (separated with 'TAB').
Each row corresponds to each gene, and each column corresponds to each cell.

##### Example of Input_file1
```
0.33	-4.95	-1.37	-4.07	...
5.01	4.45	3.82	3.02	...
.
.
.
```

##### Format of Input_file2
The Input_file2 contains the mean and variance of the initial normal distribution.

* Col1 : Index of a gene (0-origin)
* Col2 : Mean of the initial distribution for a gene
* Col3 : Variance of the initial distribution for a gene

##### Example of Input_file2
```
0	0.0	1.7
1	1.0	2.3
2	-2.0	5.9
```

##### Format of Output_file1
The Output_file1 contains the pseudo-time estimates.

* Col1 : Index of a cell (0-origin)
* Col2 : Pseudo-time of a cell

##### Example of Output_file1
```
0	0.826988
1	0.102140
2	0.758120
```

##### Format of Output_file2
The Output_file2 contains the coordinates of PCA.

* Col1 : Index of a cell (0-origin)
* Col2 - Col(D+1) : Coordinates of a cell

This file contain (C+1) lines and the last line corresponds to the root cell defined by the mean of the initial distribution.

##### Example of Output_file2
```
0	3.04	0.42	
1	-21.21	-1.52	
2	5.76	0.48
```


## Running SCOUP
Estimate the parameters of Mixute Ornstein-Uhlenbeck process.

##### Usage
```
./scoup <Options> <Input_file1> <Input_file2> <Input_file3> <Output_file1> <Output_file2> <Output_file3> <G> <C>
```

* Input_file1 : G x C matrix of expression data
* Input_file2 : Initial distribution data
* Input_file3 : Initial pseudo-time data
* Output_file1 : Optimized parameters related to genes and lineages
* Output_file2 : Optimized parameters related to cells
* Output_file3 : Log-likelihood
* G : The number of genes
* C : The number of cells

##### Options

* -k INT : The number of lineages (default is 1)
* -m INT : Upper bound of EM iteration (without pseudo-time optimization). The detailed explanation is described in the supplementary text. (default is 1,000)
* -M INT : Upper bound of EM iteration (including pseudo-time optimization) (default is 10,000).
* -a DOUBLE : Lower bound of alpha (default is 0.1)
* -A DOUBLE : Upper bound of alpha (default is 100)
* -t DOUBLE : Lower bound of pseudo-time (default is 0.001)
* -T DOUBLE : Upper bound of pseudo-time (default is 2.0)
* -s DOUBLE : Lower bound of sigma squared (default is 0.1)

##### Example of running SCOUP
```
./scoup -k 2 data/data.txt data/init.txt out/time_sp.txt out/gpara.txt out/cpara.txt out/ll.txt 500 100
```

##### Format of Input_file1
This is the expression data matrix data and is the same data as the Input_file1 of SP.

##### Format of Input_file2
This is initial distribution and is the same data as the Input_file2 of SP.

##### Format of Input_file3
This is the pseudo-time for initialization and is the same as the **Output_file1** of SP.

##### Format of Output_file1
The Output_file1 contains the optimized parameters related to genes and lineages.

* First line
	* Col1 and Col2 : Space
	* Col3 - Col(K+2) : The probability of each lineage (pi_k)
* After first line
	* Col1 : alpha_g
	* Col2 : sigma_g^2
	* Col3 - Col(K+1) : theta_{gk} 

##### Example of Output_file1
```
 	 	0.509804 	0.490196
0.501610	2.528400	-6.338714 	-2.273163
0.309094	13.046904	3.545862 	0.337260
0.223226	4.212808	-4.443503 	9.629989
2.707472	14.221109	3.959898 	-2.353994
4.361342	34.646044	1.392565 	0.789397
``` 

##### Format of Output_file2

* Col1 : Pseudo-time of a cell
* Col2 - Col(K) : Responsibility for each lineage


##### Example of Output_file2
```
0.941979	0.990196	0.009804	
2.000000	0.990196	0.009804	
2.000000	0.990196	0.009804	
1.102146	0.990196	0.009804	
0.839387	0.990196	0.009804
```

##### Format of Output_file3
The log-likelihood

##### Exapmle of Output_file3
```
```


## Running SCOUP from the middle of the activity
Re-estimate parameters from the middle of the activity.

##### Usage
```
./scoup_resume <Options> <Input_file1> <Input_file2> <Input_file3> <Input_file4> <Output_file1> <Output_file2> <Output_file3> <G> <C>
```

* Input_file1 : G x C matrix of expression data
* Input_file2 : Initial distribution data
* Input_file3 : ** Semi-optimized gene and lineage parameters (Output_file1 of scoup) **
* Input_file4 : ** Semi-optimized cell parameters (Output_file2 of scoup) **
* Output_file1 : Optimized parameters related to genes and lineages
* Output_file2 : Optimized parameters related to cells
* Output_file3 : Log-likelihood
* G : The number of genes
* C : The number of cells

##### Options
It is the same as the Options of "scoup".

##### Example of running SCOUP
```
./scoup_resume -k 2 -e 0.0001 data/data.txt data/init.txt out/gpara.txt out/cpara.txt out/gpara_2.txt out/cpara_2.txt out/ll_2.txt 500 100
```

##### Format of Input_file1
This is the same as the Input_file1 of "scoup".

##### Format of Input_file2
This is the same as the Input_file2 of "scoup".

##### Format of Input_file3
This is the parameters related to genes and lineages and is the same as the **Output_file1** of SCOUP.

##### Format of Input_file4
This is the parameters related to cells and is the same as the **Output_file2** of "scoup".

##### Format of Output_file1, 2, 3
These file are the same as the output files of SCOUP.


## Running Correlation analysis
Calculate the correlation between genes after standardization.

##### Usage
```
./cor <Options> <Input_file1> <Input_file2> <Input_file3> <Input_file4> <Output_file1> <Output_file2> <G> <C>
```

* Input_file1 : G x C matrix of expression data
* Input_file2 : Initial distribution data
* Input_file3 : Optimized gene and lineage parameters (Output_file1 of scoup)
* Input_file4 : Optimized cell parameters (Output_file2 of scoup)
* Output_file1 : Standardized expression matrix
* Output_file2 : G x G correlation matrix
* G : The number of genes
* C : The number of cells

##### Options

##### Example of running Correlation analysis
```
./cor data/data.txt data/init.txt out/gpara.txt out/cpara.txt out/nexp.txt out/cor.txt 500 100
```

##### Format of Output_file1
The Output_file1 contains the standardized expression data.

##### Format of Output_file2
The Output_file2 contains the correlation for the standardized expression data.

## License
Copyright (c) 2015 Hirotaka Matsumoto
Released under the MIT license
