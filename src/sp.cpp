#include <cstdio>
#include <vector>
#include <string>
#include <map>
#include "sp.h"

using namespace std;

int main(int argc, char* argv[]){
	FILE *fin_expression, *fin_initial, *fout_pseudo_time, *fout_mst, *fout_pca;
	if((fin_expression=fopen(argv[1], "r")) == NULL){
		printf("cannot open expression data\n");
		return 1;
	}
	if((fin_initial=fopen(argv[2], "r")) == NULL){
		printf("cannot open initial distribution data\n");
		return 1;
	}
	if((fout_pseudo_time=fopen(argv[3], "w")) == NULL){
		printf("cannot open Output_file1 (pseudo-time)\n");
		return 1;
	}
	if((fout_pca=fopen(argv[4], "w")) == NULL){
		printf("cannot open Output_file2 (PCA)\n");
		return 1;
	}

	int gene_num = atoi(argv[5]);
	int cell_num = atoi(argv[6]);
	int dim = atoi(argv[7]);

	Pseudo_Time PT(gene_num, cell_num, dim);
	PT.Set_expression(fin_expression);
	PT.Set_initial_parameter(fin_initial);

	//Prim
	PT.Prim(fout_pseudo_time, fout_pca);
	
	return 0;
}