#include <unistd.h>
#include <cstdio>
#include <vector>
#include <string>
#include <map>
#include "ou.h"

using namespace std;

int main(int argc, char* argv[]){
	int K = 1;

	int ch;
	extern char *optarg;
	extern int optind, opterr;
	while((ch = getopt(argc, argv, "k:")) != -1){
		switch(ch){
		case 'k':
			K = atoi(optarg);
			break;	
		case ':':
			printf("invalid option\n");
			return 1;
		case '?':
			printf("invalid option\n");
			return 1;
		}
	}

	FILE *fin_expression, *fin_init, *fin_gene_para, *fin_cell_para, *fout_nexp, *fout_cor;
	if((fin_expression=fopen(argv[1], "r")) == NULL){
		printf("cannot open expression data\n");
		return 1;
	}
	if((fin_init=fopen(argv[2], "r")) == NULL){
		printf("cannot open initial distribution data\n");
		return 1;
	}
	if((fin_gene_para=fopen(argv[3], "r")) == NULL){
		printf("cannot open gene and lineage parameters data\n");
		return 1;
	}
	if((fin_cell_para=fopen(argv[4], "r")) == NULL){
		printf("cannot open cell parameters data\n");
		return 1;
	}
	if((fout_nexp=fopen(argv[5], "w")) == NULL){
		printf("cannot open Output_file1 (normalized expression)\n");
		return 1;
	}
	if((fout_cor=fopen(argv[6], "w")) == NULL){
		printf("cannot open Output_file2 (correlation matrix)\n");
		return 1;
	}
		
	int gene_num = atoi(argv[7]);
	int cell_num = atoi(argv[8]);

	Continuous_OU_process OU(gene_num, cell_num, K);

	OU.Set_expression(fin_expression);
	if(OU.Set_initial_parameter(fin_init) == 1){
		return 1;
	}
	if(OU.Set_optimized_parameter(fin_gene_para, fin_cell_para) == 1){
		return 1;
	}

	//calc normalized expression
	OU.Print_correlation(fout_nexp, fout_cor);

	return 0;
}