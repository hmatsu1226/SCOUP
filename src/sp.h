#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <string>
#include <string.h>
#include <algorithm>
#include <cassert>
#include <map>
#include "node.h"

using namespace std;

extern "C" int dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);

class Pseudo_Time{
public:
	int _gene_num;
	int _cell_num;
	int _dim;

	Pseudo_Time(int g, int c, int dim){
		_gene_num = g;
		_cell_num = c;
		_dim = dim;
		genes.resize(_gene_num);
		cells.resize(_cell_num);

		for(int i=0; i<_gene_num; i++){
			genes[i].Init(i);
		}
		for(int i=0; i<_cell_num; i++){
			cells[i].Init(_gene_num);
		}
	}

	double Dist(vector<vector<double> > &pos, int id1, int id2, int dim){
		double ret = 0.0;
		for(int i=0; i<dim; i++){
			ret += (pos[id1][i] - pos[id2][i]) * (pos[id1][i] - pos[id2][i]);
		}
		return sqrtf(ret);
	}

	void Prim(FILE *ftime, FILE *fpca){
		//run PCA
		int data_num = _cell_num + 1;
		//normalization todo variance
		vector<vector<double> > normalized_data(data_num, vector<double>(_gene_num, 0));
		vector<double> ave(_gene_num, 0);
		vector<double> var(_gene_num, 0);
		//average
		for(int i=0; i<_gene_num; i++){
			for(int j=0; j<_cell_num; j++){
				ave[i] += cells[j].Get_expression(i);
			}
			ave[i] /= _cell_num;
		}
		//variance
		for(int i=0; i<_gene_num; i++){
			for(int j=0; j<_cell_num; j++){
				var[i] += (cells[j].Get_expression(i) - ave[i]) * (cells[j].Get_expression(i) - ave[i]);
			}
			var[i] /= _cell_num;
		}
		//normalization
		for(int i=0; i<_cell_num; i++){
			for(int j=0; j<_gene_num; j++){
				if(var[j] != 0){
					normalized_data[i][j] = (cells[i].Get_expression(j) - ave[j])/sqrtf(var[j]);
				}
				else{
					normalized_data[i][j] = (cells[i].Get_expression(j) - ave[j]);
				}
			}
		}
		//add root cell
		for(int i=0; i<_gene_num; i++){
			if(var[i] != 0){
				normalized_data[data_num-1][i] = (genes[i].Initial_expression() - ave[i])/sqrtf(var[i]);
			}
			else{
				normalized_data[data_num-1][i] = (genes[i].Initial_expression() - ave[i]);
			}
		}

		//calculate variance-covariance matrix
		double tmp;
		vector<double> var_cov_matrix(_gene_num*_gene_num, 0);
		for(int i=0; i<_gene_num-1; i++){
			for(int j=i; j<_gene_num; j++){
				for(int k=0; k<data_num; k++){
					tmp = normalized_data[k][i] * normalized_data[k][j];

					if(i == j){
						var_cov_matrix[i*_gene_num + j] += tmp;
					}
					else{
						var_cov_matrix[i*_gene_num + j] += tmp;
						var_cov_matrix[j*_gene_num + i] += tmp;
					}
				}
			}
		}

		for(int i=0; i<_gene_num; i++){
			for(int j=0; j<_gene_num; j++){
				var_cov_matrix[i*_gene_num + j] /= (double)(data_num-1);
			}
		}

		//PCA
		int info=0, lwork = 3*_gene_num;
		vector<double> w(_gene_num, 0);
		vector<double> work(lwork, 0);
		char jobz = 'V', uplo = 'U';
		dsyev_(&jobz, &uplo, &_gene_num, &var_cov_matrix[0], &_gene_num, &w[0], &work[0], &lwork, &info);

		vector<vector<double> > z(data_num, vector<double>(_dim, 0));
		for (int i=0; i<data_num; i++){
			for(int j=0; j<_dim; j++){
				//todo 
				z[i][j] = 0;
				for(int k=0; k<_gene_num; k++){
					z[i][j] += normalized_data[i][k] * var_cov_matrix[(_gene_num-j-1)*_gene_num+k];
				}
			}
			//printf("%d\t%lf\t%lf\n", i, z[i][0], z[i][1]);
		}


		//calculate distance on PCA
		vector<vector<double> > edges_cost(data_num, vector<double>(data_num,0));
		for(int i=0; i<data_num; i++){
			for(int j=i+1; j<data_num; j++){
				edges_cost[i][j] = Dist(z, i, j, _dim);
				edges_cost[j][i] = edges_cost[i][j];
			}
		}

		//prim
		int min_node_from, min_node_to, erase_id;
		double min, max_pseudo_time = 0;
		vector<vector<double> > mst(data_num, vector<double>(data_num,0));
		vector<int> checked_node;
		vector<int> uncheked_node(data_num-1);
		vector<double> pseudo_time(data_num, 0);
		for(int i=0; i<data_num-1; i++){
			uncheked_node[i] = i;
		}
		checked_node.push_back(data_num-1);
		for(int i=0; i<data_num-1; i++){
			min = DBL_MAX;
			for(int j=0; j<checked_node.size(); j++){
				for(int k=0; k<uncheked_node.size(); k++){
					if(min > edges_cost[checked_node[j]][uncheked_node[k]]){
						min = edges_cost[checked_node[j]][uncheked_node[k]];
						min_node_from = checked_node[j];
						min_node_to = uncheked_node[k];
						erase_id = k;
					}
				}
			}
			mst[min_node_from][min_node_to] = 1;
			checked_node.push_back(min_node_to);
			uncheked_node.erase(uncheked_node.begin() + erase_id);

			//
			pseudo_time[min_node_to] = pseudo_time[min_node_from] + min;
			if(max_pseudo_time < pseudo_time[min_node_to]){
				max_pseudo_time = pseudo_time[min_node_to];
			}
		}

		//normalize pseudo time
		for(int i=0; i<_cell_num; i++){
			pseudo_time[i] /= max_pseudo_time;
			cells[i].Add_time(pseudo_time[i]);
			fprintf(ftime, "%d\t%lf\n", i, pseudo_time[i]);
		}

		for(int i=0; i<data_num; i++){
			fprintf(fpca, "%d\t", i);
			for(int j=0; j<_dim; j++){
				fprintf(fpca, "%lf\t", z[i][j]);
			}
			fprintf(fpca, "\n");
		}
	}

	int Set_initial_parameter(FILE *fp){
		int count, id;
		double expression, dispersion;
		for(int i=0; i<_gene_num; i++){
			count = fscanf(fp, "%d\t%lf\t%lf\n", &id, &expression, &dispersion);
			if(count == EOF){
				printf("error at reading initial parameter\n");
				return 1;
			}

			genes[i].Add_initial_expression(expression);
			genes[i].Add_initial_dispersion(dispersion);
		}
		return 0;
	}

	//read expression data
	void Set_expression(FILE *fp){
		char buf[1000000];
		char *tmp;
		
		double tmp_expression;
		for(int i=0; i<_gene_num; i++){
			fgets(buf, 1000000, fp);	
			
			tmp = strtok(buf, "\t");
			if(strcmp(tmp, "NA") == 0){
				cells[0].Add_expression(i, 0, 1);
			}
			else{
				tmp_expression = atof(tmp);
				cells[0].Add_expression(i, tmp_expression, 0);
			}

			for(int j=1; j<_cell_num-1; j++){
				tmp = strtok(NULL, "\t");
				if(strcmp(tmp, "NA") == 0){
					cells[j].Add_expression(i, 0, 1);
				}
				else{
					tmp_expression = atof(tmp);
					cells[j].Add_expression(i, tmp_expression, 0);
				}
			}

			tmp = strtok(NULL, "\n");
			if(strcmp(tmp, "NA") == 0){
				cells[_cell_num-1].Add_expression(i, 0, 1);
			}
			else{
				tmp_expression = atof(tmp);
				cells[_cell_num-1].Add_expression(i, tmp_expression, 0);
			}			
		}

		return;
	}
};

