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
#include <cfloat>
#include <map>
#include "node.h"

using namespace std;

class Continuous_OU_process{
public:
	int _gene_num;
	int _cell_num;
	int _K;

	double _max_time;
	double _min_time;
	double _min_alpha;
	double _max_alpha;
	double _min_sigma_squared;

	double _init_alpha;
	double _init_siqma_squared;

	double _old_ll;

	int _max_ite1;
	int _max_ite2;
	double _thresh;

	Continuous_OU_process(int g, int c, int k, int max_ite1, int max_ite2, double alpha_min, double alpha_max, double t_min, double t_max, double sigma_squared_min, double thresh){
		_max_ite1 = max_ite1;
		_max_ite2 = max_ite2;
		_init_alpha = 5.0;
		_init_siqma_squared = 1.0;
		_old_ll = -DBL_MIN;
		_min_time = t_min;
		_max_time = t_max;
		_min_alpha = alpha_min;
		_max_alpha = alpha_max;
		_min_sigma_squared = sigma_squared_min;
		_thresh = thresh;

		_gene_num = g;
		_cell_num = c;
		_K = k;
		genes.resize(_gene_num);
		cells.resize(_cell_num);
		lineages.resize(_K);

		for(int i=0; i<_gene_num; i++){
			genes[i].Init(_K, i, _init_alpha, _init_siqma_squared);
		}
		for(int i=0; i<_cell_num; i++){
			cells[i].Init(_gene_num, _K);
		}
		for(int i=0; i<_K; i++){
			lineages[i].Init(_gene_num, 1.0/_K);
		}
	}

	Continuous_OU_process(int g, int c, int k){
		_gene_num = g;
		_cell_num = c;
		_K = k;
		genes.resize(_gene_num);
		cells.resize(_cell_num);
		lineages.resize(_K);

		for(int i=0; i<_gene_num; i++){
			genes[i].Init(_K, i, _init_alpha, _init_siqma_squared);
		}
		for(int i=0; i<_cell_num; i++){
			cells[i].Init(_gene_num, _K);
		}
		for(int i=0; i<_K; i++){
			lineages[i].Init(_gene_num, 1.0/_K);
		}
	}

	void Init_EM(){
		//Initialize with random responsibility
		Set_parameter();
		M_step();
		//update parameters
		for(int j=0; j<_K; j++){
			lineages[j].Update_parameter();
		}
		for(int j=0; j<_gene_num; j++){
			genes[j].Update_parameter();
		}
	}

	int EM(){
		double ll=0;

		int step=100, id=1;
		for(int i=0; i<_max_ite1; i++){
			E_step();
			M_step();

			//update parameters
			for(int j=0; j<_K; j++){
				lineages[j].Update_parameter();
			}
			for(int j=0; j<_gene_num; j++){
				genes[j].Update_parameter();
			}

			ll = Log_likelihood();

			//debug
			if(i%step == 0){
				printf("%d-th iteration in first EM\nlog-likelihood: %lf\n", i, ll);
			}

			_old_ll = ll;
		}

		for(int i=0; i<_max_ite2; i++){
			E_step();
			Optimize_time();
			M_step();

			//check parameter convergenece
			if(Convergence() == 1){
				break;
			}

			if(i == _max_ite2-1){
				break;
			}

			//update parameters
			for(int j=0; j<_K; j++){
				lineages[j].Update_parameter();
			}
			for(int j=0; j<_gene_num; j++){
				genes[j].Update_parameter();
			}

			ll = Log_likelihood();

			//debug
			if(i%step == 0){
				printf("%d-th iteration in second EM\nlog-likelihood: %lf\n", i, ll);
			}

			_old_ll = ll;

		}

		return 0;
	}

	void E_step(){
		//calculate the responsibility of the latent value Z
		for(int i=0; i<_cell_num; i++){
			cells[i].Calc_responsibility();
		}
	}

	void M_step(){
		//optimize theta
		double new_theta;
		long double numerator = 0.0;
		long double denominator = 0.0;
		long double e_minus_at_power;
		
		for(int i=0; i<_gene_num; i++){
			for(int j=0; j<_K; j++){
				new_theta = lineages[j].Theta(i);
				numerator = 0.0;
				denominator = 0.0;
				for(int k=0; k<_cell_num; k++){
					e_minus_at_power = exp(-genes[i].Alpha()*cells[k].Time());
					numerator += cells[k].Gamma(j) * 2 * (cells[k].Xn(i) - e_minus_at_power*X0(i,k,j) - (1 - e_minus_at_power)*lineages[j].Theta(i)) / (genes[i].Alpha()*(1+e_minus_at_power));
					//numerator += cells[k].Gamma(j) * ((cells[k].Xn(i)+genes[i].X0()-2*lineages[j].Theta(i))*(cosh(genes[i].Alpha()*cells[k].Time())-1)/sinh(genes[i].Alpha()*cells[k].Time()) + cells[k].Xn(i) - genes[i].X0()) / genes[i].Alpha();
					denominator += cells[k].Gamma(j) * cells[k].Time();
				}
				new_theta += numerator/denominator;
				lineages[j].Add_new_theta(i, new_theta);
			}
		}

		//optimize alpha
		double new_alpha;
		double f_alpha;
		double at;
		for(int i=0; i<_gene_num; i++){
			numerator = 0.0;
			denominator = 0.0;
			for(int j=0; j<_cell_num; j++){
				for(int k=0; k<_K; k++){
					at = genes[i].Alpha() * cells[j].Time();
					f_alpha = genes[i].Sigma_squared()/(genes[i].Alpha()*genes[i].Alpha());
					f_alpha -= genes[i].Sigma_squared()/genes[i].Alpha() * cells[j].Time()*cosh(at)/sinh(at);
					f_alpha += 1.0/(genes[i].Alpha()*sinh(at)) * (-cosh(at) + at/sinh(at)) * (X0_minus_theta_squared(i,j,k) + Xn_minus_theta_squared(i,j,k));
					f_alpha += 2.0/(genes[i].Alpha()*sinh(at)) * (1 - at*cosh(at)/sinh(at)) * X0_minus_theta_times_Xn_minus_theta(i,j,k);

					numerator += cells[j].Gamma(k) * (-cells[j].Time()*genes[i].Sigma_squared() - X0_minus_theta_squared(i, j, k) + Xn_minus_theta_squared(i, j, k));
					denominator += cells[j].Gamma(k) * f_alpha;
				}
			}
			new_alpha = numerator/denominator;

			if(new_alpha < _min_alpha){
				new_alpha = _min_alpha;
			}
			else if(new_alpha > _max_alpha){
				new_alpha = _max_alpha;
			}

			genes[i].Add_new_alpha(new_alpha);
		}

		//optimize sigma_squared
		int effective_cell_num;
		double new_sigma_squared;
		for(int i=0; i<_gene_num; i++){
			effective_cell_num = 0;
			new_sigma_squared = 0.0;
			for(int j=0; j<_cell_num; j++){
				effective_cell_num++;

				at = genes[i].Alpha() * cells[j].Time();
				for(int k=0; k<_K; k++){
					new_sigma_squared += cells[j].Gamma(k) * 2 * genes[i].Alpha() / (1-exp(-2*at)) * (Xn_minus_theta_squared(i,j,k) - 2*exp(-at)*X0_minus_theta_times_Xn_minus_theta(i,j,k) + exp(-2*at)*X0_minus_theta_squared(i,j,k) );
				}
			}
			new_sigma_squared /= effective_cell_num;

			if(new_sigma_squared < _min_sigma_squared){
				new_sigma_squared = _min_sigma_squared;
			}

			genes[i].Add_new_sigma_squared(new_sigma_squared);
		}

		//optimize pi
		double new_pi, sum_pi=0;
		for(int i=0; i<_K; i++){
			new_pi = 0.0;
			for(int j=0; j<_cell_num; j++){
				new_pi += cells[j].Gamma(i);
			}
			lineages[i].Add_new_pi(new_pi);

			sum_pi += new_pi;
		}
		//normalization
		for(int i=0; i<_K; i++){
			lineages[i].Normalize_new_pi(sum_pi);
		}
	}

	void Optimize_time(){
		//optimize t
		int max_ite=100;
		double at, old_time, new_time, pre_time, E, f1, f2;
		for(int i=0; i<_cell_num; i++){
			new_time = cells[i].Time();
			pre_time = cells[i].Time();

			for(int ite=0; ite<max_ite; ite++){
				old_time = new_time;
				f1 = 0.0;
				f2 = 0.0;
				for(int j=0; j<_K; j++){
					for(int k=0; k<_gene_num; k++){
						at = genes[k].Alpha() * old_time;
						E = -0.5*at*(X0_minus_theta_squared(k,i,j,old_time) + Xn_minus_theta_squared(k,i,j));
						E += at*cosh(at)*X0_minus_theta_times_Xn_minus_theta(k,i,j,old_time);
						E /= sinh(at) * sinh(at);
						
						f1 += cells[i].Gamma(j) * (at*(1-cosh(at)/sinh(at)) - 2*genes[k].Alpha()/genes[k].Sigma_squared()*E);
					}
				}

				for(int j=0; j<_K; j++){
					for(int k=0; k<_gene_num; k++){
						at = genes[k].Alpha() * old_time;
						f2 += cells[i].Gamma(j) * genes[k].Alpha() * (1-cosh(at)/sinh(at));
						f2 += cells[i].Gamma(j) * at * genes[k].Alpha() / sinh(at) / sinh(at);
						f2 += cells[i].Gamma(j) * genes[k].Alpha() * genes[k].Alpha() * (sinh(at) - 2*at*cosh(at)) * (X0_minus_theta_squared(k,i,j,old_time) + Xn_minus_theta_squared(k,i,j)) / genes[k].Sigma_squared() / sinh(at) / sinh(at) / sinh(at);
						f2 -= 2 * cells[i].Gamma(j) * genes[k].Alpha() * genes[k].Alpha() * (sinh(at)*cosh(at) - at*(1+cosh(at)*cosh(at))) * X0_minus_theta_times_Xn_minus_theta(k,i,j,old_time) / genes[k].Sigma_squared() / sinh(at) / sinh(at) / sinh(at);
					}
				}

				new_time = old_time - f1/f2;

				if(new_time < _min_time){
					new_time = _max_time;
				}
				else if(new_time > _max_time){
					new_time = _min_time;
				}

				if(fabs(old_time-new_time) < _thresh){
					break;
				}
			}

			double ll0, ll1, ll2, ll3;
			ll0 = Log_likelihood_of_cell(i, new_time);
    		ll1 = Log_likelihood_of_cell(i, _min_time);
    		ll2 = Log_likelihood_of_cell(i, _max_time);
    		ll3 = Log_likelihood_of_cell(i, pre_time);

    		if(ll3 > ll0 && ll3 > ll1 && ll3 > ll2){
    			cells[i].Add_new_time(pre_time);
    		}
    		else if(ll0 >= ll1 && ll0 >= ll2){
				cells[i].Add_new_time(new_time);
			}
			else if(ll1 > ll2){
				cells[i].Add_new_time(_min_time);
			}
			else{
				cells[i].Add_new_time(_max_time);
			}
			cells[i].Update_parameter();
		}
	}

	void Calc_gene_responsibility(){
		for(int i=0; i<_cell_num; i++){
			cells[i].Calc_responsibility_of_gene();
		}
	}

	//for null modek K=1
	void EM_for_null_model(int gene_id){
		Set_null_parameter();
		int max_ite = 10;

		for(int i=0; i<max_ite; i++){
			M_step_for_null_model(gene_id);
			genes[gene_id].Update_null_parameter();
		}
	}

	void M_step_for_null_model(int gene_id){
		//optimize theta
		double new_theta;
		double numerator = 0.0;
		double denominator = 0.0;
		double e_minus_at_power;
		int i = gene_id;
		new_theta = genes[i].Theta_null();
		numerator = 0.0;
		denominator = 0.0;
		for(int j=0; j<_cell_num; j++){
			e_minus_at_power = exp(-genes[i].Alpha()*cells[j].Time());
			numerator += 2 * (cells[j].Xn(i) - e_minus_at_power*X0(i,j) - (1 - e_minus_at_power)*genes[i].Theta_null()) / (genes[i].Alpha()*(1+e_minus_at_power));
			denominator += cells[j].Time();
		}
		new_theta += numerator/denominator;
		genes[i].Add_new_theta_null(new_theta);

		//optimize alpha
		double new_alpha;
		double f_alpha;
		double at;
		numerator = 0.0;
		denominator = 0.0;
		for(int j=0; j<_cell_num; j++){
			at = genes[i].Alpha() * cells[j].Time();
			f_alpha = genes[i].Sigma_squared()/(genes[i].Alpha()*genes[i].Alpha());
			f_alpha -= genes[i].Sigma_squared()/genes[i].Alpha() * cells[j].Time()*cosh(at)/sinh(at);
			f_alpha += 1.0/(genes[i].Alpha()*sinh(at)) * (-cosh(at) + at/sinh(at)) * (X0_minus_theta_squared(i,j) + Xn_minus_theta_squared(i,j));
			f_alpha += 2.0/(genes[i].Alpha()*sinh(at)) * (1 - at*cosh(at)/sinh(at)) * X0_minus_theta_times_Xn_minus_theta(i,j);

			numerator += (-cells[j].Time()*genes[i].Sigma_squared() - X0_minus_theta_squared(i, j) + Xn_minus_theta_squared(i, j));
			denominator += f_alpha;
		}
		//with prior distribution
		//new_alpha = (numerator*_sigma_squared_alpha - 2*_mu_alpha*genes[i].Sigma_squared()) / (denominator*_sigma_squared_alpha - 2*genes[i].Sigma_squared());
		//without prior distribution
		new_alpha = numerator/denominator;

		if(isnan(new_alpha)){
			new_alpha = genes[i].Alpha();
		}
		else if(new_alpha < _min_alpha){
			new_alpha = _min_alpha;
		}
		else if(new_alpha > _max_alpha){
			new_alpha = _max_alpha;
		}
		genes[i].Add_new_alpha(new_alpha);

		//optimize sigma_squared
		int effective_cell_num;
		double new_sigma_squared;
		effective_cell_num = 0;
		new_sigma_squared = 0.0;
		for(int j=0; j<_cell_num; j++){
			effective_cell_num++;

			at = genes[i].Alpha() * cells[j].Time();
			new_sigma_squared += 2 * genes[i].Alpha() / (1-exp(-2*at)) * (Xn_minus_theta_squared(i,j) - 2*exp(-at)*X0_minus_theta_times_Xn_minus_theta(i,j) + exp(-2*at)*X0_minus_theta_squared(i,j) );
		}
		new_sigma_squared /= effective_cell_num;
		genes[i].Add_new_sigma_squared(new_sigma_squared);
	}

	void EM_for_null_model2(int gene_id){
		int max_ite = 100;

		genes[gene_id].Add_new_theta_null(genes[gene_id].X0());
		genes[gene_id].Update_null_parameter2();

		for(int i=0; i<max_ite; i++){
			M_step_for_null_model2(gene_id);
			genes[gene_id].Updata_null_parameter2_2();
		}
	}

	void M_step_for_null_model2(int gene_id){
		double numerator, denominator;
		//optimize alpha
		double new_alpha;
		double f_alpha;
		double at;
		int i = gene_id;
		numerator = 0.0;
		denominator = 0.0;
		for(int j=0; j<_cell_num; j++){
			at = genes[i].Alpha() * cells[j].Time();
			f_alpha = genes[i].Sigma_squared()/(genes[i].Alpha()*genes[i].Alpha());
			f_alpha -= genes[i].Sigma_squared()/genes[i].Alpha() * cells[j].Time()*cosh(at)/sinh(at);
			f_alpha += 1.0/(genes[i].Alpha()*sinh(at)) * (-cosh(at) + at/sinh(at)) * (X0_minus_theta_squared(i,j) + Xn_minus_theta_squared(i,j));
			f_alpha += 2.0/(genes[i].Alpha()*sinh(at)) * (1 - at*cosh(at)/sinh(at)) * X0_minus_theta_times_Xn_minus_theta(i,j);

			numerator += (-cells[j].Time()*genes[i].Sigma_squared() - X0_minus_theta_squared(i, j) + Xn_minus_theta_squared(i, j));
			denominator += f_alpha;
		}
		//with prior distribution
		//new_alpha = (numerator*_sigma_squared_alpha - 2*_mu_alpha*genes[i].Sigma_squared()) / (denominator*_sigma_squared_alpha - 2*genes[i].Sigma_squared());
		//without prior distribution
		new_alpha = numerator/denominator;

		//todo
		if(isnan(new_alpha)){
			new_alpha = genes[i].Alpha();
		}
		else if(new_alpha < _min_alpha){
			new_alpha = _min_alpha;
		}
		else if(new_alpha > _max_alpha){
			new_alpha = _max_alpha;
		}
		genes[i].Add_new_alpha(new_alpha);

		//optimize sigma_squared
		int effective_cell_num = 0;
		double new_sigma_squared = 0;
		for(int j=0; j<_cell_num; j++){
			effective_cell_num++;

			at = genes[i].Alpha() * cells[j].Time();
			new_sigma_squared += 2 * genes[i].Alpha() / (1-exp(-2*at)) * (Xn_minus_theta_squared(i,j) - 2*exp(-at)*X0_minus_theta_times_Xn_minus_theta(i,j) + exp(-2*at)*X0_minus_theta_squared(i,j) );
		}
		new_sigma_squared /= effective_cell_num;
		genes[i].Add_new_sigma_squared(new_sigma_squared);
	}

	double X0(int i, int j, int k){
		double at = genes[i].Alpha() * cells[j].Time();
		double var = 1.0/(2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1)) + 1.0/genes[i].Initial_dispersion());
		double mu;
		if(var < 0.000001){
			mu = exp(at)*cells[j].Xn(i)+(1-exp(at))*lineages[k].Theta(i);
		}
		else{
			mu = (2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1))*(exp(at)*cells[j].Xn(i)+(1-exp(at))*lineages[k].Theta(i)) + genes[i].Initial_expression()/genes[i].Initial_dispersion()) * var;
		}
		return mu;
	}

	//for null model
	double X0(int i, int j){
		double at = genes[i].Alpha() * cells[j].Time();
		double var = 1.0/(2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1)) + 1.0/genes[i].Initial_dispersion());
		double mu = (2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1))*(exp(at)*cells[j].Xn(i)+(1-exp(at))*genes[i].Theta_null()) + genes[i].Initial_expression()/genes[i].Initial_dispersion()) * var;
		return mu;
	}

	double X0_minus_mu_squared(int i, int j, int k){
		double at = genes[i].Alpha() * cells[j].Time();
		double var = 1.0/(2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1)) + 1.0/genes[i].Initial_dispersion());
		double mu;
		if(var < 0.000001){
			mu = exp(at)*cells[j].Xn(i)+(1-exp(at))*lineages[k].Theta(i);
		}
		else{
			mu = (2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1))*(exp(at)*cells[j].Xn(i)+(1-exp(at))*lineages[k].Theta(i)) + genes[i].Initial_expression()/genes[i].Initial_dispersion()) * var;
		}
		double tmp = mu - genes[i].Initial_expression();

		return (tmp * tmp) + var;
	}

	//for null model
	double X0_minus_mu_squared(int i, int j){
		double at = genes[i].Alpha() * cells[j].Time();
		double var = 1.0/(2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1)) + 1.0/genes[i].Initial_dispersion());
		double mu = (2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1))*(exp(at)*cells[j].Xn(i)+(1-exp(at))*genes[i].Theta_null()) + genes[i].Initial_expression()/genes[i].Initial_dispersion()) * var;
		double tmp = mu - genes[i].Initial_expression();

		return (tmp * tmp) + var;
	}

	double X0_minus_theta_squared(int i, int j, int k){
		double at = genes[i].Alpha() * cells[j].Time();
		double var = 1.0/(2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1)) + 1.0/genes[i].Initial_dispersion());
		double mu;
		if(var < 0.000001){
			mu = exp(at)*cells[j].Xn(i)+(1-exp(at))*lineages[k].Theta(i);
		}
		else{
			mu = (2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1))*(exp(at)*cells[j].Xn(i)+(1-exp(at))*lineages[k].Theta(i)) + genes[i].Initial_expression()/genes[i].Initial_dispersion()) * var;
		}		
		double tmp = mu - lineages[k].Theta(i);

		return (tmp * tmp) + var;
	}

	double X0_minus_theta_squared(int i, int j, int k, double t){
		double at = genes[i].Alpha() * t;
		double var = 1.0/(2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1)) + 1.0/genes[i].Initial_dispersion());
		double mu;
		if(var < 0.000001){
			mu = exp(at)*cells[j].Xn(i)+(1-exp(at))*lineages[k].Theta(i);
		}
		else{
			mu = (2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1))*(exp(at)*cells[j].Xn(i)+(1-exp(at))*lineages[k].Theta(i)) + genes[i].Initial_expression()/genes[i].Initial_dispersion()) * var;
		}
		double tmp = mu - lineages[k].Theta(i);

		return (tmp * tmp) + var;
	}

	//for null model
	double X0_minus_theta_squared(int i, int j){
		double at = genes[i].Alpha() * cells[j].Time();
		double var = 1.0/(2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1)) + 1.0/genes[i].Initial_dispersion());
		double mu = (2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1))*(exp(at)*cells[j].Xn(i)+(1-exp(at))*genes[i].Theta_null()) + genes[i].Initial_expression()/genes[i].Initial_dispersion()) * var;
		double tmp = mu - genes[i].Theta_null();

		return (tmp * tmp) + var;
	}

	double Xn_minus_theta_squared(int i, int j, int k){
		double tmp = cells[j].Xn(i) - lineages[k].Theta(i);
		return tmp * tmp;
	}

	//for null model
	double Xn_minus_theta_squared(int i, int j){
		double tmp = cells[j].Xn(i) - genes[i].Theta_null();
		return tmp * tmp;
	}

	double X0_minus_theta_times_Xn_minus_theta(int i, int j, int k){
		double at = genes[i].Alpha() * cells[j].Time();
		double var = 1.0/(2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1)) + 1.0/genes[i].Initial_dispersion());
		double mu;
		if(var < 0.000001){
			mu = exp(at)*cells[j].Xn(i)+(1-exp(at))*lineages[k].Theta(i);
		}
		else{
			mu = (2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1))*(exp(at)*cells[j].Xn(i)+(1-exp(at))*lineages[k].Theta(i)) + genes[i].Initial_expression()/genes[i].Initial_dispersion()) * var;
		}
		double tmp1 = mu - lineages[k].Theta(i);
		double tmp2 = cells[j].Xn(i) - lineages[k].Theta(i);
		return tmp1 * tmp2;
	}

	double X0_minus_theta_times_Xn_minus_theta(int i, int j, int k, double t){
		double at = genes[i].Alpha() * t;
		double var = 1.0/(2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1)) + 1.0/genes[i].Initial_dispersion());
		double mu;
		if(var < 0.000001){
			mu = exp(at)*cells[j].Xn(i)+(1-exp(at))*lineages[k].Theta(i);
		}
		else{
			mu = (2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1))*(exp(at)*cells[j].Xn(i)+(1-exp(at))*lineages[k].Theta(i)) + genes[i].Initial_expression()/genes[i].Initial_dispersion()) * var;
		}
		double tmp1 = mu - lineages[k].Theta(i);
		double tmp2 = cells[j].Xn(i) - lineages[k].Theta(i);
		return tmp1 * tmp2;
	}

	//for null model
	double X0_minus_theta_times_Xn_minus_theta(int i, int j){
		double at = genes[i].Alpha() * cells[j].Time();
		double var = 1.0/(2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1)) + 1.0/genes[i].Initial_dispersion());
		double mu = (2*genes[i].Alpha()/(genes[i].Sigma_squared()*(exp(2*at)-1))*(exp(at)*cells[j].Xn(i)+(1-exp(at))*genes[i].Theta_null()) + genes[i].Initial_expression()/genes[i].Initial_dispersion()) * var;
		double tmp1 = mu - genes[i].Theta_null();
		double tmp2 = cells[j].Xn(i) - genes[i].Theta_null();
		return tmp1 * tmp2;
	}

	double Log_likelihood(){
		double ll = 0.0;
		double tmp1, tmp2;
		for(int i=0; i<_cell_num; i++){
			for(int j=0; j<_K; j++){
				tmp1 = log(lineages[j].Pi());
				for(int k=0; k<_gene_num; k++){
					if(cells[i]._missing[k] != 1){
						tmp1 += genes[k].LogOU(cells[i].Xn(k), cells[i].Time(), j);
					}
				}

				if(j == 0)
					tmp2 = tmp1;
				else
					tmp2 = logsumexp(tmp2, tmp1);
			}
		
			ll += tmp2;
		}
		return ll;
	}

	double Log_likelihood_null_model(){
		double ll = 0.0;
		double tmp1;
		for(int i=0; i<_cell_num; i++){
			tmp1 = 0.0;
			for(int k=0; k<_gene_num; k++){
				if(cells[i]._missing[k] != 1){
					tmp1 += genes[k].LogOU_null_model(cells[i].Xn(k), cells[i].Time());
				}
			}
		
			ll += tmp1;
		}
		return ll;
	}

	double Log_likelihood_of_cell(int id){
		double tmp1, tmp2;

		for(int j=0; j<_K; j++){
			tmp1 = log(lineages[j].Pi());
			for(int k=0; k<_gene_num; k++){
				if(cells[id]._missing[k] != 1){
					tmp1 += genes[k].LogOU(cells[id].Xn(k), cells[id].Time(), j);
				}
			}

			if(j == 0)
				tmp2 = tmp1;
			else
				tmp2 = logsumexp(tmp2, tmp1);
		}

		return tmp2;
	}

	double Log_likelihood_of_cell(int id, double tmp_time){
		double tmp1, tmp2;
		
		for(int j=0; j<_K; j++){
			tmp1 = log(lineages[j].Pi());
			for(int k=0; k<_gene_num; k++){
				if(cells[id]._missing[k] != 1){
					tmp1 += genes[k].LogOU(cells[id].Xn(k), tmp_time, j);
				}
			}

			if(j == 0)
				tmp2 = tmp1;
			else
				tmp2 = logsumexp(tmp2, tmp1);
		}

		return tmp2;
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

	void Set_initial_parameter2(){
		srand((unsigned)time(NULL));
		int sample_size = _cell_num/_K;
		double tmp_theta;

		for(int i=0; i<_gene_num; i++){
			for(int j=0; j<_K; j++){
				tmp_theta = 0.0;
				for(int k=0; k<sample_size; k++){
					tmp_theta += cells[k].Xn(i);
				}
				tmp_theta /= sample_size;
				lineages[j].Add_theta(i, tmp_theta);
			}
		}
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

	int Set_time(FILE *fp){
		int count, id;
		double time;
		for(int i=0; i<_cell_num; i++){
			count = fscanf(fp, "%d\t%lf\n", &id, &time);
			if(count == EOF){
				printf("error at reading initial parameter\n");
				return 1;
			}
			
			if(time < _min_time){
				time = _min_time;
			}
			cells[i].Add_time(time);
		}
		return 0;
	}

	void Set_null_parameter(){
		double ave;
		for(int i=0; i<_gene_num; i++){
			ave = 0.0;
			for(int j=0; j<_K; j++){
				ave += lineages[j].Theta(i);
			}
			ave /= _K;

			genes[i].Add_theta_null(ave);
		}
	}

	//Initialization of mixutre OU
	void Set_parameter(){
		for(int i=0; i<_gene_num; i++){
			for(int j=0; j<_K; j++){
				lineages[j].Add_theta(i, genes[i].Theta_null());
			}
		}

		//random responsibility
		srand((unsigned) time(NULL));
		for(int i=0; i<_cell_num; i++){
			cells[i].Random_responsibility();
		}
	}

	int Set_optimized_parameter(FILE *fin_gene, FILE *fin_cell){
		char buf[100000];
		char *tmp;

		//read pi
		fgets(buf, 1000000, fin_gene);
		tmp = strtok(buf, "\t");
		tmp = strtok(NULL, "\t");
		for(int k=0; k<_K-1; k++){
			tmp = strtok(NULL, "\t");
			lineages[k].Add_new_pi(atof(tmp));
		}
		tmp = strtok(NULL, "\n");
		lineages[_K-1].Add_new_pi(atof(tmp));

		//read gene parameters
		for(int g=0; g<_gene_num; g++){
			fgets(buf, 1000000, fin_gene);

			//read alpha
			tmp = strtok(buf, "\t");
			genes[g].Add_new_alpha(atof(tmp));
			//read sigma squared
			tmp = strtok(NULL, "\t");
			genes[g].Add_new_sigma_squared(atof(tmp));
			for(int k=0; k<_K-1; k++){
				tmp = strtok(NULL, "\t");
				lineages[k].Add_new_theta(g, atof(tmp));
			}
			tmp = strtok(NULL, "\n");
			lineages[_K-1].Add_new_theta(g, atof(tmp));

			genes[g].Update_parameter();
		}

		for(int k=0; k<_K; k++){
			lineages[k].Update_parameter();
		}

		//read cell time
		for(int c=0; c<_cell_num; c++){
			fgets(buf, 1000000, fin_cell);
			tmp = strtok(buf, "\t");
			cells[c].Add_new_time(atof(tmp));
			cells[c].Update_parameter();
		}

		return 0;
	}

	int Set_optimized_null_parameter(FILE *fp){
		char buf[100000];
		char *tmp;

		for(int i=0; i<_gene_num; i++){
			fgets(buf, 1000000, fp);
			tmp = strtok(buf, "\t");
			tmp = strtok(NULL, "\n");
			genes[i].Add_theta_null(atof(tmp));
		}

		return 0;
		//todo
		for(int i=0; i<_gene_num; i++){
			fgets(buf, 1000000, fp);
			tmp = strtok(buf, "\t");
			tmp = strtok(NULL, "\n");
			genes[i].Add_new_alpha(atof(tmp));

			fgets(buf, 1000000, fp);
			tmp = strtok(buf, "\t");
			tmp = strtok(NULL, "\n");
			genes[i].Add_new_sigma_squared(atof(tmp));

			genes[i].Update_parameter();
		}

		for(int i=0; i<_cell_num; i++){
			fgets(buf, 1000000, fp);
			tmp = strtok(buf, "\t");
			tmp = strtok(NULL, "\n");
			cells[i].Add_new_time(atof(tmp));
			cells[i].Update_parameter();
		}

		return 0;
	}

	int Convergence(){
		for(int k=0; k<_K; k++){
			if(lineages[k].Convergence(_max_time, _thresh) == 0){
				//printf("at lineage %d\n", k);
				return 0;
			}
		}

		for(int g=0; g<_gene_num; g++){
			if(genes[g].Convergence(_max_time, _thresh) == 0){
				//printf("at gene %d\n", g);
				return 0;
			}
		}

		for(int c=0; c<_cell_num; c++){
			if(cells[c].Convergence(_thresh) == 0){
				//printf("at cell %d\n", c);
				return 0;
			}
		}

		return 1;
	}

	//for debug
	void Print_responsibility(){
		for(int i=0; i<_cell_num; i++){
			printf("t=%lf\t", cells[i].Time());
			for(int j=0; j<_K; j++){
				printf("%LF\t", cells[i]._gamma[j]);
			}
			printf("\n");
		}
	}

	void Print_cell_parameter(FILE *fp){
		for(int i=0; i<_cell_num; i++){
			fprintf(fp, "%lf\t", cells[i].Time());
			for(int j=0; j<_K; j++){
				fprintf(fp, "%.3LF\t", cells[i]._gamma[j]);
			}
			fprintf(fp, "\n");
		}
	}

	void Print_gene_parameter(){
		for(int i=0; i<_K; i++){
			printf("pi lineage[%d] is %lf\n", i, lineages[i].Pi());
			for(int j=0; j<_gene_num; j++){
				printf("theta gene[%d] lineage[%d] is %lf\n", i, j, lineages[i].Theta(j));
			}
		}

		for(int i=0; i<_gene_num; i++){
			printf("alpha[%d]\t%lf\n", i, genes[i].Alpha());
			printf("sigma_squared[%d]\t%lf\n", i, genes[i].Sigma_squared());
			printf("theta\t%lf\n", lineages[0].Theta(i));
		}

		return;
	}

	void Print_gene_parameter(FILE *fp){
		fprintf(fp, " \t \t");
		for(int k=0; k<_K; k++){
			fprintf(fp, "%lf", lineages[k].Pi());
			if(k != _K-1)	fprintf(fp, " \t");
			else	fprintf(fp, "\n");
		}

		for(int g=0; g<_gene_num; g++){
			fprintf(fp, "%lf\t%lf\t", genes[g].Alpha(), genes[g].Sigma_squared());
			for(int k=0; k<_K; k++){
				fprintf(fp, "%lf", lineages[k].Theta(g));
				if(k != _K-1)	fprintf(fp, " \t");
				else	fprintf(fp, "\n");
			}
		}
	}

	void Print_gene_nullparameter(FILE *fp){
		for(int g=0; g<_gene_num; g++){
			fprintf(fp, "%lf\t%lf\t%lf\n", genes[g].Alpha(), genes[g].Sigma_squared(), genes[g].Theta_null());
		}
	}

	void Print_ll(FILE *fp){
		fprintf(fp, "%lf\n", Log_likelihood());
	}

	void Print_correlation(FILE *fp, FILE *fcor){
		int id;
		double at, mean, variance, normalized_value;
		vector<vector<double> > expr(_gene_num, vector<double>(_cell_num, 0));

		for(int i=0; i<_cell_num; i++){
			for(int j=0; j<_gene_num; j++){
				//todo mixture
				id = 0;

				at = genes[j].Alpha() * cells[i].Time();
				mean = exp(-at)*genes[j].Initial_expression() + (1 - exp(-at))*lineages[id].Theta(j);
				variance = genes[j].Sigma_squared()*(1-exp(-2*at))/(2*genes[j].Alpha()) + exp(-2*at)*genes[j].Initial_dispersion();
				normalized_value = (cells[i].Xn(j) - mean)/(sqrtf(variance));
				
				expr[j][i] = normalized_value;
			}
		}

		//print normalized expression
		for(int g=0; g<_gene_num; g++){
			for(int c=0; c<_cell_num; c++){
				fprintf(fp, "%lf", expr[g][c]);
				if(c != _cell_num-1){
					fprintf(fp, "\t");
				}
				else{
					fprintf(fp, "\n");
				}
			}
		}

		//correlation
		int effective_size;
		double mean1, mean2, cov, var1, var2, cor;
		for(int i=0; i<_gene_num; i++){
			for(int j=0; j<_gene_num; j++){
				effective_size = 0;
				mean1 = 0;
				mean2 = 0;
				for(int c=0; c<_cell_num; c++){
					mean1 += expr[i][c];
					mean2 += expr[j][c];
					effective_size++;
				}
				mean1 /= effective_size;
				mean2 /= effective_size;

				cov = 0;
				var1 = 0;
				var2 = 0;
				for(int c=0; c<_cell_num; c++){
					cov += (expr[i][c] - mean1) * (expr[j][c] - mean2);
					var1 += (expr[i][c] - mean1) * (expr[i][c] - mean1);
					var2 += (expr[j][c] - mean2) * (expr[j][c] - mean2);
				}

				if(effective_size <= 1){
					cor = 0;
				}
				else if(i == j){
					cor = 0;
				}
				else{
					cor = cov/(sqrtf(var1) * sqrtf(var2));
				}
				fprintf(fcor, "%lf", cor);

				if(j != _gene_num-1){
					fprintf(fcor, "\t");
				}
				else{
					fprintf(fcor, "\n");
				}
			}
		}
	}

	/*
	//for debug
	void Check_parameter_optimization(){
		srand((unsigned) time(NULL));

		for(int i=0; i<5; i++){
			string tmp="out/check_time_";
			tmp += to_string(i);
			tmp += ".txt";
			FILE *fout = fopen(tmp.c_str(), "w");
			int cell_id = rand()%100;
			double opt_time = cells[cell_id].Time();
			fprintf(fout, "%lf\t%lf\n", opt_time, Log_likelihood());
			for(double t=0.1; t<2.0; t+= 0.1){
				cells[cell_id].Add_new_time(t);
				cells[cell_id].Update_parameter();
				fprintf(fout, "%lf\t%lf\n", t, Log_likelihood());
			}

			cells[cell_id].Add_new_time(opt_time);
			cells[cell_id].Update_parameter();
			fclose(fout);
		}

		for(int i=0; i<5; i++){
			string tmp="out/check_alpha_theta_";
			tmp += to_string(i);
			tmp += ".txt";
			FILE *fout = fopen(tmp.c_str(), "w");
			int gene_id = rand()%500;
			double opt_alpha = genes[gene_id].Alpha();
			double opt_theta = lineages[0].Theta(gene_id);
			fprintf(fout, "%lf\t%lf\t%lf\n", opt_alpha, opt_theta, Log_likelihood());
			for(double a=0.1; a<10; a+= 0.5){
				for(double theta=-10; theta<11; theta+=1){
					genes[gene_id]._alpha = a;
					lineages[0]._theta[gene_id] = theta;
					fprintf(fout, "%lf\t%lf\t%lf\n", a, theta, Log_likelihood());
				}
			}

			genes[gene_id]._alpha = opt_alpha;
			lineages[0]._theta[gene_id] = opt_theta;
		}

		for(int i=0; i<5; i++){
			string tmp="out/check_sigma_theta_";
			tmp += to_string(i);
			tmp += ".txt";
			FILE *fout = fopen(tmp.c_str(), "w");
			int gene_id = rand()%500;
			double opt_sigma = genes[gene_id].Sigma_squared();
			double opt_theta = lineages[0].Theta(gene_id);
			fprintf(fout, "%lf\t%lf\t%lf\n", opt_sigma, opt_theta, Log_likelihood());
			for(double sigma=0.1; sigma<100; sigma+= 5.0){
				for(double theta=-10; theta<11; theta+=1){
					genes[gene_id]._sigma_squared = sigma;
					lineages[0]._theta[gene_id] = theta;
					fprintf(fout, "%lf\t%lf\t%lf\n", sigma, theta, Log_likelihood());
				}
			}

			genes[gene_id]._sigma_squared = opt_sigma;
			lineages[0]._theta[gene_id] = opt_theta;
		}

		for(int i=0; i<5; i++){
			string tmp="out/check_alpha_sigma_";
			tmp += to_string(i);
			tmp += ".txt";
			FILE *fout = fopen(tmp.c_str(), "w");
			int gene_id = rand()%500;
			double opt_alpha = genes[gene_id].Alpha();
			double opt_sigma = genes[gene_id].Sigma_squared();
			fprintf(fout, "%lf\t%lf\t%lf\n", opt_alpha, opt_sigma, Log_likelihood());
			for(double a=0.1; a<10; a+= 0.5){
				for(double sigma=0.1; sigma<100; sigma+= 5.0){
					genes[gene_id]._alpha = a;
					genes[gene_id]._sigma_squared = sigma;
					fprintf(fout, "%lf\t%lf\t%lf\n", a, sigma, Log_likelihood());
				}
			}

			genes[gene_id]._alpha = opt_alpha;
			genes[gene_id]._sigma_squared = opt_sigma;
		}
	}
	*/
};

