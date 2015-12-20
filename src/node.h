#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <cassert>
#include <cfloat>

using namespace std;

class GENE;
class CELL;
class LINEAGE;

vector<GENE> genes;
vector<CELL> cells;
vector<LINEAGE> lineages;

double Gaussian(double x, double mean, double sigma_squared){
	return 1.0/sqrt(2 * M_PI * sigma_squared) * exp(-(x-mean)*(x-mean)/(2*sigma_squared));
}

double Log_Gaussian(double x, double mean, double sigma_squared){
	return -0.5 * log(2 * M_PI * sigma_squared) - (x-mean)*(x-mean)/(2*sigma_squared);
}

double logsumexp(double x, double y){
	if(isinf(x)){
		printf("inf at logsumexp %lf\n", x);
		return y;
	}
	else if(isinf(y)){
		printf("inf at logsumexp %lf\n", y);
		return x;
	}
	else if(x > y) return (x + log1p(exp(-x+y)));
	else return (y + log1p(exp(-y+x)));
}

class LINEAGE{
public:
	int _gene_num;
	double _pi;
	double _new_pi;
	vector<double> _theta;
	vector<double> _new_theta;

	void Init(int g, double pi){
		_gene_num = g;
		_theta.resize(g, 0);
		_new_theta.resize(g, 0);
		_pi = pi;
	}

	double Pi(){
		return _pi;
	}

	void Add_theta(int gene_id, double theta){
		_theta[gene_id] = theta;
	}

	void Add_new_theta(int gene_id, double theta){
		_new_theta[gene_id] = theta;
	}

	void Add_new_pi(double pi){
		_new_pi = pi;
	}

	void Normalize_new_pi(double sum){
		_new_pi /= sum;
	}

	double Theta(int gene_id){
		return _theta[gene_id];
	}

	void Update_parameter(){
		_pi = _new_pi;
		for(int i=0; i<_gene_num; i++){
			_theta[i] = _new_theta[i];
		}
	}

	int Convergence(double t, double thresh);
};

class GENE{
public:
	int _gene_id;
	int _K;
	double _alpha;
	double _sigma_squared;
	double _new_alpha;
	double _new_sigma_squared;
	double _tmp_alpha;
	double _tmp_sigma_squared;

	double _initial_expression;
	double _initial_sigma_squared;
	double _new_initial_expression;
	double _new_initial_sigma_squared;

	double _theta_null;
	double _new_theta_null;

	double _sigma_squared_MG;

	void Init(int k, int gene_id, double alpha, double sigma){
		_gene_id = gene_id;
		_K = k;
		_alpha = alpha;
		_sigma_squared = 1.0;
	}

	void Init(int gene_id){
		_gene_id = gene_id;
	}

	void Add_initial_expression(double expression){
		_initial_expression = expression;
	}

	void Add_new_initial_expression(double expression){
		_new_initial_expression = expression;
	}

	void Add_initial_dispersion(double dispersion){
		_initial_sigma_squared = dispersion;
	}

	void Add_new_initial_dispersion(double dispersion){
		_new_initial_sigma_squared = dispersion;
	}

	void Add_new_alpha(double alpha){
		_new_alpha = alpha;
	}

	void Add_new_sigma_squared(double sigma){
		_new_sigma_squared = sigma;
	}

	void Add_theta_null(double theta){
		_theta_null = theta;
	}

	void Add_new_theta_null(double theta){
		_new_theta_null = theta;
	}

	void Add_sigma_squared_MG(double sigma){
		_sigma_squared_MG = sigma;
	}

	double OU(double x, double t, int k){
		double mean = exp(-_alpha * t)*_initial_expression + (1 - exp(-_alpha*t))*lineages[k].Theta(_gene_id);
		double variance = _sigma_squared*(1-exp(-2*_alpha*t))/(2*_alpha);
		
		double ret = Gaussian(x, mean, variance);
		return ret;
	}

	double LogOU(double x, double t, int k){
		double mean = exp(-_alpha * t)*_initial_expression + (1 - exp(-_alpha*t))*lineages[k].Theta(_gene_id);
		double variance = _sigma_squared*(1-exp(-2*_alpha*t))/(2*_alpha) + exp(-2*_alpha*t)*_initial_sigma_squared;
		
		double ret = Log_Gaussian(x, mean, variance);
		return ret;
	}

	double LogOU_with_new_parameter(double x, double t, int k){
		double mean = exp(-_new_alpha * t)*_initial_expression + (1 - exp(-_new_alpha*t))*lineages[k]._new_theta[_gene_id];
		double variance = _new_sigma_squared*(1-exp(-2*_new_alpha*t))/(2*_new_alpha) + exp(-2*_new_alpha*t)*_initial_sigma_squared;
		
		double ret = Log_Gaussian(x, mean, variance);
		if(isinf(ret)){
			//for debug
			printf("-inf at OU %lf %lf %lf\n", x, mean, variance);
			return -DBL_MIN;
		}
		else{
			return ret;
		}
	}

	//todo
	double LogOU_with_new_theta(double x, double t, int k){
		double mean = exp(-_alpha * t)*_initial_expression + (1 - exp(-_alpha*t))*lineages[k]._new_theta[_gene_id];
		double variance = _sigma_squared*(1-exp(-2*_alpha*t))/(2*_alpha) + exp(-2*_alpha*t)*_initial_sigma_squared;
		
		double ret = Log_Gaussian(x, mean, variance);
		if(isinf(ret)){
			//for debug
			printf("-inf at OU %lf %lf %lf\n", x, mean, variance);
			return -DBL_MIN;
		}
		else{
			return ret;
		}
	}

	double LogOU_with_new_sigma(double x, double t, int k, double sigma){
		double mean = exp(-_alpha * t)*_initial_expression + (1 - exp(-_alpha*t))*lineages[k].Theta(_gene_id);
		double variance = sigma*(1-exp(-2*_alpha*t))/(2*_alpha) + exp(-2*_alpha*t)*_initial_sigma_squared;
		
		double ret = Log_Gaussian(x, mean, variance);
		if(isinf(ret)){
			//for debug
			printf("-inf at OU %lf %lf %lf\n", x, mean, variance);
			return -DBL_MIN;
		}
		else{
			return ret;
		}
	}

	double LogOU_with_new_alpha(double x, double t, int k, double alpha){
		double mean = exp(-alpha * t)*_initial_expression + (1 - exp(-alpha*t))*lineages[k].Theta(_gene_id);
		double variance = _sigma_squared*(1-exp(-2*alpha*t))/(2*alpha) + exp(-2*alpha*t)*_initial_sigma_squared;
		
		double ret = Log_Gaussian(x, mean, variance);
		if(isinf(ret)){
			//for debug
			printf("-inf at OU %lf %lf %lf\n", x, mean, variance);
			return -DBL_MIN;
		}
		else{
			return ret;
		}
	}

	double LogOU_null_model(double x, double t){
		double mean = exp(-_alpha * t)*_initial_expression + (1 - exp(-_alpha*t))*_theta_null;
		double variance = _sigma_squared*(1-exp(-2*_alpha*t))/(2*_alpha) + exp(-2*_alpha*t)*_initial_sigma_squared;
		
		double ret = Log_Gaussian(x, mean, variance);
		if(isinf(ret)){
			//for debug
			printf("-inf at OU %lf %lf %lf\n", x, mean, variance);
			return -DBL_MIN;
		}
		else{
			return ret;
		}
	}

	double LogOU_null_model2(double x, double t){
		double mean = exp(-_alpha * t)*_initial_expression + (1 - exp(-_alpha*t))*_initial_expression;
		double variance = _sigma_squared*(1-exp(-2*_alpha*t))/(2*_alpha) + exp(-2*_alpha*t)*_initial_sigma_squared;
		
		double ret = Log_Gaussian(x, mean, variance);
		if(isinf(ret)){
			//for debug
			printf("-inf at OU %lf %lf %lf\n", x, mean, variance);
			return -DBL_MIN;
		}
		else{
			return ret;
		}
	}

	double Mean_of_OU(double x, double t, int k){
		double mean = exp(-_alpha * t)*_initial_expression + (1 - exp(-_alpha*t))*lineages[k].Theta(_gene_id);
		return mean;
	}

	double Expectation(double t, int k){
		return exp(-_alpha * t)*_initial_expression + (1 - exp(-_alpha*t))*lineages[k].Theta(_gene_id);
	}

	double Alpha(){
		return _alpha;
	}

	double Sigma_squared(){
		return _sigma_squared;
	}

	double Initial_expression(){
		return _initial_expression;
	}

	double Initial_dispersion(){
		return _initial_sigma_squared;
	}

	double X0(){
		return _initial_expression;
	}

	double Theta_null(){
		return _theta_null;
	}

	void Update_parameter(){
		_alpha = _new_alpha;
		_sigma_squared = _new_sigma_squared;
	}

	void Update_null_parameter(){
		_theta_null = _new_theta_null;
		_alpha = _new_alpha;
		_sigma_squared = _new_sigma_squared;
	}

	void Update_null_parameter2(){
		_theta_null = _new_theta_null;
	}

	void Updata_null_parameter2_2(){
		_alpha = _new_alpha;
		_sigma_squared = _new_sigma_squared;
	}

	void Store_parameter(){
		_tmp_alpha = _alpha;
		_tmp_sigma_squared = _sigma_squared;
	}

	void Restore_parameter(){
		_alpha = _tmp_alpha;
		_sigma_squared = _tmp_sigma_squared;
	}

	int Convergence(double t, double thresh){
		double tmp1, tmp2;

		tmp1 = (_sigma_squared*(1-exp(-_alpha*t)))/(2*_alpha);
		tmp2 = (_new_sigma_squared*(1-exp(-_new_alpha*t)))/(2*_new_alpha);
		if(fabs(tmp1 - tmp2) > thresh){
			return 0;
		}

		return 1;
	}
};

class CELL{
public:
	int _gene_num;
	int _K;
	int _flag_differentiation;

	double _old_time;
	double _time;
	double _new_time;
	vector<double> _expression;
	vector<int> _missing;
	vector<long double> _gamma;
	vector<vector<long double> > _gamma_gene;
	vector<double> _time_gene;

	void Init(int g, int k){
		_gene_num = g;
		_K = k;
		_expression.resize(g, 0);
		_missing.resize(g, 0);
		_gamma.resize(k, 0);
		_gamma_gene.resize(_gene_num);
		for(int i=0; i<_gene_num; i++){
			_gamma_gene[i].resize(_K, 0.0);
		}
		_time_gene.resize(g, 0);
	}

	void Init(int g){
		_gene_num = g;
		_expression.resize(g, 0);
		_missing.resize(g, 0);
	}

	void Init_for_gene_time(){
		for(int i=0; i<_gene_num; i++){
			_time_gene[i] = _time;
		}
	}

	void Add_expression(int id, double expression, int missing_flag){
		if(missing_flag == 0){
			_expression[id] = expression;
		}
		else{
			_expression[id] = 0;
			_missing[id] = 1;
		}
	}

	void Add_time(double t){
		_time = t;
	}

	void Add_new_time(double t){
		_new_time = t;
	}

	void Add_new_gene_time(int id, double t){
		_time_gene[id] = t;
	}

	void Add_flag_dif(int flag){
		_flag_differentiation = flag;
	}

	int Flag_dif(){
		return _flag_differentiation;
	}

	double Get_expression(int id){
		return _expression[id];
	}

	double Time(){
		return _time;
	}

	double Gene_time(int id){
		return _time_gene[id];
	}

	double Gamma(int id){
		return _gamma[id];
	}
	double Gamma_gene(int gene, int state){
		return _gamma_gene[gene][state];
	}

	double Xn(int id){
		return _expression[id];
	}

	void Calc_responsibility(){
		//calculate in log
		for(int i=0; i<_K; i++){
			_gamma[i] = log(lineages[i].Pi());
			for(int j=0; j<_gene_num; j++){
				//todo
				if(_missing[j] == 0){
					_gamma[i] += genes[j].LogOU(_expression[j], _time, i);
				}
			}
		}
		
		//logsum version
		long double logsum = _gamma[0];
		for(int i=1; i<_K; i++){
			logsum = logsumexp(logsum, _gamma[i]);
		}
 
		//todo
		if(isinf(logsum)){
			printf("aa %LF %LF %LF\n", logsum, _gamma[0], _gamma[1]);

			int max_id = -1;
			double max = -DBL_MIN;
			for(int i=0; i<_K; i++){
				if(max < _gamma[i]){
					max_id = i;
					max = _gamma[i];
				}
			}
			for(int i=0; i<_K; i++){
				if(i == max_id){
					_gamma[i] = 1.0;
				}
				else{
					_gamma[i] = 0.0;
				}
			}
		}
		else{
			for(int i=0; i<_K; i++){
				_gamma[i] = expl(_gamma[i] - logsum);
			}
		}

		//todo
		//add pseudocount to avoid overfitting
		double sum=0.0;
		for(int i=0; i<_K; i++){
			_gamma[i] += 0.01;
			sum += _gamma[i];
		}
		for(int i=0; i<_K; i++){
			_gamma[i] /= sum;
		}
	}

	/*
	void Calc_missing_value(){
		for(int i=0; i<_gene_num; i++){
			//todo
			if(_missing[i] == 0){
				continue;
			}

			//todo
			double mean = 0;
			for(int j=0; j<_K; j++){
				mean += _gamma[j] * genes[i].Mean_of_OU(_expression[i], _time, j);
			}
			_expression[i] = mean;
		}
	}
	*/

	void Random_responsibility(){
		srand((unsigned) time(NULL));
		//calculate in log
		double sum = 0;
		for(int i=0; i<_K; i++){
			_gamma[i] = ((double)rand())/RAND_MAX;
			sum += _gamma[i];
		}
		
		for(int i=0; i<_K; i++){
			_gamma[i] /= sum;
		}
	}

	void Calc_responsibility_of_gene(){
		for(int i=0; i<_gene_num; i++){
			//calculate in log
			for(int j=0; j<_K; j++){
				_gamma_gene[i][j] = 0;
				//todo
				if(_missing[i] == 0){
					_gamma_gene[i][j] += genes[i].LogOU(_expression[i], _time, j);
				}
			}
			
			//logsum version
			double logsum = _gamma_gene[i][0];
			for(int j=1; j<_K; j++){
				logsum = logsumexp(logsum, _gamma_gene[i][j]);
			}
			
			for(int j=0; j<_K; j++){
				_gamma_gene[i][j] = exp(_gamma_gene[i][j] - logsum);
			}		
		}
	}

	double Expectation(int g, double t, int k){
		if(t > _time){
			double dt = t - _time;
			double at = genes[g].Alpha() * dt;

			//todo
			return exp(-at)*Xn(g) + (1-exp(-at))*lineages[k].Theta(g);
		}
		else{
			double mean1, mean2, var1, var2;

			double dt = _time - t;
			double at = genes[g].Alpha() * dt;
			mean1 = exp(at)*Xn(g) + (1-exp(at))*lineages[k].Theta(g);
			var1 = exp(2*at)*genes[g].Sigma_squared()*(1-exp(-2*at))/(2*genes[g].Alpha());

			at = genes[g].Alpha() * t;
			mean2 = exp(-at)*genes[g].Initial_expression() + (1-exp(-at))*lineages[k].Theta(g);
			var2 = genes[g].Sigma_squared()*(1-exp(-2*at))/(2*genes[g].Alpha()) + exp(-2*at)*genes[g].Initial_dispersion();

			return (var2*mean1 + var1*mean2)/(var1 + var2);
		}
	}

	double Var(int g, double t, int k){
		if(t > _time){
			double dt = t - _time;
			double at = genes[g].Alpha() * dt;

			return genes[g].Sigma_squared()*(1-exp(-2*at))/(2*genes[g].Alpha());
		}
		else{
			double var1, var2;

			double dt = _time - t;
			double at = genes[g].Alpha() * dt;
			var1 = exp(2*at)*genes[g].Sigma_squared()*(1-exp(-2*at))/(2*genes[g].Alpha());

			at = genes[g].Alpha() * t;
			var2 = genes[g].Sigma_squared()*(1-exp(-2*at))/(2*genes[g].Alpha()) + exp(-2*at)*genes[g].Initial_dispersion();

			return (var1*var2)/(var1 + var2);
		}
	}

	void Update_parameter(){
		_old_time = _time;
		_time = _new_time;
	}

	int Convergence(double thresh){
		if(fabs(_old_time - _new_time) > thresh){
			return 0;
		}

		return 1;
	}
};

int LINEAGE::Convergence(double t, double thresh){
	if(fabs(_pi - _new_pi) > thresh){
		return 0;
	}

	double tmp1, tmp2;
	for(int g=0; g<_gene_num; g++){
		tmp1 = exp(- genes[g]._new_alpha * t) * _new_theta[g];
		tmp2 = exp(- genes[g]._alpha * t) * _theta[g];
		if(fabs(tmp1 - tmp2) > thresh){
			return 0;
		}
	}
	return 1;
}

