#include <vector>
#include <algorithm>

class Cutoff{
	public:
		double NPDF(double, double, double);
		double EPDF(double, double);
		void normalize();
		void readin(const std::string&);
		void zero();
		void init(int);
		void resize();
		int findValue(const std::string&);

	private:
		double mu_best;
		double mu_step;
		double mu_old = 0;
		double sigma_best;
		double sigma_step;
		double sigma_old = 0;
		double lambda_best;
		double lambda_step;
		double lambda_old = 0;
		double pnorm_best;
		double pnorm_step = 0.5;
		double pnorm_old = 0;
		double KLD_best = 10000;
		double KLD;
		double M = 0.5;
		double XE = 0;
		double XN = 0;
		double X2N = 0;
		double CE = 0;
		double CN = 0;
		std::vector<double> X_i;
		std::vector<double> P_i;
		int N = 0;
};
