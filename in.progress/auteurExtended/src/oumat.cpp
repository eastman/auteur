#include <Rcpp.h>
#include <vector>

/* internal C++ */
void process_OU_covariance ( std::vector<double> const &VCV, 
							 double const &alphasq, 
							 double const &sigmasq, 
							 int const &ntip,
							 std::vector<double> &V) 
{	
	unsigned int i, j, nt=ntip;
	double sij, ti, tj, elti, eltj, W;
	
	for (i = 0; i < nt; i++) {
		for (j = 0; j <= i; j++) {
			ti = VCV[i+i*nt];
			sij = VCV[i+j*nt];
			tj = VCV[j+j*nt];
			elti = exp(-alphasq*(ti-sij));
			eltj = exp(-alphasq*(tj-sij));
			W = elti*sigmasq*eltj/(alphasq+alphasq);
			V[i+nt*j] += W;
			if (j != i) 
			{
				V[j+nt*i] += W;
			}
		}
	}
}


/* C++ | R INTERFACE; main function */
RcppExport SEXP oumat (
						   SEXP TIMES, 
						   SEXP ALPHA, 
						   SEXP SIGMA,
						   SEXP TIPS) 
{
	/* 
	 * TIMES: array form of VCV matrix (under Brownian motion): split times for pairs of species 
	 * ALPHA: herein, alpha is expected given as alpha^2
	 * SIGMA: herein, sigma is expected given as sigma^2
	 * TIPS: number of species
	 */
	
	try {
		Rcpp::List phylo(TIMES);
		std::vector<double> VCV =  phylo["vmat"];
		double alphasq = Rcpp::as<double> (ALPHA);
		double sigmasq = Rcpp::as<double> (SIGMA);
		unsigned int ntip = Rcpp::as<int> (TIPS);
		
		std::vector<double> V;
		V.assign(ntip*ntip,0);
		
		process_OU_covariance(VCV, alphasq, sigmasq, ntip, V);
		
		/* PREPARE OUTPUT FOR R */
		return Rcpp::List::create(Rcpp::Named("oumat", V));
		
	} catch( std::exception &ex ) {		
		forward_exception_to_r( ex );
    } catch(...) { 
		::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; 
}

