#include <Rcpp.h>
#include <vector>

/* internal C++: computes the contribution of each ancestral edge for a species to the weights matrix */
void process_epochs ( std::vector<double> const &epochs, 
					  double const &alphasq, 
					  std::vector<double> &y ) 
{
	int i, np = epochs.size();
	double t;
	
	for (i = 0; i < np; i++) {
		t = epochs[0]-epochs[i];
		y.push_back( exp(-alphasq*t) );
	}
	for (i = 0; i < np-1; i++) {
		y[i] -= y[i+1];
	}
}

/* internal C++: builds a binary array for each species, signifying in which selective regime the ancestor of the species is bounded */
void process_beta ( std::vector<int> &beta,
				    std::vector<int> const &ancestors,
				    std::vector<double> const &regimes,
				    std::vector<double> const &optima )
{
	int np, r, k, nreg=optima.size();

	np=ancestors.size();
	beta.assign(np*nreg, 0);
	
	for(r=0; r<np; r++) {
		for(k=0; k<nreg; k++){
			if(regimes.at( ancestors.at(r)-1 )==optima.at(k))
			{
				beta.at(r+k*np)=1;
			}
		}
	}
}


/* C++ | R INTERFACE; main function */
RcppExport SEXP ouweights ( SEXP PHYLO ) 
{
	/* 
	 * PHYLO: 
	 *	- epochs: list of branching times for each segment from tip to root [T=t1>t2>...>troot=0] for each species
	 *	- ancestors: list of nodes from tip to root for each species 
	 *	- regimes: array of optima associated with each edge; order from node 1 to node 2N-2 
	 *	- optima: array of unique selective regimes across the tree
	 *	- alpha: constraint factor for the Ornstein-Uhlenbeck process; herein, alpha is expected given as alpha^2
	 *  - tips: number of species in tree
	 */
			
	try {
		Rcpp::List phylo(PHYLO);
		
		Rcpp::List epochs = phylo["epochs"];
		Rcpp::List ancestors = phylo["ancestors"];
		std::vector<double> regimes =  phylo["regimes"];
		std::vector<double> optima = phylo["optima"];
		double alphasq = Rcpp::as<double> ( phylo["alpha"] );
		int ntip = Rcpp::as<int> ( phylo["tips"] );
		int nreg = optima.size();
		
		int t, r, q;
		int segments;
		std::vector<double> y;
		std::vector<double> W;
		W.assign(ntip*nreg,0);
		std::vector<int> tbeta;

		for(t=0; t<ntip; t++) {
			std::vector<int> tancestors=ancestors[t];
			process_beta(tbeta, tancestors, regimes, optima); 
			std::vector<double> tepochs=epochs[t];
			segments=tepochs.size();
			y.reserve(segments);
			process_epochs(tepochs, alphasq, y);
			for(r=0; r<nreg; r++) {
				for(q=0; q<segments; q++) {
					W[t+r*ntip] += y[q]*tbeta[q+r*segments];
				}
			}
			std::vector<double>().swap(y);
			std::vector<double>().swap(tepochs);
			std::vector<int>().swap(tancestors);
			std::vector<int>().swap(tbeta);
		}
		
		/* PREPARE OUTPUT FOR R */
		return Rcpp::List::create(Rcpp::Named("weights", W));

	} catch( std::exception &ex ) {		
		forward_exception_to_r( ex );
    } catch(...) { 
		::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; 
}

