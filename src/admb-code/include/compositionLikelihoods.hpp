// negLogLikelihood.hpp
#include <admodel.h>

#ifndef negLogLikelihood_H
#define negLogLikelihood_H

/**
 * @file negLogLikelihood.hpp
 * @defgroup Likelihoods
 * @author Steven Martell
 * @namespace acl
 * @date   Feb 10, 2014
 * @title Selectivity functions
 * @details  Uses abstract base class for computing negative loglikelihoods
 * TODO:  Still need to fix templates so it works with df1b2variables.
 */
namespace acl
{

	/**
	 * @brief Base class for composition likelihoods.
	 * @details Virtual methods for negative loglikelihood and standardized residuals.
	 * 
	 * @tparam DATA Data type argument.
	 * @tparam DVAR Variable type argument.
	 */
	template<class DATA, class DVAR>
	class negLogLikelihood
	{
	private:
		int r1,r2;		/// index for first and last row.
		int c1,c2;		/// index for rist and last column.
		ivector m_jmin; /// index for ragged start columns.
		ivector m_jmax;	/// index for ragged end columns.
		imatrix m_jagg; /// index for aggregation among bins.

		DATA m_O;  		/// observed composition object
		DATA m_rO; 		/// observed ragged composition object

		DVAR m_P;       /// predicted composition object
		DVAR m_rP;      /// predicted ragged composition object
		

	public:
		// constructors
		negLogLikelihood(){}
		negLogLikelihood(const DATA& _O, const DVAR& _P)
		:m_O(_O), m_P(_P)
		{
			r1 = m_O.rowmin();
			r2 = m_O.rowmax();
			c1 = m_O.colmin();
			c2 = m_O.colmax();

			getRaggedVectors();
		}

		// pure virtual methods
		virtual ~negLogLikelihood(){};
		virtual const dvariable nloglike() const = 0;
		virtual const DVAR      residual() const = 0;
		
		// virtual methods
		inline
		const DATA tailCompression(const DATA &_M) const;
		
		// methods
		inline
		void getRaggedVectors(double pmin = 0.0);

		inline
		void aggregate(const double pmin = 0.0);

		// getters & setters
		void set_O(DATA & O) { this -> m_O = O; }
		DATA get_O()    const{ return m_O;      }

	 	DVAR get_P() const  {return m_P; }
	 	void set_P(DVAR &P) {this->m_P=P;}

	 	DATA get_rO() const {return m_rO; }
	 	void set_rO(DATA &X){this->m_rO=X;}

	 	DVAR get_rP() const {return m_rP; }
	 	void set_rP(DVAR &X){this->m_rP=X;}

	 	imatrix get_jagg() const {return m_jagg;}

		template <typename T>
		inline
		const T compress(const T& _M) const
		{
			T R;
			T M = _M;
			R.allocate(r1,r2,m_jmin,m_jmax);
			R.initialize();
			// fill ragged array R
			for(int i = r1; i <= r2; i++ )
			{
				M(i) /= sum(M(i));
				R(i)(m_jmin(i),m_jmax(i)) = M(i)(m_jmin(i),m_jmax(i));

				// add cumulative sum to tails.
				R(i)(m_jmin(i)) = sum(M(i)(c1,m_jmin(i)));
				R(i)(m_jmax(i)) = sum(M(i)(m_jmax(i),c2));
			}
			
			return R;
		}
	
	};

	template<class DATA, class DVAR>
	inline
	void acl::negLogLikelihood<DATA,DVAR>::aggregate(const double pmin)
	{
		m_jmin.allocate(r1,r2);
		m_jmax.allocate(r1,r2);
		ivector n(r1,r2);
		n.initialize();

		// get number of observations > pmin each year.		
		for(int i = r1; i <= r2; i++ )
		{
			dvector oo = m_O(i)/sum(m_O(i));
			for(int j = c1; j <= c2; j++ )
			{
				if( oo(j) > pmin ) n(i)++;
			}		
		}

		m_rO.allocate(r1,r2,1,n);
		m_rP.allocate(r1,r2,1,n);
		m_jagg.allocate(r1,r2,1,n);
		m_rO.initialize();
		m_rP.initialize();
		m_jagg.initialize();

		for(int i = r1; i <= r2; i++ )
		{
			int k = 1;
			dvector     oo = m_O(i)/sum(m_O(i));
			dvar_vector pp = m_P(i)/sum(m_P(i));
			for(int j = c1; j <= c2; j++ )
			{
				m_rO(i)(k) += oo(j);
				m_rP(i)(k) += pp(j);
				if(k<=n(i)) m_jagg(i,k) = j;
				if( oo(j)>pmin && k<n(i)) k++;
			}		
		}

	}

	/**
	 * @brief Get column indicies for each row to compress matrix into ragged array.
	 * @details Determine the start and end positions of each row in which to compute
	 * the likelihoods where values < pmin are pooled into the tails of the distribution.
	 * 
	 * @param pmin value to assume for the minium proportion.
	 * @return NULL Modifies private member variables: m_jmin and m_jmax.
	 */
	template<class DATA, class DVAR>
	inline
	void acl::negLogLikelihood<DATA,DVAR>::getRaggedVectors(double pmin)
	{
		m_jmin.allocate(r1,r2);
		m_jmax.allocate(r1,r2);
		
		for(int i = r1; i <= r2; i++ )
		{				
			dvector o = m_O(i);

			m_jmin(i) = c1;			// index for min column
			m_jmax(i) = c2;			// index for max column
			dvector cumsum = o/sum(o);
			for(int j = c1; j < c2; j++ )
			{
				// min bin and plus bin column index.
				 cumsum(j)       <= pmin ? m_jmin(i)++ : NULL;
				 1.0 - cumsum(j) <= pmin ? m_jmax(i)-- : NULL;
				 cumsum(j+1)     += cumsum(j);

				 // o(j) > pmin ? n(i)++ : NULL;
			}
		}
		
	}
	















	// |--------------------------------------------------------------------------------|
	// | MULTIVARIATE LOGISTIC NEGATIVE LOGLIKELIHOOD derived class                     |
	// |--------------------------------------------------------------------------------|
	
	/**
	 * @brief multivariate logistic negative log likelihood
	 * @details The negative log likeihood for the multivariate 
	 * logistic distribution based on the MLE of the variance.
	 * 
	 * @param NU Matrix of residuals
	 * @return the negative loglikelihood
	 */
	template<class T, class DVAR>
	inline
	const T dmvlogistic(const DVAR& NU)
	{
		int n = size_count(NU) - (NU.rowmax()-NU.rowmin()+1);
		T var = 1.0/n * norm2(NU);
		return n * log(var);
	}

	template<class DATA, class DVAR>
	inline
	const DVAR dmvlogisticResidual(const DATA& O, const DVAR& P)
	{
		DVAR nu;
		nu.allocate(P);
		int i;
		for( i = O.rowmin(); i<= O.rowmax(); i++)
		{
			nu(i) = log(O(i)) - log(P(i));
			nu(i) -= mean(nu(i));
		}
		return nu;
	}

	/**
	 * @brief Negative loglikelihood for compostion data using multivaritae logistic 
	 * @author Steve Martell
	 * @param T usually a dmatrix
	 */
	 template<class T, class DATA, class DVAR>
	 class multivariteLogistic: public negLogLikelihood<DATA,DVAR>
	 {
	 private:

	 	DVAR m_nu;		/// logistic residuals based on ragged object.
	 	DVAR m_NU;		/// residuals for rectangular matrix

	 public:
	 	// constructor
	 	multivariteLogistic(const DATA &_O, const DVAR &_P, const double eps=0)
	 	:negLogLikelihood<DATA,DVAR>(_O,_P) 
	 	{
	 		// tail compression
	 		DATA tmp = this->compress(this->get_O()) + eps;
	 		this->set_rO(tmp);
	 		DVAR vmp = this->compress(this->get_P()) + eps;
	 		this->set_rP(vmp);


	 		// residuals
	 		this->aggregate(eps);
	 		DVAR tnu = acl::dmvlogisticResidual(this->get_rO(),this->get_rP());
	 		set_nu(tnu);

	 		// residuals to return.
	 		m_NU.allocate(_P);
	 		m_NU.initialize();
	 		imatrix idx = this->get_jagg();
	 		for(int i = _O.rowmin(); i <= _O.rowmax(); i++ )
	 		{
	 			for(int j = 1; j <= size_count(m_nu(i)); j++ )
	 			{
	 				int jj = idx(i,j); /* code */
		 			m_NU(i)(jj) = m_nu(i)(j);
	 			}
	 		}
	 	}


	 	const T nloglike() const
	 	{
	 		return acl::dmvlogistic<T,DVAR>(this->get_nu());
	 	}

	 	const DVAR residual() const
	 	{
	 		// This should return the residuals for the original container.
	 		// Not the residuals for the ragged object.

	 		return this->get_NU();
	 	}
	 	
	 	DVAR get_nu() const {return m_nu; }
	 	DVAR get_NU() const {return m_NU; }
	 	void set_nu(DVAR &R){this->m_nu=R;}


	 };


	 // |-------------------------------------------------------------------------------|
	 // | MULTINOMIAL DISTRIBUTION                                                      |
	 // |-------------------------------------------------------------------------------|

	 template <class T, class DATA, class DVAR>
	 class multinomial: public negLogLikelihood<DATA,DVAR>
	 {
	 public:
	 	// constructor
	 	multinomial(const DATA &_O, const DVAR &_P, const double eps=0)
	 	:negLogLikelihood<DATA,DVAR>(_O,_P)
	 	{

	 	}
	 	
	 	const T nloglike() const
	 	{
	 		return 0;
	 	}
	 	
	 	const DVAR residual() const
	 	{
	 		return 0;
	 	}
	 };

}




#endif