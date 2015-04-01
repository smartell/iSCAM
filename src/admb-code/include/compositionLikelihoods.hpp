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
 */
namespace acl
{

	// template<class T>
	// inline
	// const T tailCompression(const T &_M, const double pmin=0.0)
	// {
	// 	int r1,r2,c1,c2;
	// 	r1 = _M.rowmin();
	// 	r2 = _M.rowmax();
	// 	c1 = _M.colmin();
	// 	c2 = _M.colmax();

	// 	ivector jmin(r1,r2);
	// 	ivector jmax(r1,r2);
	// 	for(int i = r1; i <= r2; i++ )
	// 	{				
	// 		auto o = _M(i);

	// 		jmin(i) = c1+1;			// index for min column
	// 		jmax(i) = c2;			// index for max column
	// 		auto cumsum = o/sum(o);
	// 		for(int j = c1+1; j <= c2; j++ )
	// 		{
	// 			 cumsum(j) += cumsum(j-1);
	// 			 cumsum(j) <= pmin ? jmin(i)++ : NULL;
	// 			 j != c2 ? 1.0 - cumsum(j) < pmin ? jmax(i)-- : NULL : NULL;
	// 		}
	// 	}

	// 	// Now compress the matrix.
	// 	T M_;
	// 	T M = _M;
	// 	M_.allocate(r1,r2,jmin,jmax);
	// 	M_.initialize();

	// 	// fill ragged array M_
	// 	for(int i = r1; i <= r2; i++ )
	// 	{
	// 		M(i) /= sum(M(i));
	// 		M_(i)(jmin(i),jmax(i)) = M(i)(jmin(i),jmax(i));

	// 		// add cumulative sum to tails.
	// 		M_(i)(jmin(i)) = sum(M(i)(c1,jmin(i)));
	// 		M_(i)(jmax(i)) = sum(M(i)(jmax(i),c2));
	// 	}
	// 	return M_;

	// }

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

		DATA m_O;  		/// observed composition object
		DVAR m_P;       /// predicted composition object

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
		// inline
		// const DVAR tailCompression(const DVAR &_M) const;

		// methods
		inline
		void getRaggedVectors(double pmin = 0.0);

		// getters & setters
		void set_O(DATA & O) { this -> m_O = O; }
		DATA get_O()    const{ return m_O;      }

	 	DVAR get_P() const  {return m_P; }
	 	void set_P(DVAR &P) {this->m_P=P;}

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
				 cumsum(j) <= pmin ? m_jmin(i)++ : NULL;
				 j != c2 ? 1.0 - cumsum(j) <= pmin ? m_jmax(i)-- : NULL : NULL;
				 
				 cumsum(j+1) += cumsum(j);
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
		T var = 1.0/size_count(NU)*norm2(NU);
		return size_count(NU) * log(var);
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
	 	DATA m_rO; 		/// ragged observed  composition object.
	 	DVAR m_rP; 		/// ragged predicted composition object.
	 	DVAR m_nu;		/// logistic residuals.

	 public:
	 	// constructor
	 	// todo: add eps value to constructor.
	 	multivariteLogistic(const DATA &_O, const DVAR &_P, const double eps=1.e-8)
	 	:negLogLikelihood<DATA,DVAR>(_O,_P) 
	 	{
	 		// tail compression
	 		DATA tmp = this->compress(this->get_O());
	 		set_rO(tmp+eps);
	 		DVAR vmp = this->compress(this->get_P());
	 		set_rP(vmp+eps);


	 		// residuals
	 		DVAR tnu = acl::dmvlogisticResidual(this->get_rO(),this->get_rP());
	 		set_nu(tnu);
	 	}


	 	const T nloglike() const
	 	{
	 		return acl::dmvlogistic<T,DVAR>(this->get_nu());
	 	}

	 	const DVAR residual() const
	 	{
	 		return acl::dmvlogisticResidual(this->get_O(), this->get_P());
	 	}
	 	
	 	DATA get_rO() const {return m_rO; }
	 	DVAR get_rP() const {return m_rP; }
	 	DVAR get_nu() const {return m_nu; }

	 	void set_rO(DATA &X){this->m_rO=X;}
	 	void set_rP(DVAR &X){this->m_rP=X;}
	 	void set_nu(DVAR &R){this->m_nu=R;}

	 };


}




#endif