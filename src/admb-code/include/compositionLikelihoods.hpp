// compositionLikelihoods.hpp
#include <admodel.h>

#ifndef COMPOSITIONLIKELIHOODS_H
#define COMPOSITIONLIKELIHOODS_H

/**
 * @file compositionLikelihoods.hpp
 * @defgroup Likelihoods
 * @author Steven Martell
 * @namespace acl
 * @date   Feb 10, 2014
 * @title Selectivity functions
 * @details  Uses abstract base class for computing negative loglikelihoods
 */
namespace acl
{

	template<class T>
	inline
	const T tailCompression(const T &_M, const double pmin=0.0)
	{
		int r1,r2,c1,c2;
		r1 = _M.rowmin();
		r2 = _M.rowmax();
		c1 = _M.colmin();
		c2 = _M.colmax();

		ivector jmin(r1,r2);
		ivector jmax(r1,r2);
		for(int i = r1; i <= r2; i++ )
		{				
			auto o = _M(i);

			jmin(i) = c1+1;			// index for min column
			jmax(i) = c2;			// index for max column
			auto cumsum = o/sum(o);
			for(int j = c1+1; j <= c2; j++ )
			{
				 cumsum(j) += cumsum(j-1);
				 cumsum(j) <= pmin ? jmin(i)++ : NULL;
				 j != c2 ? 1.0 - cumsum(j) < pmin ? jmax(i)-- : NULL : NULL;
			}
		}

		// Now compress the matrix.
		T M_;
		T M = _M;
		M_.allocate(r1,r2,jmin,jmax);
		M_.initialize();

		// fill ragged array M_
		for(int i = r1; i <= r2; i++ )
		{
			M(i) /= sum(M(i));
			M_(i)(jmin(i),jmax(i)) = M(i)(jmin(i),jmax(i));

			// add cumulative sum to tails.
			M_(i)(jmin(i)) = sum(M(i)(c1,jmin(i)));
			M_(i)(jmax(i)) = sum(M(i)(jmax(i),c2));
		}
		return M_;

	}

	/**
	 * @brief Base class for composition likelihoods.
	 * @details Virtual methods for negative loglikelihood and standardized residuals.
	 * 
	 * @tparam T [description]
	 */
	template<class T>
	class compositionLikelihoods
	{
	private:
		int r1,r2;		/// index for first and last row.
		int c1,c2;		/// index for rist and last column.
		ivector m_jmin; /// index for ragged start columns.
		ivector m_jmax;	/// index for ragged end columns.

		T m_O;  		/// observed composition object

	public:
		// constructors
		compositionLikelihoods(){}
		compositionLikelihoods(const T& _O)
		:m_O(_O)
		{
			r1 = m_O.rowmin();
			r2 = m_O.rowmax();
			c1 = m_O.colmin();
			c2 = m_O.colmax();
		}

		// virtual methods
		virtual ~compositionLikelihoods(){};
		virtual const T nloglike(const T & _O) const = 0;
		virtual const T residual(const T & _O) const = 0;
		
		// getters & setters
		void set_O(T & O) { this -> m_O = O; }
		T    get_O() const{ return m_O;      }
	
	};

	

	// |--------------------------------------------------------------------------------|
	// | MULTIVARIATE LOGISTIC NEGATIVE LOGLIKELIHOOD derived class                     |
	// |--------------------------------------------------------------------------------|
	template<class T>
	inline
	const T dmvlogistic(const T& O, const T& P)
	{
		T nll;
		return nll;
	}

	template<class T>
	inline
	const T dmvlogisticResidual(const T& O, const T& P)
	{
		T nu;
		return nu;
	}

	/**
	 * @brief Negative loglikelihood for compostion data using multivaritae logistic 
	 * @author Steve Martell
	 * @param T usually a dmatrix
	 */
	 template<class DATA, class DVAR>
	 class multivariteLogistic: public compositionLikelihoods<T>
	 {
	 private:
	 	T m_P; 		/// predicted composition object.

	 public:
	 	multivariteLogistic(const DATA &_O, const DVAR &_P)
	 	:compositionLikelihoods<T>(_O),m_P(_P) {}

	 	const T nloglike(const T &O) const
	 	{
	 		T rO = tailCompression(O);
	 		return acl::dmvlogistic(O, this->get_P());
	 	}

	 	const T residual(const T &O) const
	 	{
	 		return acl::dmvlogisticResidual(O, this->get_P());
	 	}
	 	
	 	T get_P() const {return m_P; }
	 	void set_P(T P) {this->m_P=P;}

	 };


}




#endif