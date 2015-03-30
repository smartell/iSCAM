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
	class compositionLikelihoods
	{
	private:
		int r1,r2;		/// index for first and last row.
		int c1,c2;		/// index for rist and last column.
		ivector m_jmin; /// index for ragged start columns.
		ivector m_jmax;	/// index for ragged end columns.

		T m_O;  		/// observed composition object

	public:
		virtual const T nloglike(const T & _O) const = 0;
		virtual const T residual(const T & _O) const = 0;
		virtual ~compositionLikelihoods(){};
		
		void set_O(T & O) { this -> m_O = O; }
		T    get_O() const{ return m_O;      }

		void tail_compression_indexes();
	};

	void acl::test();
	


	// |--------------------------------------------------------------------------------|
	// | MULTIVARIATE LOGISTIC NEGATIVE LOGLIKELIHOOD derived class                     |
	// |--------------------------------------------------------------------------------|
	template<class T>
	inline
	const T dmvlogistic(const T& O, const T& P)
	{

	}

	/**
	 * @brief Negative loglikelihood for compostion data using multivaritae logistic 
	 * @author Steve Martell
	 * @param T usually a dmatrix
	 */
	 template<class T>
	 class multivariteLogistic: public compositionLikelihoods<T>
	 {
	 private:
	 	T m_P; 		/// predicted composition object.

	 public:
	 	multivariteLogistic(T _P): m_P(_P) {}

	 	const T nloglike(const T &O) const
	 	{
	 		return acl::dmvlogistic(O, this->get_P());
	 	}
	 	
	 	T get_P() const {return m_P; }
	 	void set_P(T P) {this->m_P=P;}

	 };


}




#endif