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
		T m_x;
	public:
		virtual const T nloglike(const T & x) const = 0;
		virtual const T residual(const T & x) const = 0;
		virtual ~compositionLikelihoods(){};
		
		void set_x(T & x) { this -> m_x = x; }
		T    get_x() const{ return m_x;      }

	};


}

#endif