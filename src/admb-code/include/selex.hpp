//selex.hpp

#ifndef SELEX_HPP
#define SELEX_HPP

#include <admodel.h>
#include <iostream>
#include <assert.h>

/**
 * @defgroup Selectivities
 * @Selectivities All of the alternative selectivity functions in the SLX namespace are
 * derived from the slx::Selex base class.  
 * 
 * @author Steven Martell
 * @date   Feb 10, 2014
 * 
 * <br> Available Selectivity options are: <br><br>
 * <br>Selectivity              FUNCTIONS                Class name
 * <br>Logistic                 plogis                   LogisticCurve
 * <br>Nonparametric            nonparametric            SelectivityCoefficients
 * <br>
 */
namespace slx {
	
	#define Register(x) slxInterface<REAL_T>::template Init<x>();

	/**
	 * A Generic Selectivity Interface class.
	 */
	template<typename REAL_T>
	class slxInterface
	{
	private:

		typedef REAL_T(*EVAL_FUNCTION_PTR)(void* const);

		EVAL_FUNCTION_PTR eval_ptr;

		template<typename T> static
		REAL_T Evaluate_T(void* const pObj)
		{
			return static_cast<T*> (pObj)->Evaluate();
		}

	protected:

		// set a pointer to the evaluate function.
		template<typename T> void Init()
		{
			eval_ptr = (EVAL_FUNCTION_PTR) & Evaluate_T<T>;
		}

	public:

		slxInterface()
		: eval_ptr(0){

		}

		virtual ~ slxInterface(){

		}

		inline REAL_T Evaluate(){
			assert(eval_ptr);    // ensure Init() was called.
			return (*eval_ptr)(this);
		}

	};







	/**
	 * Logistic Selectivity Class that inherits from the slxInterface
	 */
	template<typename REAL_T>
	class slx_Logistic: public slxInterface<REAL_T>
	{
	private:
		friend class slxInterface<REAL_T>;

		// Evaluate should return a dvar_vector or a df1b2_vector
		inline REAL_T Evaluate() {
			//cout<<"Evaluating slx_Logistic"<<endl;
			dvariable mu = exp(m_log_mu);
			dvariable sd = exp(m_log_sd);
			
			return log( 1.0/(1.0+exp(-(m_x-mu)/sd)) );
		}

	protected:
		dvariable m_log_mu;
		dvariable m_log_sd;
		dvector   m_x;

	public:
		template<typename T1,typename T2>
		slx_Logistic(T1 &x, T2 &log_mu, T2 &log_sd)
		:m_log_mu(log_mu),m_log_sd(log_sd),m_x(x)
		{
			Register(slx_Logistic<REAL_T>);
		}
	};







	/**
	 * Selectivity coefficients.
	 */

	template<typename REAL_T>
	class slx_Coefficients: public slxInterface<REAL_T>
	{
	private:
		friend class slxInterface<REAL_T>;

		inline REAL_T Evaluate() {
			// cout<<"Evaluating slx_Coefficients"<<endl;
			int s1 = m_log_sel_coeffs.indexmin();
			int s2 = m_log_sel_coeffs.indexmax();
			int x1 = m_x.indexmin();
			int x2 = m_x.indexmax();
			REAL_T log_sel(x1,x2);

			for(int i = s1; i <= s2; i++ )
			{
				log_sel(i) = m_log_sel_coeffs(i);
			}
			log_sel(x1,s1) = log_sel(s1);
			log_sel(s2,x2) = log_sel(s2);
			return log_sel;
		}

	protected:
		dvector m_x;
		REAL_T m_log_sel_coeffs;

	public:

		template<typename T1,typename T2>
		slx_Coefficients(T1 &x, T2 &log_sel_coeffs)
		:m_x(x),m_log_sel_coeffs(log_sel_coeffs)
		{
			Register(slx_Coefficients<REAL_T>);
		}
	};






	/**
	 * Cubic spline
	 */
	template<typename REAL_T>
	class slx_CubicSpline: public slxInterface<REAL_T>
	{
	private:
		friend class slxInterface<REAL_T>;

		inline REAL_T Evaluate() {
			//cout<<"Evaluating CubicSpline"<<endl;
			int nodes = size_count(m_log_spline_knots);
			dvector ia(1,nodes);
			dvector fa = (m_x-min(m_x))/(max(m_x)-min(m_x));
			ia.fill_seqadd(0,1./(nodes-1));
			vcubic_spline_function ffa(ia,m_log_spline_knots);
			return ffa(fa);
		}

	protected:
		dvector m_x;
		REAL_T m_log_spline_knots;
	public:
		template<typename T1, typename T2>
		slx_CubicSpline(T1 &x, T2 &log_spline_coeffs)
		:m_x(x),m_log_spline_knots(log_spline_coeffs)
		{
			Register(slx_CubicSpline<REAL_T>);
		}
	};



	/**
	 * BiCubic spline
	 */
	template<typename REAL_T>
	class slx_BiCubicSpline: public slxInterface<REAL_T>
	{
	private:
		friend class slxInterface<REAL_T>;

		inline REAL_T Evaluate() {
			//cout<<"Evaluating BiCubicSpline"<<endl;
			//cout<<m_log_sel.rowmax()<<endl;
			//cout<<m_log_spline_knots.colmin()<<endl;
			bicubic_spline(m_x,m_y,m_log_spline_knots,m_log_sel);
			
			return(m_log_sel);
		}

	protected:
		dvector m_x;
		dvector m_y;
		dvar_matrix m_log_sel;
		REAL_T m_log_spline_knots;
	public:
		template<typename T1, typename T2>
		slx_BiCubicSpline(const T1 &x, const T1 &y, T2 &log_spline_coeffs, T2 &A)
		:m_x(x),m_y(y),m_log_sel(A),m_log_spline_knots(log_spline_coeffs)
		{
			Register(slx_BiCubicSpline<REAL_T>);
		}
	};












	template<class T,class T2>
	const T plogis(const T &x, const T2 &mean, const T2 & sd)
	{
		return T2(1.0)/(T2(1.0)+exp(-(x-mean)/sd));
	}
	
	/**
	 * @ingroup Selectivities
	 * @brief An abstract class for Selectivity functions.	
	 * @details Classes that derive from this class overload the pure virtual functions:<br><br>
	 * const T Selectivity(const T &x) const <br>
	 * 
	 * @tparam x Independent variable (ie. age or size) for calculating selectivity.
	 */

	template<class T>
	class Selex
	{
	private:
		T m_x;

	public:
		virtual  const T Selectivity(const T &x) const = 0;

		virtual  const T logSelectivity(const T &x) const = 0;

		virtual  const T logSelexMeanOne(const T &x) const = 0;
		
		virtual ~Selex(){}

		void Set_x(T & x) { this-> m_x = x; }
		T    Get_x() const{ return m_x;     }
	};



	/**
	 * @brief Logistic curve
	 * @details Uses the logistic curve for a two parametere parametric function
	 * 
	 * @param  [description]
	 * @tparam T data vector or dvar vector
	 * @tparam T2 double or dvariable for mean and standard deviation of the logisit curve
	 */
	template<class T,class T2>
	class LogisticCurve: public Selex<T>
	{
	private:
		T2 m_mean;
		T2 m_std;

	public:
		LogisticCurve(T2 mean = T2(0), T2 std = T2(1))
		: m_mean(mean), m_std(std) {}

		T2 GetMean() const { return m_mean; }
		T2 GetStd()  const { return m_std;  }

		void SetMean(T2 mean) { this->m_mean = mean;}
		void SetStd(T2 std)   { this->m_std  = std; }

		const T Selectivity(const T &x) const
		{
			return slx::plogis<T>(x, this->GetMean(), this->GetStd());
		}

		const T logSelectivity(const T &x) const
		{
			return log(slx::plogis<T>(x, this->GetMean(), this->GetStd()));
		}

		const T logSelexMeanOne(const T &x) const
		{
			T y = log(slx::plogis<T>(x, this->GetMean(), this->GetStd()));
			y  -= log(mean(mfexp(y)));
			return y;
		}

	};





	/**
	 * @brief Nonparametric selectivity coefficients
	 * @details Assumes that the last age/size class has the same selectivity coefficient
	 * as the terminal sel_coeffs.
	 * 
	 * @param x Independent variable
	 * @param sel_coeffs Vector of estimated selectivity coefficients.
	 * @return Selectivity coefficients.
	 */
	template<class T>
	const T nonparametric(const T &x, const T &sel_coeffs)
	{
		int x1 = x.indexmin();
		int x2 = x.indexmax();
		int y2 = sel_coeffs.indexmax();
		T y(x1,x2);
		for(int i = x1; i < y2; i++ )
		{
			y(i) = sel_coeffs(i);
		}
		y(y2,x2) = sel_coeffs(y2);
		return y;
	}
	/**
	 * @brief Selectivity coefficients
	 * @details Age or size-specific selectivity coefficients for n-1 age/size classes
	 * 
	 * @tparam T vector of coefficients
	 */
	template<class T>
	class SelectivityCoefficients: public Selex<T>
	{
	private:
		T m_sel_coeffs;

	public:
		SelectivityCoefficients(T params = T(1))
		:m_sel_coeffs(params) {}

		T GetSelCoeffs() const { return m_sel_coeffs;    }
		void SetSelCoeffs(T x) { this->m_sel_coeffs = x; }

		const T Selectivity(const T &x) const
		{
			// Call the age specific function
			return slx::nonparametric(x, this->GetSelCoeffs());
		}

		const T logSelectivity(const T &x) const
		{
			// Call the age specific function
			return log(slx::nonparametric(x, this->GetSelCoeffs()));
		}

		const T logSelexMeanOne(const T &x) const
		{
			T y = log(slx::nonparametric(x, this->GetSelCoeffs()));
			y  -= log(mean(mfexp(y)));
			return y;
		}
	};

}//slx


#endif /* SELEX_HPP */   	