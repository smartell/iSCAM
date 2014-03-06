#ifndef _MSY_H
#define _MSY_H

#define MAXITER 300
#define TOL     1.e-04

#include <admodel.h>
#include <fvar.hpp>


namespace rfp {

	/**
	 * @brief Base class for reference point calculations
	 * @details The base class for MSY and SPR-based reference point 
	 * calculations.
	 * 
	 * @param fe vector of fishing mortality rates
	 * @return getFmsy is a pure virtual function.
	 */
	template<typename T, typename T1, typename T2>
	class referencePoints 
	{
	private:
		T1 m_fe;

	public:
		

		// Pure virtual functions... set = 0;
		virtual const T1 getFmsy(const T1 &fe) const = 0;
		// virtual const Type getBmsy(const Type &be) const = 0;

		virtual ~referencePoints(){}

		void setFe(T1 & fe) { this -> m_fe = fe; }
		T1 getFe()    const { return m_fe;     }
		
	};



	/**
	 * @brief MSY-based reference points
	 * @details Class object for computing MSY-based reference points.
	 * 
	 * @author Steven Martell
	 * 
	 * @tparam T variable
	 * @tparam T1 vector
	 * @tparam T2 matrix
	 * @tparam T3 3d_array
	 */
	template<class T, class T1, class T2, class T3>
	class msy : public referencePoints<T,T1,T2>
	{  
	private:
		// Indexes for dimensions
		int m_sage;
		int m_nage;
		int m_nGear;
		int m_nGrp;
		

		T m_ro;
		T m_h;
		T m_rho;	    /// Fraction of mortality that occurs before spawning.
		T m_phie;		/// Spawning biomass per recruit.

		T1 m_fe;		/// Fishing mortality rate
		T1 m_ye;		/// Equilibrium yield.


		T2 m_Ma;		/// Natural mortality rate matrix.
		T2 m_Wa;		/// Weight-at-age matrix.
		T2 m_Fa;		/// Fecundity-at-age matrix.

		T3 m_Va;		/// Selectivity-at-age.

		void calcPhie();
		void calcSurvivorship(const T1 &fe);
		
		//void calcEquilibrium(const T1 &fe);

	public:
		// default constructor
		msy(const T  ro ,
		    const T  h  ,
		    const T  rho,
		    const T2 ma ,
		    const T2 wa ,
		    const T2 fa ,
		    const T3 V )
		:m_ro(ro),m_h(h),m_rho(rho),m_Ma(ma),m_Wa(wa),m_Fa(fa),m_Va(V) 
		{
			if(m_Ma.indexmin() != m_Fa.indexmin() || m_Ma.indexmax() != m_Fa.indexmax())
			{
			cerr<<"Indexes do not mach in calcPhie"<<endl;
			exit(1);
			}

			m_nGrp = m_Ma.rowmax() - m_Ma.rowmin() + 1;
			m_nGear = m_Va.indexmax() - m_Va.indexmin() + 1;
			cout<<"NGEAR "<<m_nGear<<endl;
			
			m_sage = m_Ma.colmin();
			m_nage = m_Ma.colmax();

			calcPhie();

			cout<<"In constructor\n"<<m_phie<<endl;
		}

		virtual const T1 getFmsy(const T1 &fe) const
		{
			calcSurvivorship(fe);

			return 0;
		}
		
	};


	template<class T, class T1, class T2, class T3>
	void msy<T,T1,T2,T3>::calcSurvivorship(const T1 &fe)
	{
		int h;
		T1  za(m_sage,m_nage);
		for( h = 1; h <= m_nGrp; h++ )
		{
			za = m_Ma(h);
		}

	}

	/**
	 * @brief Calculate spawning biomass per recruit.
	 * @details Calculate spawning biomass per recruit based on survivorship and maturity
	 * at age.
	 * 
	 * @author Steven Martell
	 * 
	 * @param rho Fraction of natural mortality that occurs prior to spawning.
	 * @param Ma Natural mortality rate at age and sex
	 * @param Fa Maturity-at-age and sex
	 * @tparam T1 Vector
	 * @tparam T2 Matrix
	 */
	template<class T, class T1, class T2, class T3>
	void msy<T,T1,T2,T3>::calcPhie()
	{
		int i,j;
		
		m_phie = 0;
		T2 lx(1,m_nGrp,m_sage,m_nage);
		lx.initialize();
		for( i = 1; i <= m_nGrp; i++ )
		{
			 for( j = m_sage; j <= m_nage; j++ )
			 {
			 	lx(i,j) = exp(-m_Ma(i,j)*(j-m_sage) - m_rho*m_Ma(i,j));
			 	if(j==m_nage) lx(i,j) /= 1.0 -exp(-m_Ma(i,j));
			 }
			 m_phie += 1./(m_nGrp) * lx(i) * m_Fa(i);
		}
		cout<<"From Bar "<<m_phie<<endl;
		
	}




	// template<class Type>
	// class spr : public referencePoints<Type>
	// {
	// public:
	// 	spr();
	// 	~spr();
	// };
} //rfp


#endif