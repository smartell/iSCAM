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
	class msy //: public referencePoints<T,T1,T2>
	{  
	private:
		// Indexes for dimensions
		int m_sage;
		int m_nage;
		int m_nGear;
		int m_nGrp;
		

		T m_ro;
		T m_bo;
		T m_h;
		T m_rho;	    /// Fraction of mortality that occurs before spawning.
		T m_phie;		/// Spawning biomass per recruit in unfished conditions.
		T m_phif;		/// Spawning biomass per recruit in fished conditions.

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
			//m_Va.allocate(*V);
			//m_Va = *V;

			if(m_Ma.indexmin() != m_Fa.indexmin() || m_Ma.indexmax() != m_Fa.indexmax())
			{
			cerr<<"Indexes do not mach in calcPhie"<<endl;
			exit(1);
			}

			m_nGrp = m_Ma.rowmax() - m_Ma.rowmin() + 1;
			m_nGear = m_Va(1).rowmax();
			cout<<"NGEAR "<<m_nGear<<endl;
			
			m_sage = m_Ma.colmin();
			m_nage = m_Ma.colmax();

			calcPhie();

			cout<<"In constructor\n"<<m_phie<<endl;
		}

		
		virtual const T1 getFmsy(const T1 &fe);
		
	};

	template<class T, class T1, class T2, class T3>
	const T1 msy<T,T1,T2,T3>::getFmsy(const T1 & fe)
	{
		calcSurvivorship(fe);
		return(0);	
	}

	template<class T, class T1, class T2, class T3>
	void msy<T,T1,T2,T3>::calcSurvivorship(const T1 &fe)
	{
		cout<<"Working on this routine"<<endl;
		int j,h,k;
		T phif = 0.0;
		

		T1 pza(m_sage,m_nage);
		T1 psa(m_sage,m_nage);
		T1 poa(m_sage,m_nage);
		T2  za(1,m_nGrp,m_sage,m_nage);
		T2  sa(1,m_nGrp,m_sage,m_nage);
		T2  oa(1,m_nGrp,m_sage,m_nage);
		T2  lz(1,m_nGrp,m_sage,m_nage);
		T2  lw(1,m_nGrp,m_sage,m_nage);

		T2   qa(1,m_nGear,m_sage,m_nage);
		T2  dlz(1,m_nGear,m_sage,m_nage);
		T2 d2lz(1,m_nGear,m_sage,m_nage);
		T2  dlw(1,m_nGear,m_sage,m_nage);
		T2 d2lw(1,m_nGear,m_sage,m_nage);
		dlz.initialize();
		dlw.initialize();
		d2lz.initialize();
		d2lw.initialize();

		T3   qa_m(1,m_nGrp,1,m_nGear,m_sage,m_nage);
		T3  dlz_m(1,m_nGrp,1,m_nGear,m_sage,m_nage);
		T3 d2lz_m(1,m_nGrp,1,m_nGear,m_sage,m_nage);
		T3  dlw_m(1,m_nGrp,1,m_nGear,m_sage,m_nage);
		T3 d2lw_m(1,m_nGrp,1,m_nGear,m_sage,m_nage);

		cout<<m_Va<<endl;
		for( h = 1; h <= m_nGrp; h++ )
		{
			za(h) = m_Ma(h);
			for( k = 1; k <= m_nGear; k++ )
			{
				za(h) = za(h) + fe(k) * m_Va(h)(k);
				cout<<h<<" "<<k<<endl;
			}
			sa(h) = exp(-za(h));
			oa(h) = 1.0 - sa(h);
			pza   = m_rho*za(h);
			psa   = exp(-pza);
			poa   = 1.0 - elem_prod(sa(h),psa);
			for(k=1;k<=m_nGear;k++)
			{
				qa(k)      = elem_div(elem_prod(elem_prod(m_Va(h)(k),m_Wa(h)),oa(h)),za(h));
				qa_m(h)(k) = qa(k);

				dlw(k,m_sage)      = -psa(m_sage)*m_rho*m_Va(h)(k)(m_sage);
				dlw_m(h,k,m_sage)  = dlw(k,m_sage);
				d2lw(k,m_sage)     = psa(m_sage)*square(m_rho)*square(m_Va(h)(k)(m_sage));
				d2lw_m(h,k,m_sage) = d2lw(k,m_sage);

				dlz(k,m_sage)      = 0;
				dlz_m(h,k,m_sage)  = 0;
				d2lz(k,m_sage)     = 0;
				d2lz_m(h,k,m_sage) = 0;
			}

			// Survivorship
			lz(h)(m_sage) = 1.0/m_nGrp;
			lw(h)(m_sage) = 1.0/m_nGrp * psa(m_sage);
			for( j = m_sage+1; j <= m_nage; j++ )
			{
				lz(h,j) = lz(h,j-1) * sa(h,j-1);
				lw(h,j) = lw(h,j-1) * psa(j);
				if( j == m_nage )
				{
					lz(h,j) = lz(h,j)/oa(h,j);
					lw(h,j) = lz(h,j-1)*sa(h,j-1)*psa(j)/oa(h,j);
				}

				for( k = 1; k <= m_nGear; k++ )
				{
					// derivatives for survivorship
					dlz(k)(j)  = sa(h)(j-1)*( dlz(k)(j-1)-lz(h)(j-1)*m_Va(h)(k)(j-1));
					d2lz(k)(j) = sa(h)(j-1)*(d2lz(k)(j-1)+lz(h)(j-1)*square(m_Va(h)(k)(j-1)));
					
					// derivatives for spawning survivorship
					dlw(k)(j)  = -lz(h)(j)*m_rho*m_Va(h)(k)(j)*psa(j); 
					d2lw(k)(j) =  lz(h)(j)*square(m_rho)*square(m_Va(h)(k)(j))*psa(j);

					if( j == m_nage ) // + group derivatives
					{
						dlz(k)(j)  = dlz(k)(j)/oa(h)(j) 
						             - lz(h)(j-1)*sa(h)(j-1)*m_Va(h)(k)(j)*sa(h)(j)
						             /square(oa(h)(j));
						
						dlw(k)(j)  = -lz(h)(j-1)*sa(h)(j-1)*m_rho*m_Va(h)(k)(j)/oa(h)(j)
									- lz(h)(j-1)*psa(j)*m_Va(h)(k)(j)*sa(h)(j)
									/square(oa(h)(j));
						
						T V1  	   = m_Va(h)(k)(j-1);
						T V2  	   = m_Va(h)(k)(j);
						T oa2 	   = oa(h)(j)*oa(h)(j);
						
						d2lz(k)(j) = d2lz(k)(j)/oa(h)(j) 
									+ 2*lz(h)(j-1)*V1*sa(h)(j-1)*V2*sa(h)(j)/oa2
									+ 2*lz(h)(j-1)*sa(h)(j-1)*V2*V2*sa(h)(j)*sa(h)(j)
									/(oa(h)(j)*oa2)
									+ lz(h)(j-1)*sa(h)(j-1)*V2*V2*sa(h)(j)/oa2;
						
						d2lw(k)(j) = lz(h)(j-1)*square(m_rho)*square(V2)*psa(j)/oa(h)(j)
									+ 2*lz(h)(j-1)*m_rho*square(V2)*psa(j)*sa(h)(j)/oa2
									+ 2*lz(h)(j-1)*psa(j)*square(V2)*square(sa(h)(j))
									/(oa(h)(j)*oa2)
									+ lz(h)(j-1)*psa(j)*square(V2)*sa(h)(j)/oa2;
					}
					dlz_m(h,k,j)  =  dlz(k)(j);
					d2lz_m(h,k,j) = d2lz(k)(j);
					dlw_m(h,k,j)  =  dlw(k)(j);
					d2lw_m(h,k,j) = d2lw(k)(j);
				} // m_nGear
			} // m_nage

			// Spawning biomass per recruit in fished conditions.
			phif += lz(h) * m_Fa(h);

		} // m_nGrp
		m_phif = phif;

		cout<<m_phif<<endl;
		cout<<"End of CalcSurvivorship"<<endl;
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
		m_bo = m_ro * m_phie;
		cout<<"Bo = "<<m_bo<<endl;
		
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