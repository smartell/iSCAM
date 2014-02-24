// milka.h

#ifndef MILKA_H
#define MILKA_H

#include <admodel.h>

/**
 * @defgroup Milka Operating model for iSCAM
 * @Milka Classes and functions for the operating model
 * 
 * @author  Steven Martell
 * @deprecated  Feb 23, 2014
 * 
 * @details The namespace is mse, short for management strategy evaluation.
 * <br> The following class objects are: <br><br>
 * <br> Class
 */

namespace mse {
	class omData{
	private:
		// |------------------|
		// | Model dimensions |
		// |------------------|
		int m_nStock;
		int m_nArea;
		int m_nSex;
		int m_nSyr;
		int m_nNyr;
		int m_nSage;
		int m_nNage;
		int m_nGear;
		int m_nFleet;

	public:
		omData();
	   ~omData();

	   void set_nStock(const int &n) { m_nStock = n; }
	   void set_nArea (const int &n) { m_nArea  = n; }
	   void set_nSex  (const int &n) { m_nSex   = n; }
	   void set_nSyr  (const int &n) { m_nSyr   = n; }
	   void set_nNyr  (const int &n) { m_nNyr   = n; }
	   void set_nSage (const int &n) { m_nSage  = n; }
	   void set_nNage (const int &n) { m_nNage  = n; }
	   void set_nGear (const int &n) { m_nGear  = n; }
	   void set_nFleet(const int &n) { m_nFleet = n; }

	   friend class OperatingModel;
	};

	class omVariables{

	};

	class OperatingModel{
	private:
		omData       m_data;
		omVariables  m_vars;
	protected:
	public:

	};
} // mse namespace


#endif