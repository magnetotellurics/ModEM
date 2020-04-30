#include <string>
#include <vector>
#include <iostream>
//
#include "DataEntryBase.h"
//
#ifndef DATA_ENTRY_CSEM_H_
#define DATA_ENTRY_CSEM_H_
//
// Dipole Moment(Am) Azi Dip X Y Z
//
class DataEntryCSEM : public DataEntryBase
{
	public:
		// TYPE, DIPOLE, PERIOD, MOMENT, AZI, DIP, TX_X, TX_Y TX_Z, CODE, X, Y, Z, Component Real Imag, Error
		DataEntryCSEM( std::string, std::string, double, double, double, double, double, double, double,
				std::string, double, double, double, std::string, double, double, double  );
		~DataEntryCSEM( void );
		//
		void show( void );

	//
	// GETTERS & SETTERS
	//
	double getAzimuth( void ) const;
	void setAzimuth( double );
	//
	double getDip( void ) const;
	void setDip( double );
	//
	const std::string& getDipole( void ) const;
	void setDipole( const std::string& );
	//
	double getMoment( void ) const;
	void setMoment( double );
	//
	double getTxX( void ) const;
	void setTxX( double );
	//
	double getTxY( void ) const;
	void setTxY( double );
	//
	double getTxZ( void ) const;
	void setTxZ( double );
	//
	protected:
		std::string dipole;
		double moment, azimuth, dip, tx_x, tx_y, tx_z;
};

#endif
