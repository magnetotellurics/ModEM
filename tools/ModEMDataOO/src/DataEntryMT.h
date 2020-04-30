#include <string>
#include <vector>
#include <iostream>
//
#include "DataEntryBase.h"
//
#ifndef DATA_ENTRY_MT_H_
#define DATA_ENTRY_MT_H_
//
// GG_Lat GG_Lon
//
class DataEntryMT : public DataEntryBase
{
	public:
		// TYPE, PERIOD, CODE, LAT, LONG, X, Y, Z, COMPONENT, REAL, IMAGINARY, ERROR
		DataEntryMT( std::string, double, std::string, double, double, double, double, double, std::string, double, double, double );
		~DataEntryMT( void );
		//
		void show( void );
		//
		// GETTERS & SETTERS
		//
		double getLat( void ) const;
		void setLat( double );
		//
		double getLong( void ) const;
		void setLong( double );

	protected:
		double latitude, longitude;
};

#endif
