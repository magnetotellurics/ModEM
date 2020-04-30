#include "DataEntryCSEM.h"

#include <stddef.h>
#include <iterator>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
//
// TYPE, DIPOLE, PERIOD, MOMENT, AZI, DIP, TX_X, TX_Y TX_Z, CODE, X, Y, Z, Component Real Imag, Error
DataEntryCSEM :: DataEntryCSEM( std::string type, std::string dipole, double period, double moment, double azimuth,
		double dip, double tx_x, double tx_y, double tx_z, std::string code, double x, double y, double z,
		std::string component, double real, double imaginary, double error )
:DataEntryBase( type, period, code, x, y, z, component, real, imaginary, error )
{
	this->dipole = dipole;
	this->moment = moment;
	this->azimuth = azimuth;
	this->dip = dip;
	this->tx_x = tx_x;
	this->tx_y = tx_y;
	this->tx_z = tx_z;
	this->moment = moment;
}
DataEntryCSEM :: ~DataEntryCSEM( void ){}
//
void DataEntryCSEM :: show( void )
{
	std:: cout << "MT DATA ENTRY:" <<
	"\n\tperiod:\t" << period <<
	"\n\tcode:\t" << code <<
	"\n\tlatitude:\t" << latitude <<
	"\n\tlongitude:\t" << longitude <<
	"\n\t(x,y,z):\t(" << x << ", " << y << ", " << z << ")" <<
	"\n\tcomponent:\t" << component <<
	"\n\treal:\t" << real <<
	"\n\timaginary:\t" << imaginary <<
	"\n\terror:\t" << error  << std::endl;
}
//
// GETTERS & SETTERS
//
double DataEntryCSEM :: getAzimuth( void ) const
{
	return azimuth;
}
//
void DataEntryCSEM :: setAzimuth( double azimuth )
{
	this->azimuth = azimuth;
}
//
double DataEntryCSEM :: getDip( void ) const
{
	return dip;
}
//
void DataEntryCSEM :: setDip( double dip )
{
	this->dip = dip;
}
//
const std::string& DataEntryCSEM :: getDipole( void ) const
{
	return dipole;
}
//
void DataEntryCSEM :: setDipole( const std::string &dipole )
{
	this->dipole = dipole;
}
//
double DataEntryCSEM :: getMoment( void ) const
{
	return moment;
}
//
void DataEntryCSEM :: setMoment( double moment )
{
	this->moment = moment;
}
//
double DataEntryCSEM :: getTxX( void ) const
{
	return tx_x;
}
//
void DataEntryCSEM :: setTxX( double txX )
{
	tx_x = txX;
}
//
double DataEntryCSEM :: getTxY( void ) const
{
	return tx_y;
}
//
void DataEntryCSEM :: setTxY( double txY )
{
	tx_y = txY;
}
//
double DataEntryCSEM :: getTxZ( void ) const
{
	return tx_z;
}
//
void DataEntryCSEM :: setTxZ( double txZ )
{
	tx_z = txZ;
}
