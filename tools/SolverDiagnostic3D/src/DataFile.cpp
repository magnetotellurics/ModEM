#include "DataFile.h"

#include <stddef.h>
#include <iterator>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>

DataFile :: DataFile( std::string path )
:TextFile( path )
{
	loadValues();
	//
	createHtml();
}
DataFile :: ~DataFile( void ){}

int DataFile :: existPeriod( double target_value )
{
	for( int pid = 0; pid < (int) _periods.size(); pid++ )
		if( _periods[ pid ]->period == target_value )
			return pid;
	return -1;
}
void DataFile :: loadValues( void )
{
	int max_iteration = 0;
	//
	lx = ly = lz = hx = hy = hz = 0.0;
	//
	for( size_t lid = 0; lid < _lines.size(); lid++ )
	{
		std::vector< std::string > _line_parts = split( _lines[ lid ] );
		//
		Period * period = NULL;
		//
		double period_value = log10( atof( _line_parts[ 1 ].c_str() ) );
		//
		int period_index 	= atoi( _line_parts[ 2 ].c_str() );
		int polarization_index 	= atoi( _line_parts[ 3 ].c_str() );
		int iteration 	= atoi( _line_parts[ 4 ].c_str() );
		double residual = log10( atof( _line_parts[ 5 ].c_str() ) );
		//
		//
		lz = ( ( lz == 0.0 ) ? period_value : ( period_value < lz ) ? period_value : lz );
		hz = ( ( hz == 0.0 ) ? period_value : ( period_value > hz ) ? period_value : hz );
		//
		ly = ( ( ly == 0.0 ) ? residual : ( residual < ly ) ? residual : ly );
		hy = ( ( hy == 0.0 ) ? residual : ( residual > hy ) ? residual : hy );
		//
		lx = ( ( lx == 0.0 ) ? iteration : ( iteration < lx ) ? iteration : lx );
		hx = ( ( hx == 0.0 ) ? iteration : ( iteration > hx ) ? iteration : hx );
		//
		//
		if( atoi( _line_parts[ 4 ].c_str() ) > max_iteration )
			max_iteration = atoi( _line_parts[ 4 ].c_str() );
		//
		int founded_index = existPeriod( period_value );
		//
		if( founded_index == -1 )
			period = new Period;
		else
			period = _periods[ founded_index ];
		//
		period->index = period_index;
		period->period = period_value;
		//
		switch( polarization_index )
		{
			case 1:
				period->_pol_1_residuals.push_back( residual );
				break;
			case 2:
				period->_pol_2_residuals.push_back( residual );
				break;
			default:
				std::cout << "Wrong polarization: " << period_index << std::endl;
		}
		//
		if( founded_index == -1 )
			_periods.push_back( period );
		else
			_periods[ founded_index ] = period;
		//
	}
}
void DataFile :: createHtml( void )
{
	std::stringstream html_content;
	html_content << "<!DOCTYPE html>\n" <<
	"<html>\n" <<
	"<head>\n" <<
	"<meta charset='utf-8'>\n" <<
	"<title>" << path << "</title>\n" <<
	"<script src='build/mathbox-bundle.js'></script>\n" <<
	"<link rel='stylesheet' href='build/mathbox.css'>\n" <<
	"<meta name='viewport' content='initial-scale=1, maximum-scale=1'>\n" <<
	"</head>\n" <<
	"<body>\n" <<
	"<script>\n" <<
	"mathbox = mathBox({\n" <<
	"plugins: [ 'core', 'controls', 'cursor' ],\n" <<
	"controls: {\n" <<
	"klass: THREE.OrbitControls\n" <<
	"},\n" <<
	"});\n" <<
	"three = mathbox.three;\n" <<

	"three.camera.position.set( 5, 5, 10 );\n" <<
	"three.controls.maxDistance = 10;\n" <<
	"three.renderer.setClearColor( new THREE.Color(0xFAFAF8), 1.0 );\n" <<

	"view_pol_1 = mathbox.cartesian({\n" <<
	"range: [[0, 2], [0, 1], [0, 1]],\n" <<
	"scale: [3, 1, 1],\n" <<
	"position: [5, 0, 0],\n" <<
	"});\n" <<

	"view_pol_2 = mathbox.cartesian({\n" <<
	"range: [[0, 2], [0, 1], [0, 1]],\n" <<
	"scale: [3, 1, 1],\n" <<
	"position: [-5, 0, 0],\n" <<
	"});\n" <<

	"var dataMaximums = [" << hx << ", " << hy << ", " << lz << "];\n" <<
	"var dataMinimums = [" << lx << ", " << ly << ", " << hz << "];\n" <<
	"var dataRanges = [0,1,2,3].map(function(i){ return dataMaximums[i] - dataMinimums[i] })\n" <<
	"var dataScaledMinimums = [0,1,2,3].map(function(i){ return dataMinimums[i] / dataRanges[i] })\n" <<
	"var colors = {\n" <<
	"x: 0x000,   // red\n" <<
	"y: 0x000,   // yellow\n" <<
	"z: 0x000,   // blue\n" <<
	"xy: 0xFF851B,  // orange\n" <<
	"xz: 0xB10DC9,  // purple\n" <<
	"yz: 0x2ECC40,  // green\n" <<
	"xyz: 0x654321, // brown\n" <<
	"}\n" <<

	"function interpolate(lo, hi, n){\n" <<
	"n--; // go to end of range\n" <<
	"var vals = [];\n" <<
	"for (var i = 0; i <= n; i++){\n" <<
	"vals.push(Math.round(10 * (lo + (hi - lo)*(i/n)))/10);\n" <<
	"}\n" <<
	"return vals;\n" <<
	"}\n" <<

	"view_pol_1.scale({\n" <<
	"divide: 5,\n" <<
	"origin: [0,0,1,0],\n" <<
	"axis: 'x',\n" <<
	"}).text({\n" <<
	"live: false,\n" <<
	"data: interpolate(dataMinimums[0], dataMaximums[0], 5)\n" <<
	"}).label({\n" <<
	"color: colors.x,\n" <<
	"})\n" <<

	"view_pol_1.scale({\n" <<
	"divide: 3,\n" <<
	"origin: [0,0,1,0],\n" <<
	"axis: 'y',\n" <<
	"}).text({\n" <<
	"live: false,\n" <<
	"data: interpolate(dataMinimums[1], dataMaximums[1], 3)\n" <<
	"}).label({\n" <<
	"color: colors.y,\n" <<
	"offset: [-16, 0]\n" <<
	"})\n" <<

	"view_pol_1.scale({\n" <<
	"divide: 3,\n" <<
	"origin: [2,0,0,0],\n" <<
	"axis: 'z',\n" <<
	"}).text({\n" <<
	"live: false,\n" <<
	"data: interpolate(dataMinimums[2], dataMaximums[2], 3)\n" <<
	"}).label({\n" <<
	"color: colors.z,\n" <<
	"offset: [16, 0]\n" <<
	"})\n" <<

	"view_pol_2.scale({\n" <<
	"divide: 5,\n" <<
	"origin: [0,0,1,0],\n" <<
	"axis: 'x',\n" <<
	"}).text({\n" <<
	"live: false,\n" <<
	"data: interpolate(dataMinimums[0], dataMaximums[0], 5)\n" <<
	"}).label({\n" <<
	"color: colors.x,\n" <<
	"})\n" <<

	"view_pol_2.scale({\n" <<
	"divide: 3,\n" <<
	"origin: [0,0,1,0],\n" <<
	"axis: 'y',\n" <<
	"}).text({\n" <<
	"live: false,\n" <<
	"data: interpolate(dataMinimums[1], dataMaximums[1], 3)\n" <<
	"}).label({\n" <<
	"color: colors.y,\n" <<
	"offset: [-16, 0]\n" <<
	"})\n" <<

	"view_pol_2.scale({\n" <<
	"divide: 3,\n" <<
	"origin: [2,0,0,0],\n" <<
	"axis: 'z',\n" <<
	"}).text({\n" <<
	"live: false,\n" <<
	"data: interpolate(dataMinimums[2], dataMaximums[2], 3)\n" <<
	"}).label({\n" <<
	"color: colors.z,\n" <<
	"offset: [16, 0]\n" <<
	"})\n" <<

	"view_pol_1.grid({\n" <<
	"axes: 'xy',\n" <<
	"divideX: 3,\n" <<
	"divideY: 3,\n" <<
	"width: 5,\n" <<
	"})\n" <<
	".grid({\n" <<
	"axes: 'xz',\n" <<
	"divideX: 3,\n" <<
	"divideY: 3,\n" <<
	"width: 5,\n" <<
	"})\n" <<
	".grid({\n" <<
	"axes: 'yz',\n" <<
	"divideX: 3,\n" <<
	"divideY: 3,\n" <<
	"width: 5,\n" <<
	"})\n" <<

	"view_pol_2.grid({\n" <<
	"axes: 'xy',\n" <<
	"divideX: 3,\n" <<
	"divideY: 3,\n" <<
	"width: 5,\n" <<
	"})\n" <<
	".grid({\n" <<
	"axes: 'xz',\n" <<
	"divideX: 3,\n" <<
	"divideY: 3,\n" <<
	"width: 5,\n" <<
	"})\n" <<
	".grid({\n" <<
	"axes: 'yz',\n" <<
	"divideX: 3,\n" <<
	"divideY: 3,\n" <<
	"width: 5,\n" <<
	"})\n" <<

	"view_pol_1.array({\n" <<
	"data: [[dataMaximums[0]/2, dataMinimums[1] - 1, dataMaximums[2] - 1], [dataMinimums[0] - 20, dataMaximums[1] + 1, dataMinimums[2]/2], [dataMaximums[0] + 20, dataMinimums[1], dataMinimums[2]/2]],\n" <<
	"channels: 3,\n" <<
	"live: false,\n" <<
	"})\n" <<
	".transform({\n" <<
	"scale: dataRanges.slice(0,3).map(function(d,i){return i ? 1/d : 2/d}),\n" <<
	"position: dataScaledMinimums.slice(0,3).map(function(d,i){return i ? -d : -2*d}),\n" <<
	"})\n" <<
	".text({\n" <<
	"data: ['Iteration', 'LOG10 [Error]', 'LOG10[Period (s)]'],\n" <<
	"font: 'Helvetica', weight: 'bold', style: 'normal'\n" <<
	"})\n" <<
	".label({\n" <<
	"color: 0x000,\n" <<
	"size: 12\n" <<
	"});	\n" <<

	"view_pol_2.array({\n" <<
	"data: [[dataMaximums[0]/2, dataMinimums[1] - 1, dataMaximums[2] - 1], [dataMinimums[0] - 20, dataMaximums[1] + 1, dataMinimums[2]/2], [dataMaximums[0] + 20, dataMinimums[1], dataMinimums[2]/2]],\n" <<
	"channels: 3,\n" <<
	"live: false,\n" <<
	"})\n" <<
	".transform({\n" <<
	"scale: dataRanges.slice(0,3).map(function(d,i){return i ? 1/d : 2/d}),\n" <<
	"position: dataScaledMinimums.slice(0,3).map(function(d,i){return i ? -d : -2*d}),\n" <<
	"})\n" <<
	".text({\n" <<
	"data: ['Iteration', 'LOG10 [Error]', 'LOG10[Period (s)]'],\n" <<
	"font: 'Helvetica', weight: 'bold', style: 'normal'\n" <<
	"})\n" <<
	".label({\n" <<
	"color: 0x000,\n" <<
	"size: 12\n" <<
	"});\n";

	for( int pol = 1; pol <= 2; pol++ )
	{
		for( size_t pid = 0; pid < _periods.size(); pid++ )
		{
			html_content << "view_pol_" << pol << ".array({\n" <<
			"items: 1,\n" <<
			"channels: 3,\n" <<
			"live: false,\n" <<
			"id: 'data_" << pid+1 << "_" << pol << "',\n" <<
			"})\n" <<
			".swizzle({\n" <<
			"order: 'xyz'\n" <<
			"})\n" <<
			".transform({\n" <<
			"scale: dataRanges.slice(0,3).map(function(d,i){return i ? 1/d : 2/d}),\n" <<
			"position: dataScaledMinimums.slice(0,3).map(function(d,i){return i ? -d : -2*d}),\n" <<
			"})\n" <<
			".lerp({\n" <<
			"width: 40,\n" <<
			"})\n" <<
			".line({\n" <<
			"color: 0xFF4136,\n" <<
			"width: 10,\n" <<
			"})\n";
			//
			std::stringstream pol_stream_data;
			for( size_t poid = 0; poid < ( pol == 1 ? _periods[ pid ]->_pol_1_residuals.size() : _periods[ pid ]->_pol_2_residuals.size() ); poid++ )
				pol_stream_data << "[" << poid + 1 << ", " << ( pol == 1 ? _periods[ pid ]->_pol_1_residuals[ poid ] : _periods[ pid ]->_pol_2_residuals[ poid ] ) << ", " << _periods[ pid ]->period << "],\n"; 
			//
			html_content << "view_pol_" << pol << ".select('#data_" << pid+1 << "_" << pol << "').set('data', [\n" <<
			pol_stream_data.str() << "])\n";
		}
	}
	//
	html_content << "</script>\n" <<
	"</body>\n" <<
	"</html>\n";
	//
	std::stringstream output_file_name;
	output_file_name << "HTML_OUTPUT_" << name.substr( 0, name.find( "." ) ) << ".html";
	std::ofstream html_file;
	html_file.open ( output_file_name.str() );
	html_file << html_content.str();
	html_file.close();
}

