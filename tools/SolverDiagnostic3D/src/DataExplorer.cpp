#include "DataFile.h"
//
//
int main( int argc, char ** argv )
{
	if( argc == 3 )
		DataFile* data_file_1 = new DataFile( argv[ 1 ], argv[ 2 ] );
	else
		std::cout << "usage: SolverDiagnostic3D <diagnostic_file_path> <math_box_bundle_path>" << std::endl;
}
