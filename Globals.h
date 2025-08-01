
#pragma once
#include <string>
#include "TreeDataIndepend.h"
#include "TDGlobals.h"

namespace Globals {

	using namespace std;
	using namespace TreeDataIndepend;
	using namespace TDGlobals;

//#define WIN
#ifdef WIN
	string slash = "\\";
#endif
#ifndef WIN
	string slash = "/";
#endif

	inline string trim( string s ) {
		bool rep = true;
		while ( s.length( ) > 0 && rep ) {
			rep = false;
			if ( s [ 0 ] == ' ' ) {
				s = s.substr( 1 , s.length( ) - 1 );
				rep = true;
			}
			else {
				if ( s [ s.length( ) - 1 ] == ' ' ) {
					s = s.substr( 0 , s.length( ) - 1 );
					rep = true;
				}
			}
		}
		return s;
	}
	inline string Adjustl( string s ) {
		bool rep = true;
		while ( s.length( ) > 0 && rep ) {
			rep = false;
			if ( s [ 0 ] == ' ' ) {
				s = s.substr( 1 , s.length( ) - 1 );
				rep = true;
			}
		}
		return s;
	}
	inline string Adjustr( string s ) {
		bool rep = true;
		while ( s.length( ) > 0 && rep ) {
			rep = false;
			if ( s [ s.length( ) - 1 ] == ' ' ) {
				s = s.substr( 0 , s.length( ) - 1 );
				rep = true;
			}
		}
		return s;
	}
	inline long Ichar( string s ) {
		long R = 0;
		if ( s.length( ) ) {
			R = long( s [ 0 ] );
		}
		return R;
	}
	inline string Achar( long l ) {
		char rb [ 2 ];
		rb [ 1 ] = '\0';
		rb [ 0 ] = char( l );
		return string( rb );
	}
	inline long MAKEDIRQQ( string s ) {
		return system( ( "mkdir " + s ).c_str( ) );
	}
	inline long ClearFile( string s ) {
#ifdef WIN
		return system( ( "del " + s ).c_str( ) );
#endif
#ifndef WIN
		return system( ( "rm " + s ).c_str( ) );
#endif
	}
	void ReadFile( string fn, double *x, int n){
		string tmp;
		ifstream fin( fn, ifstream::in );
		for (int i = 0; i < n; i++){
			fin >> x[i];
			getline(fin , tmp);
		}
		fin.close();
	}
    
	const static double PI = 3.1415926;
	const static double GRAV = 982.2;  // ��/c**2
	const static double ANGLE = PI / 2;
	const static double SMALLPAR = 0.05;
	const static long Z_VESSELS = 111;

	//** ���� ��������� ������� � �����
	const static long BIFURCATION = 1;	// inner, junction
	const static long FLOW = 2;				// boundary
    const static long VALVE = 22;    
	const static long INCOMING = 1;
	const static long OUTGOING = -1;

	const static long LYMPHATIC = 0; // ADDED BY RUFINA
	const static long PULMART = 3; //1; //4
	const static long PULMVEN = 4; //2; //2
	const static long SYSART = 1; //3; //3
	const static long SYSVEN = 2; //4; //1
	const static long BRONCHIAL = 5;

	const static long ZadachaPars = 30;
	const static long GlobalPars = 8;
	const static long LymphPars = 6;
	const static long BRANCH_ORDER = 4;

	double T , dt , Smax , time_calc , write_time , write_delta , write_dump , last_write , last_write_dump;
	long Metod , NKO , Nz , Norg , chng;
	long use_MATLAB_out , use_GNUPLOT_out;

	vector<Derevo> TreeLst;
	Derevo BronchTree;
	MultiKnotList MKnots;


	double start_time , last_elapsed_time;
	long finAlvExtP;
	string SharedDirectory;
}
