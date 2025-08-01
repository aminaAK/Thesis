#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

namespace CalcUtils {
	using namespace std;
	//________________________________________________________________________
	//					Точка на плоскости и ее координаты
	//________________________________________________________________________
	class Dot {
	public:
		double X , Y , Z;
		Dot( double x = 0 , double y = 0 , double z = 0 ) {
			X = x;
			Y = y;
			Z = z;
		}
		Dot( const Dot& d ) {
			X = d.X;
			Y = d.Y;
			Z = d.Z;
		}
		Dot operator+( const Dot& d )const {
			Dot R;
			R.X = X + d.X;
			R.Y = Y + d.Y;
			R.Z = Z + d.Z;
			return R;
		}
		Dot operator*( const double d )const {
			Dot R;
			R.X *= d;
			R.Y *= d;
			R.Z *= d;
			return R;
		}
		Dot operator-( const Dot& d )const {
			return *this + ( d * ( -1.0 ) );
		}
	};

	//________________________________________________________________________
	//			Двумерная вещественная матрица и операции с ней
	//________________________________________________________________________
	class Matrix {
		vector< vector<double> > f;
		long row;
		long column;
		double det( Matrix Arr ) const {
			if ( row != column || ( row*column ) == 0 ) {
				cout << "Warning! Cant calc Matrix det!" << endl;
				getchar( );
				abort( );
			}
			long i , j;
			double detInternal = 0;
			Matrix matr;
			if ( Arr.rows( ) == 1 ) {
				detInternal = Arr( 0 , 0 );
			}
			else {
				if ( Arr.rows( ) == 2 ) {
					detInternal = Arr( 0 , 0 ) * Arr( 1 , 1 ) - Arr( 0 , 1 ) * Arr( 1 , 0 );
				}
				else {
					matr.SetSize( Arr.rows( ) - 1 , Arr.rows( ) - 1 );
					for ( i = 0; i < Arr.rows( ); i++ ) {
						for ( j = 0; j < Arr.rows( ) - 1; j++ ) {
							if ( j < i ) {
								matr.GetFLink( ) [ j ] = Arr.GetFLink( ) [ j ];
								matr.GetFLink( ) [ j ].pop_back( ); // отрезаем лишнее
							}
							else {
								matr.GetFLink( ) [ j ] = Arr.GetFLink( ) [ j + 1 ];
								matr.GetFLink( ) [ j ].pop_back( ); // отрезаем лишнее
							}
						}
						detInternal += pow( -1. , double( i + j ) ) * det( matr ) * Arr( i , Arr.rows( ) - 1 );
					}
				}
			}
			return detInternal;
		}
	public:
		void Identity( ) {
			if ( row == column ) {
				for ( long i = 0; i < row; i++ ) {
					for ( long j = 0; j < column; j++ ) {
						f [ i ] [ j ] = long( i == j );
					}
				}
			}
		}
		void Zero( ) {
			for ( long i = 0; i < row; i++ ) {
				for ( long j = 0; j < column; j++ ) {
					f [ i ] [ j ] = 0.0;
				}
			}
		}
		vector< vector<double> > GetFCopy( ) const {
			return f;
		}
		vector< vector<double> >& GetFLink( ) {
			return f;
		}
		long rows( ) const {
			return row;
		}
		long columns( ) const {
			return column;
		}
		Matrix( const long ro = 1, const long col = 1 ) {
			row = ro;
			column = col;
			f.resize( ro );
			for ( int i = 0; i < col; ++i ){
				f [ i ].resize( col );
			}
		}
		Matrix( const Matrix& m ) {
			row = m.rows( );
			column = m.columns( );
			f = m.GetFCopy( );
		}
		void SetSize( long RowsSize , long ColumnsSize ) {
			if ( row != RowsSize || column != ColumnsSize ) {
				row = RowsSize;
				column = ColumnsSize;
				f.resize( row );
				for ( long i = 0; i < f.size( ); i++ ) {
					f [ i ].resize( column );
					for ( long j = 0; j < f [ i ].size( ); j++ ) {
						f [ i ] [ j ] = 0.0;
					}
				}
			}
		}
		Matrix(const vector< vector<double> >& arg){
			row = arg.size();
			column = 0;
			if(row > 0){
				column = arg[0].size();
			}
			f = arg;
		}
		double at( const long ro , const long col ) const {
			return f [ ro ] [ col ];
		}
		double& operator()( const long ro , const long col ) {
			return f [ ro ] [ col ];
		}

		vector<double> MulVr (vector<double> vec){  // calc vector M*v

            vector<double> out;
            double n;

            n = vec.size();
            out.resize(n,0);

            if ( n != column ) {
				cout << "Warning! Matrix and vector size failure!" << endl;
				getchar( );
				abort( );
			}

			for ( long i = 0; i < column; i++ ) {

				out[i] = 0; // just in case
				for ( long j = 0; j < column; j++ )
                    out[i] += vec[j]*f[i][j];
            }

            return out;
		}



		Matrix operator+( const Matrix& m ) const {
			if ( m.rows( ) != rows( ) || m.columns( ) != columns( ) ) {
				cout << "Warning! Matrix sizes failure!" << endl;
				getchar( );
				abort( );
			}
			Matrix MR( *this );
			for ( long i = 0; i < row; i++ ) {
				for ( long j = 0; j < column; j++ ) {
					MR( i , j ) += m.at( i , j );
				}
			}
			return MR;
		}
		Matrix operator-( const Matrix& m ) const {
			if ( m.rows( ) != rows( ) || m.columns( ) != columns( ) ) {
				cout << "Warning! Matrixes size failure!" << endl;
				getchar( );
				abort( );
			}
			Matrix MR( *this );
			for ( long i = 0; i < row; i++ ) {
				for ( long j = 0; j < column; j++ ) {
					MR( i , j ) -= m.at( i , j );
				}
			}
			return MR;
		}
		Matrix operator*( const double d ) const {
			Matrix MR( *this );
			for ( long i = 0; i < row; i++ ) {
				for ( long j = 0; j < column; j++ ) {
					MR( i , j ) *= d;
				}
			}
			return MR;
		}
		Matrix operator/( const double d ) const {
			if ( d == 0.0 ) {
				cout << "Warning! Division by zero!" << endl;
				getchar( );
				abort( );
			}
			Matrix MR( *this );
			for ( long i = 0; i < row; i++ ) {
				for ( long j = 0; j < column; j++ ) {
					MR( i , j ) /= d;
				}
			}
			return MR;
		}

		Matrix operator*( const Matrix& m ) const {
			if ( m.rows( ) != column ) {
				cout << "Warning! Matrixes size failure!" << endl;
				getchar( );
				abort( );
			}
			Matrix MR;
			MR.SetSize(row, m.columns());
			for ( long i = 0; i < row; i++ ) {
				for ( long j = 0; j < m.columns(); j++ ) {
					MR( i , j ) = 0;
					for ( long k = 0; k < column; k++ ) {
						MR( i , j ) += f [ i ] [ k ] * m.at( k , j );
					}
				}
			}
			return MR;
		}
		Matrix InvMatrix( ) const {
			if ( row != column || ( row*column ) == 0 ) {
				cout << "Warning! Cant reverse Matrix!" << endl;
				getchar( );
				abort( );
			}
			long size = row;
			Matrix inmat( *this ) , rm( *this );
			for ( long i = 0; i < size; i++ ) {
				for ( long j = 0; j < size; j++ ) {
					rm( i , j ) = long( i == j );
				}
			}
			vector<long> pos;
			pos.resize( size );
			for ( long i = 0; i < size; i++ ) {
				pos [ i ] = i;
			}
			for ( long i = 0; i < size; i++ ) {
				// начиная с i+1 строки удаляем нахер i-й столбец
				long bir = i;
				for ( long j = 0; j < size; j++ ) {
					if ( fabs( inmat( i , j ) ) > fabs( inmat( i , bir ) ) ) {
						bir = j;
					}
				}
				// переставляем столбцы i-й и bir-тый местами.
				double swa;
				{
					for ( long j = 0; j < size; j++ ) {
						swa = inmat( j , i );
						inmat( j , i ) = inmat( j , bir );
						inmat( j , bir ) = swa;
					}
				}
				{
					for ( long j = 0; j < size; j++ ) {
						swa = rm( j , i );
						rm( j , i ) = rm( j , bir );
						rm( j , bir ) = swa;
					}
				}
				{
					long iswa = pos [ i ];
					pos [ i ] = pos [ bir ];
					pos [ bir ] = iswa;
				}
				for ( long j = i + 1; j < size; j++ ) {
					double ct = inmat( j , i ) / inmat( i , i );
					for ( long k = 0; k < size; k++ ) {
						inmat( j , k ) -= ct * inmat( i , k );
						rm( j , k ) -= ct * rm( i , k );
					}
				}
			}
			// провели прямой ход преобразования
			// теперь обратный ход...
			for ( long i = size - 1; i > 0; i-- ) {
				for ( long j = i - 1; j > -1; j-- ) {
					double ct = inmat( j , i ) / inmat( i , i );
					for ( long k = 0; k < size; k++ ) {
						inmat( j , k ) -= ct * inmat( i , k );
						rm( j , k ) -= ct * rm( i , k );
					}
				}
			}
			for ( long i = 0; i < size; i++ ) {
				for ( long j = 0; j < size; j++ ) {
					rm( i , j ) /= inmat( i , i );
				}
				inmat( i , i ) = 1;
			}
			for ( long i = 0; i < size; i++ ) {
				bool done = false;
				for ( long j = i; ( j < size ) && ( !done ); j++ ) {
					if ( pos [ j ] == i ) {
						done = true;
						// надо поменять столбцы i, j местами, затем строки i, j местами
						double swa;
						long bir = j;
						{
							for ( long j = 0; j < size; j++ ) {
								swa = inmat( j , i );
								inmat( j , i ) = inmat( j , bir );
								inmat( j , bir ) = swa;
							}
						}
						{
							for ( long j = 0; j < size; j++ ) {
								swa = rm( j , i );
								rm( j , i ) = rm( j , bir );
								rm( j , bir ) = swa;
							}
						}
						{
							long iswa = pos [ i ];
							pos [ i ] = pos [ bir ];
							pos [ bir ] = iswa;
						}
						vector<double> swal;
						{
							swal = inmat.GetFLink( ) [ i ];
							inmat.GetFLink( ) [ i ] = inmat.GetFLink( ) [ bir ];
							inmat.GetFLink( ) [ bir ] = swal;
						}
						{
							swal = rm.GetFLink( ) [ i ];
							rm.GetFLink( ) [ i ] = rm.GetFLink( ) [ bir ];
							rm.GetFLink( ) [ bir ] = swal;
						}
					}
				}
			}
			return rm;
		}
		double det( ) const {
			Matrix Arr( *this );
			return det( Arr );
		}
	};
};
