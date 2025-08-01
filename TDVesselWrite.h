#include "TaskData.h"

namespace TDVesselWrite {
	using namespace TaskData;
	using namespace Globals;

	void TDWrite(Zadacha& Z);

	inline int WriteCalcData ( Zadacha& Z) {

	    int IsDataWrite, n_wr;

	    if ( (T >= write_time)&&((last_write < 0.0000000001 )||(T - last_write > write_delta)))
            IsDataWrite = 1;
        else
            IsDataWrite = 0;

        if (IsDataWrite == 1){

            n_wr = int((T - write_time)/write_delta);
            last_write = n_wr*write_delta + write_time;

            //steps_skipped = 0;
            Globals::Nz++;

            TDWrite(Z); // writing time is inside
            //

        }


		return IsDataWrite;
	}
	void TDWriteCurrentTime(string filename) {
		ofstream fout;
		fout.open (filename, ofstream::app);
		fout << Globals::T << endl;
		fout.close();
	}
	void TDWriteCurrentTimeSteps(string filename) {
		ofstream fout;
		fout.open (filename, ofstream::out);
		fout << Globals::Nz << endl;
		fout.close();
	}
//------- Функции записи в файл задачно-зависимых переменных -------------------------
	void WriteResultLymph ( Derevo& Tr , Zadacha& Z ) {
		double q;
		ofstream fout;
		for ( int i = 0; i < Tr.Nbr; ++i ) {
			fout.open (Tr.B[i].TD.fnameTD[0], ofstream::app);
			for (int k = 0; k < Tr.B[i].pts; ++k) {
				fout << Tr.B[i].VB[0][k] * Tr.B[i].VB[1][k] << " ";
			}
			fout << endl;
			fout.close();
		}
	}
	void WriteResultBronch ( Derevo& Tr , Zadacha& Z ) { }
	void WriteResultArtery ( Derevo& Tr , Zadacha& Z ) { }
	void WriteResultVein ( Derevo& Tr , Zadacha& Z ) { }

	void WriteResultTD ( Derevo& Tr , Zadacha& Z ) {
		if ( Tr.ID == LYMPHATIC ) WriteResultLymph(Tr, Z);
		if ( ( Tr.ID == PULMART ) || ( Tr.ID == SYSART ) ) WriteResultArtery(Tr, Z);
		if ( ( Tr.ID == PULMVEN ) || ( Tr.ID == SYSVEN ) ) WriteResultVein(Tr, Z);
		if ( Tr.ID == BRONCHIAL ) WriteResultBronch(Tr, Z);
	}

//-------Подпрограмма записи в файл переменных s,u,p и прочих------------------------------
// может работать долго, потому что вычисляет давление
	void WriteResultV(Derevo& Tr , Zadacha& Z) {
		double p, Tcur;
		ofstream fout;
		vector<double> out(3);
        
        Tcur = T - Z.T_last_Hbeat;
		//int cur_track = 0; // если пишутся не все ветви
		for ( int i = 0; i < Tr.Nbr; i++ ) {
			/*if (Z.N_towrite) {
				if ( i == Z.Num_towrite[cur_track] ) {
					cur_track++;
				} else {
					continue;
				}
			}*/
			for ( int j = 0; j < Z.Cor; ++j ) { 				// s, u
				fout.open (Tr.B[i].TD.fnameVar[j], ofstream::app);
				//cout << Tr.B[i].TD.fnameVar[j] << endl;
				//cout << Tr.B[i].TD.fnameVar[Z.Cor] << endl;
				for (int k = 0; k < Tr.B[i].pts; ++k) {
					fout << Tr.B[i].VB[j][k] << endl;
					//cout << Tr.B[i].VB[1][Tr.B[i].pts-1] << ' '<<i<< " writen"<< endl;
				}
				fout.close();
			}

            //if ((i == 16) or (i == 17) or (i == 18) or (i == 23)){
            if ((i == 11) or (i == 12)){
                fout.open (Tr.B[i].TD.fnameVar[Z.Cor], ofstream::app); // Z.Cor = 2
                for (int k = 0; k < Tr.B[i].pts; k++) { 			// p
                    out =  Tr.B[i].MusclePump(Tr.B[i].VB[0][k],Tr.B[i].VB[1][k],Tcur);
                    p = out[0]/1333.2; // mmHg
                    fout << p << endl;
                }
            }
            else {
                fout.open (Tr.B[i].TD.fnameVar[Z.Cor], ofstream::app); // Z.Cor = 2
                for (int k = 0; k < Tr.B[i].pts; k++) {             // p
                    out =  Tr.B[i].URSOB(Tr.B[i].VB[0][k],Tr.B[i].VB[1][k]);
                    p = out[0]/1333.2; // mmHg
                    fout << p << endl;
                }
                    
            }
			fout.close();

			fout.open (Tr.B[i].TD.fnameVar[Z.Cor + 1], ofstream::app);
			for (int k = 0; k < Tr.B[i].pts; k++) { 			// q
				p = Tr.B[i].VB[0][k]*Tr.B[i].VB[1][k];
				fout << p << endl;
			}
			fout.close();

			fout.open (Tr.B[i].TD.fnameVar[Z.Cor + 2], ofstream::app);
			for (int k = 0; k < Tr.B[i].pts; k++) { 			// pave
				p = Tr.B[i].TD.Pave_next/1333.2;
				fout << p << endl;
			}
			fout.close();

        }
		//WriteResultTD(Tr , Z);
		TDGlobals::isFirstTime = 0;
	}

	void WriteBranchResultForGNUPLOT (Derevo& Tr , Zadacha& Z) {}

	void TDWrite(Zadacha& Z) {
	//-------------------------- MATLAB output ---------------------------------
		if ( Globals::use_MATLAB_out ) {
			for ( int i = 0; i < Z.Ntr; ++i ) {
				WriteResultV(TreeLst[i], Z);
			}
			//if ( Z.useLungs ) {
			//	WriteResultV(BronchTree, Z);
			//}
			TDWriteCurrentTimeSteps(trim(SharedDirectory) + "tsteps.tres");
			//TDWriteCurrentTime(trim(SharedDirectory) + "time.tres");
		}
	//-------------------------- GNUPLOT output ---------------------------------
		/*if ( Globals::use_GNUPLOT_out ) {
			for ( int i = 0; i < Z.Ntr; ++i ) {
				WriteBranchResultForGNUPLOT (TreeLst[i] , Z );
			}
			if ( Z.useLungs ) {
				WriteBranchResultForGNUPLOT ( BronchTree, Z );
			}
		}*/
	}
//===========================================================================================
//  сохранение текущих значений для использования в дальнейших вычислениях
//===========================================================================================
	void SaveDump(Derevo& Tr , Zadacha& Z) {
		string fnDump;
		int LResult;

		LResult = MAKEDIRQQ( trim( SharedDirectory ) + trim(Tr.dirname) );
		LResult = MAKEDIRQQ( trim( SharedDirectory ) + trim(Tr.dirname) + slash + "result" );
		fnDump = trim(SharedDirectory) + trim(Tr.dirname) + slash + "result" + slash + "tree" + trim( Achar( Ichar( "0" ) + Tr.ID ) )  +  ".dmp";
        
		cout << "Dump file : " << fnDump << endl;
		ofstream fout( fnDump , ofstream::out );

		fout <<  Tr.ID << endl;
		for ( int i = 0; i < Tr.Nbr; ++i ) {
			fout <<  Tr.B[ i ].ID << " " <<  Tr.B[ i ].pts << endl;;		// номер ветви, количество точек разбиения
			for ( int j = 0; j < Z.Cor; ++j ) {
				for (int k = 0; k < Tr.B[i].pts; ++k) {
					fout << Tr.B[i].VB[j][k] << " ";
				}
				fout << endl;
			}
			for ( int j = 0; j < Z.Cor; ++j ) {
				for (int k = 0; k < Tr.B[i].pts; ++k) {
					fout << Tr.B[i].VBO[j][k] << " ";
				}
				fout << endl;
			}
		}
		fout.close();
	}

	void ConvertResults(Derevo& Tr , Zadacha& Z) {}

	/*void WriteFFR(Derevo& Tr , Zadacha& Z) {

	    double Pd,Pa, FFR;
	    long id_next, addBr;

	    ofstream fou;
        fou.open(Tr.FFRfile);



        for (long i = 0; i < Z.Nsten; i++){

            if (Tr.B[Z.IDsten[i] - 1].stenType == 1){

                id_next = (*(*Tr.B[Z.IDsten[i] - 1].Kn2).Bou[0]).ID;

                Pd = Tr.B[id_next - 1].TD.Pave_next;
                Pa = Tr.B[0].TD.Pave_next;
                FFR = Pd/Pa;

                if (id_next <= (Tr.NbrL + 2))
                    addBr = 2;
                else
                    addBr = Tr.NbrL + 2;

                fou << "BrID " <<   (Z.IDsten[i] - addBr)  << endl;
                fou << "FFR = " <<  FFR << endl;
                fou << endl;

            }

        }
        cout << "FFR calculated" << endl;

        fou.close();

	}*/
};
