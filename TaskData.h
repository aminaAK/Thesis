#include <iostream>
#include <fstream>
#include "Globals.h"
//#include <omp.h>

namespace TaskData {
	using namespace Globals;

	class Zadacha {
	public:
		long Cor , ID , Ntr , useDump , Nmatter , useLungs , useExtImpactToLungs , useOrgans , useFlowAveraging , calcMax , useGravity , useSaveDump , writePmax;
		double TQaver0; // ������ ���������� ��������������� ������������ ������
		double TQaver; // ������ ���������� ������������ ������
		double Kur , P_barorec , h_period_new , T_last_Hbeat , h_period_curr , N_heart_cycles , Pnapol; //� ����� : �������� ��������� , ����� ����������� ��������� ������ , ����� �������� ������� , ���� �������� ���������� �������
		double Tc_st , S_st; //stenosis

		double corPres, vasCoef, Qratio, diamCoef, Pmean, Pveins;
		double HR, SV, Res, PWV, PowM, NL, NR;
		double Tsl, Tsr, Tfl, Tfr, TsH, TfH;  //inside heart period

		long Nave_cycles , N_towrite, Nsten;

		long* Num_towrite; // ������ ������ ��� ������������ ������ ���� ����������� � ������� �����������
		long* IDsten; // stenoses IDs
		long Debug; // 1 or 0

		void IdentifyTask( ) {
			//����������� ���������� ������
			string filename = Globals::SharedDirectory + "ini" + slash + "zadacha.ini";
            double tmp[ZadachaPars];
			ReadFile(filename, tmp, ZadachaPars);
			ID = tmp[0];
			Cor = tmp[1];				//���������� ����������
			Ntr = tmp[2];				//���������� ��������
			Kur = tmp[3];
			Nmatter = tmp[4];			//���������� ����������� �������
			useDump = tmp[5];
			useSaveDump = tmp[6];
			useLungs = tmp[7];
			useExtImpactToLungs = tmp[8];
			useOrgans = tmp[9];
			useFlowAveraging = tmp[10];
			TQaver0 = tmp[11];          // overwritten later
			TQaver = tmp[12];           // overwritten later
			calcMax = tmp[13];
			diamCoef = tmp[14];
			vasCoef = tmp[15];
			corPres = tmp[16];
			Res = tmp[17];
			PWV = tmp[18];
			PowM = tmp[19];
			Qratio = tmp[20];
			HR = tmp[21];
			SV = tmp[22];
			Pmean = tmp[23];
			Pveins = tmp[24];
			Debug = tmp[25];

			// ������ ������ ��� ������������ ������ ���� ����������� � ������� �����������
			filename = Globals::SharedDirectory + "ini" + slash + "towrite.ini";
			fstream fin;
			fin.open( filename , ifstream::in );
			fin >> N_towrite; // amount of tracked vessels
			getline(fin , filename);
			if ( N_towrite > 0 ) {
				Num_towrite = new long [ N_towrite ];
				for ( long i = 0; i < N_towrite; i++ ) {
					fin >> Num_towrite [ i ];
				}
			}
			fin.close( );

			this->Nave_cycles = 1;
			this->T_last_Hbeat = 0.0;
			this->h_period_curr = 60.0/this->HR;
			this->h_period_new = this->h_period_curr;
			this->TQaver0 = this->h_period_curr;
			this->TQaver = this->h_period_curr;
			this->N_heart_cycles = 0;
			this->writePmax = 0;
		}

		void PrintTask() {
			cout << "Task initializing" << endl;
			cout << "ID: BLOOD" << endl;
			//cout << "Number of equations: " << this->Cor << endl;
			//cout << "Number of networks: " << this->Ntr << endl;
			cout << "Courant number: " << this->Kur << endl;
			//cout << "Number of substances: " << this->Nmatter << endl;

			if ( this->useDump == 1 ) {
				cout << "Continued simulations" << endl;
			}
			else {
				cout << "New simulations" << endl;
			}
			if ( this->useLungs == 1 ) {
				cout << "Lungs: Yes" << endl;
			}
			else {
				cout << "Lungs: No" << endl;
			}
			cout << "======================================" << endl;
		}
	};

	//�������� ���������� ����������
	inline void LoadGlobals( ) {
		ifstream fin( "/Users/amina/Documents/thesis/thesis/paths.ini" , ifstream::in );
		fin >> Globals::SharedDirectory;
		fin.close( );

		string filename = Globals::SharedDirectory + "ini" + slash + "globals.ini";
		double tmp[GlobalPars];
		ReadFile(filename, tmp, GlobalPars);
		//cout << filename << endl;
		Globals::time_calc = tmp[0];
		//cout << tmp[0] << endl;
		Globals::NKO = tmp[1];
		Globals::Metod = tmp[2];
		Globals::write_time = tmp[3];
		Globals::write_delta = tmp[4];
		Globals::write_dump = tmp[5];
		Globals::use_MATLAB_out = tmp[6];
		Globals::use_GNUPLOT_out = tmp[7];

		cout << "MainDir: " << Globals::SharedDirectory << endl;
		cout << "Total time (sec): " << Globals::time_calc << endl;
		//cout << "Start results writing (write_time): " << Globals::write_time << endl;
		//cout << "Steps to skip while writing (NKO): " << Globals::NKO << endl;
		//cout << "Numerical method ID (Metod): " << Globals::Metod << endl;

		/*if ( Globals::use_MATLAB_out > 0 ) {
			cout << "MATLAB output enabled" << endl;
		}
		else {
			cout << "MATLAB output disabled" << endl;
		}

		if ( Globals::use_GNUPLOT_out > 0 ) {
			cout << "GNUPLOT output enabled" << endl;
		}
		else {
			cout << "GNUPLOT output disabled" << endl;
		}*/


		//Globals::start_time = omp_get_wtime( );
		Globals::last_elapsed_time = 0.;
	}
};
