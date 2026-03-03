//########################################################################
//	Модуль определения задача-независимых элементов дерева или сети (графа)
// Uzel   - узел (точка ветвления, граничная точка)
// Vetv   - ветвь (ребро)
// Derevo - по сути сеть (граф), дерево - частный случай, содержит списки
// ветвей, узлов, внещних воздействий и т.п., а также имена файлов с
// параметрами ветвей и узлов
//########################################################################
#include "TreeDataDepend.h"
#include <vector>
//#include <math.h>
namespace TreeDataIndepend {

	using namespace TreeDataDepend;
	//________________________________________________________________________
	//					Узел
	//	С		- Точка на координатной плоскости
	//	ID	- номер узла
	// Nou	- кол-во выходящих ветвей
	// Nin	- кол-во входящих ветвей
	// IG	- тип узла (ветвление, вход, выход и т.п.)
	// TD	- задача-зависимый набор параметров узла
	// определяет тип вхождения ветви в узел (входящая или исходящая)
	//________________________________________________________________________
	class Vetv;
	class Uzel {
	public:
		Dot C;
		long ID , Nou , Nin , IG;
		vector<Vetv*> Bou;
		vector<Vetv*> Bin;
		//bool LYMPHNODE;
		TDUzelParameter TD;

	};

	//________________________________________________________________________
	//					Ветвь
	//	ID				- номер ветви
	// pts				- количество точек на ветви
	// dx				- шаг на ветви
	// len				- длина ветви
	//	width			- ширина ветви
	//	Kn1, Kn2	- указатели на прилегающие к ребру узлы; Kn1 - начало, Kn2 - конец
	//	VB				- значение вектора переменных на ветви
	// VBO				- значение вектора переменных на ветви на предыдущем шаге
	// TD				- задача-зависимый набор параметров ветви
	// gravity           - работает ли гравитация (1 или 0)
	//________________________________________________________________________
	class Vetv {
	public:
		long ID , pts , myTreeID , InvertPoints, group, stenType;
		Uzel *Kn1 , *Kn2;
		vector< vector<double> > VB, VBO;
		double len , width , dx;
		TDVetvParameter TD;


		// S - cross-section, u - velocity
		vector<double> URSOB(double S,double u, long ID, double T, long x) const {

		    vector<double> output(3);       // [0] - pressure, [1] - derivative Pr, [2] - zvuk
		    double rho_w,rho;                   //density of a wall
            double  Pmax, Tstep, Pcor, PI,eps,a,PWV, Snew, b, epsi;                   //density of a wall
            PI = 3.14159;
            //Tstep = 1.0;
            Tstep = 1.33;
            rho_w = 1.0;
            rho = 1.0;
            //Pmax = 1*1333.2;
            eps = 6.00;
            a = 5;
            b=4;
            epsi=5;
            
            if (abs(x - b) <= epsi) Pmax = 200*1333.2*exp(-(epsi*epsi)/(epsi*epsi - pow(abs(x - b),2)));
            else Pmax = 0;
            
            if (abs(T - a) <= eps) Pcor = Pmax*exp(-(eps*eps)/(eps*eps - pow(abs(T - a),2)));
            else Pcor = 0;
            //Pcor = 20*1333.2 * (sin(2 * PI * T / 10)+1);
            
		    rho_w = 1.0;
		    rho = 1.0;
            if (ID == 2){
                //cout<< "s "<< S << " s0 " << TD.S0<<endl;
                if (S > TD.S0){
                    output[0] = TD.c*TD.c*rho_w*(exp(S/TD.S0 - 1) - 1) + Pcor;
                    output[1] = TD.c*TD.c*rho_w*exp(S/TD.S0 - 1)/TD.S0;
                    
                } else {
                    output[0] = TD.c*TD.c*rho_w*log(abs(S/TD.S0)) + Pcor;
                    output[1] = TD.c*TD.c*rho_w/S;
                    
                }
            }
            else {
                //            if (S <= 0) {
                //                cout << "URSOB: negative S = " << S << " at branch " << ID << ", point " << x << endl;
                //                S = 1e-6;
                //            }
                
                if (S > TD.S0){
                    
                    output[0] = TD.c*TD.c*rho_w*(exp(S/TD.S0 - 1) - 1);
                    output[1] = TD.c*TD.c*rho_w*exp(S/TD.S0 - 1)/TD.S0;
                    
                } else {
                    
                    output[0] = TD.c*TD.c*rho_w*log(abs(S/TD.S0));
                    output[1] = TD.c*TD.c*rho_w/S;
                }
            }
            
            if ( output[1] < 0){
                cout<<"dP/dS "<< S << " at branch " << ID << ", point " << x << endl;
            }
            if (S <= 0) {
                cout << "URSOB: negative S = " << S << " at branch " << ID << ", point " << x << endl;
                S = 1e-6;
            }
            
                
            output[2] = sqrt(S*output[1]/rho);
                
                
            if ( output[2] < 0.00001){
                output[2] = 0.00001;
                cout << "Warning! PWV is too small" << endl;
            }
        

			return output;
		}


		vector<double> FPRCH(double S,double u, long i) const {

		    vector<double> FPR(2);

		    FPR[0] = 0;
		    FPR[1] = (-8.0)*(3.1415926)*(0.04)*u/S;

		    return FPR;
		}
	};

	void printVetv ( Vetv& B ) {
		cout << "ID = " << B.ID << " pts = " << B.pts << " myTreeID = " << B.myTreeID << endl;
		cout << "len = " << B.len << " width = " << B.width << " dx = " << B.dx << endl;
		cout << "Kn1 = " << B.Kn1->ID << " Kn2 = " << B.Kn2->ID << endl;
		cout << "--------------------------------------" << endl;
	}
	void printUzel ( Uzel& K ) {
			cout << "ID = " << K.ID << " IG = " << K.IG <<  endl;
			cout << "Nin = " << K.Nin << " Nou = " << K.Nou << endl;
			cout << "Bin : ";
			for (int i = 0; i < K.Bin.size(); ++i) cout << K.Bin[i]->ID << " ";
			cout << endl <<  "Bou : ";
			for (int i = 0; i < K.Bou.size(); ++i) cout << K.Bou[i]->ID << " ";
			cout << endl << "--------------------------------------" << endl;
		}

	//________________________________________________________________________
	//					Дерево
	//	ID	- номер дерева
	// Nbr	- количество ветвей
	// Nkn	- количество узлов
	//	B	- список ветвей дерева
	//	K	- список узлов дерева
	//________________________________________________________________________
	class Derevo {
	public:
		long ID , Nbr , Nkn , Nimp , Norg, NbrL, NbrR, NknL, NknR, Nlad, Nlcx, Nrca, iLad, iLcx, iRca, iCa;
		string treefilename , dirname;
		string knotfilename , TDknotfilename , TDknotExternalFilename;
		string branchfilename , TDbranchfilename;
		string TDImpactfilename;
		string inputdata;
		string FFRfile;
		Vetv *B;
		Uzel *K;
		vector<TDExternalImpact> I;
	};

	// Группа узлов одного дерева для мультиузла
	class MultiKnotTreeGroup {
	public:
		long TrID;
		long sz;
		vector<Uzel> KntLst;
	};

	// Мультиузел  - осуществляет склейку нескольких концевых узлов разных графов (или одного и того же графа) в один
	class MultiKnot {
	public:
		long sz;
		long szin;
		long szou;
		vector<MultiKnotTreeGroup> GrpLst;
	};

	// Список мультиузлов
	class MultiKnotList {
	public:
		long sz;
		vector<MultiKnot> Lst;
		string multiKnotsfilename;
	};
};
