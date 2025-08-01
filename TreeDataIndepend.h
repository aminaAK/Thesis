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
	//					Uzel(Knot)
	//	C - point on the coordinate plane
	//	ID	- knot number
	// Nou	- number of outgoing branches
	// Nin	- number of ingoing branches
	// IG	- type of knot
	// TD	- task-depended knot parameter set
	// determines the type of entry of a branch into a knot
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
	//					Branch
	//	ID				- branch number
	// pts				- number of points on the branch
	// dx				- branch step
	// len				- branch lenght
	//	width			- branch width
	//	Kn1, Kn2	- pointers to knots adjacent to the edge; Kn1 - beginning, Kn2 - end
	//	VB				- the value of the vector of variables on the branch
	// VBO				- the value of the vector of variables on the branch on the previous step
	// TD				- task-depended branch parameter set
	// gravity           - does gravity work ( 1 - yes, 0 - no)
	//________________________________________________________________________
	class Vetv {
	public:
		long ID , pts , myTreeID , InvertPoints, group, stenType;
		Uzel *Kn1 , *Kn2;
		vector< vector<double> > VB, VBO;
		double len , width , dx;
		TDVetvParameter TD;


		// S - cross-section, u - velocity
		vector<double> URSOB(double S,double u) const {

		    vector<double> output(3);       // [0] - pressure, [1] - derivative Pr, [2] - zvuk
		    double rho_w,rho;                   //density of a wall

		    rho_w = 1.0;
		    rho = 1.0;
            
            //cout<<S<<' '<<ID<<" URSOB"<<endl;
		    if (S > TD.S0){

                output[0] = TD.c*TD.c*rho_w*(exp(S/TD.S0 - 1) - 1);
                output[1] = TD.c*TD.c*rho_w*exp(S/TD.S0 - 1)/TD.S0;
                //cout << output[0] << ' '<< " P "<<endl;

		    } else {

                output[0] = TD.c*TD.c*rho_w*log(abs(S/TD.S0));
                output[1] = TD.c*TD.c*rho_w/S;
                //cout << output[0] << ' '<< " P "<<endl;
		    }

		    output[2] = sqrt(S*output[1]/rho);

		    if ( output[2] < 0.00001){
                output[2] = 0.00001;
                cout << "Warning! PWV is too small" << endl;
		    }

			return output;
		}
        vector<double> MusclePump(double S,double u, double Tcur) const {

            vector<double> pump(3);       // [0] - pressure, [1] - derivative Pr, [2] - zvuk
            double rho_w,rho, Pmax, Tstep, Pcor, PI,eps,a,PWV;                   //density of a wall
            PI = 3.14159;
            //Tstep = 1.0;
            Tstep = 1.33;
            rho_w = 1.0;
            rho = 1.0;
            Pmax = 100*1333.2;
            //Pmax = 350*1333.2; // Pmax = 10 kPa 1mmHg = 133.32 Pa 308 329 350
            eps = 1.00;
            a = 3; //T = 2.00
            //eps = 0.32;// 1.33- Amina added Pmax*(-np.sin(np.pi*2.8*Tcur/Tstep+np.pi/3))
            //a = 0.755; //1.338
            //Pcor = Pmax*(-sin(2.5*PI*Tcur/Tstep+PI*0.29)+0.9);
            //eps = 0.27;// 1.00 0.23
            //a = 0.6; //1.00
            //PWV = 100;
            //cout<<S<<' '<<ID<<" pump"<<endl;
            if (abs(Tcur - a) <= eps) Pcor = Pmax*exp(-(eps*eps)/(eps*eps - pow(abs(Tcur - a),2)));
            //if (abs(Tcur - a) <= eps) Pcor = Pmax*exp(-(eps*eps)/(eps*eps - pow(abs(Tcur - a),2)))*2*eps;
            else Pcor = 0;
            //cout<<Tcur<<endl;
            if (S > TD.S0){
                pump[0] = TD.c*TD.c*rho_w*(exp(S/TD.S0 - 1) - 1) - Pcor;
                //S = TD.S0*(1 + log(1 + (pump[0] - Pcor)/(TD.c*TD.c)));
                pump[1] = TD.c*TD.c*rho_w*exp(S/TD.S0 - 1)/TD.S0;
                //pump[0] =  PWV * PWV *rho_w*(exp(S/TD.S0 - 1) - 1) + Pcor;
                //pump[1] =  PWV * PWV *rho_w*exp(S/TD.S0 - 1)/TD.S0;
                //cout << pump[0]/1333.2 << ' ' << S << " P S"<<endl;

            } else {

                //pump[0] = PWV *PWV *rho_w*log(abs(S/TD.S0)) + Pcor;
                //pump[1] = PWV *PWV *rho_w/S;
                pump[0] = TD.c*TD.c*rho_w*log(abs(S/TD.S0)) - Pcor;
                //S = TD.S0*(exp((pump[0]-Pcor))/(TD.c*TD.c));
                pump[1] = TD.c*TD.c*rho_w/S;
                //cout << pump[0]/1333.2 << ' ' << S << " P S"<<endl;
            }

            pump[2] = sqrt(S*pump[1]/rho);

            if ( pump[2] < 0.00001){
                pump[2] = 0.00001;
                cout << "Warning! PWV is too small" << endl;
            }

            return pump;
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
