#include "TreeInit.h"
#include "tuple"
#include "Subroutines.h"

namespace Predictor {
	using namespace TreeInit;
	using namespace Subroutines;


	inline vector<double> FV (vector<double> V, Zadacha& Z , double P, Vetv& br){

         double rho;
         vector<double> F(2);

         rho = 1.0;

         F[0] = V[0]*V[1];
         F[1] = V[1]*V[1]/2 + P/rho;

         return F;

	}

    inline vector<double> SumV(vector<double> A, vector<double> B)
    {

        int n;
        vector<double> Res;

        if ( A.size()!= B.size() )
        {
            cout << "Warning! Vectors size failure!" << endl;
            getchar( );
            abort( );
        }

        n = A.size();
        Res.resize(n);

        for ( int i = 0; i < n; i++ )
            Res[i] = A[i] + B[i];

        return Res;

	}

	tuple<Matrix,Matrix> SH_Hybrid2nd(vector<double> V, Zadacha& Z, Vetv& br){


            double zvuk, P, Pr, gamma;
            vector<double> out(3);
            Matrix B_om(Z.Cor,Z.Cor), D_om(Z.Cor,Z.Cor), B_h(Z.Cor,Z.Cor), D_h(Z.Cor,Z.Cor), OM(Z.Cor,Z.Cor), OMO(Z.Cor,Z.Cor);
            vector<double> S(Z.Cor);

            out = br.URSOB(V[0],V[1]);
            P = out[0];
            Pr = out[1];
            zvuk = out[2];
            //cout << P << " SH_Hybrid2nd"<<endl;

            tie(OMO,OM) = OMEGAB(V, Z, zvuk, P, br);

            gamma = 1.0;
            tie(S, B_h, D_h) = SZH_Hybr(V, Z, zvuk, P, gamma, br);

            if ((abs(S[0])/br.dx) > Smax)
                Smax = abs(S[0])/br.dx;
            if ((abs(S[1])/br.dx) > Smax)
                Smax = abs(S[1])/br.dx;

            B_om = (OMO*B_h)*OM;
            D_om = (OMO*D_h)*OM;

            return {B_om,D_om};

	}
    tuple<Matrix,Matrix> SH_Hybrid2ndpump(vector<double> V, Zadacha& Z, Vetv& br){


            double zvuk, P, Pr, gamma, Tcur;
            vector<double> out(3);
            Matrix B_om(Z.Cor,Z.Cor), D_om(Z.Cor,Z.Cor), B_h(Z.Cor,Z.Cor), D_h(Z.Cor,Z.Cor), OM(Z.Cor,Z.Cor), OMO(Z.Cor,Z.Cor);
            vector<double> S(Z.Cor);
            Tcur = T - Z.T_last_Hbeat;
            out = br.MusclePump(V[0],V[1], Tcur);
            
            P = out[0];
            Pr = out[1];
            zvuk = out[2];
            //cout << P << " SH_Hybrid2ndpump"<<endl;

            tie(OMO,OM) = OMEGAB(V, Z, zvuk, P, br);

            gamma = 1.0;
            tie(S, B_h, D_h) = SZH_Hybr(V, Z, zvuk, P, gamma, br);

            if ((abs(S[0])/br.dx) > Smax)
                Smax = abs(S[0])/br.dx;
            if ((abs(S[1])/br.dx) > Smax)
                Smax = abs(S[1])/br.dx;

            B_om = (OMO*B_h)*OM;
            D_om = (OMO*D_h)*OM;

            return {B_om,D_om};

    }




	void Hybrid2nd (Zadacha& Z , Derevo& Tr , Vetv& br){
        // Hybrid 1-2 order scheme

        long Np;
        double zvuk, P, Pr, Eps;

        vector< vector<double> >  Vint(Z.Cor,vector<double>(br.pts)), V1(Z.Cor,vector<double>(br.pts));

        vector<double> FM(Z.Cor),FP(Z.Cor),UM(Z.Cor),UP(Z.Cor),VM(Z.Cor),VP(Z.Cor),FPR(Z.Cor),V(Z.Cor),V_cur(Z.Cor);

        vector<double> out(3);

        Matrix SL(Z.Cor,Z.Cor), SP(Z.Cor,Z.Cor), OM(Z.Cor,Z.Cor), OMO(Z.Cor,Z.Cor), B_om(Z.Cor,Z.Cor), D_om(Z.Cor,Z.Cor);
        /*if (br.ID ==1){ cout << "Hybrid2nd"<<endl;
            cout<< br.VBO[1][br.pts-1] << " VBO "<< br.VB[1][br.pts-1] << " VB "<< br.ID <<endl;
        }*/
        Np = br.pts;
        br.VBO = br.VB;
        /*if (br.ID ==1){
            cout<< br.VBO[1][br.pts-1] << " VBO "<< br.VB[1][br.pts-1] << " VB "<< br.ID <<endl;
            cout <<endl;
        }*/
        

        // ---------------------- 1st cycle -------------------

        Vint = br.VB;

        for ( int i = 1; i < (Np - 1); i++ ){

            VM[0] = br.VB[0][i-1];
            VM[1] = br.VB[1][i-1];

            out = br.URSOB(VM[0],VM[1]);
            P = out[0];

            FM = FV(VM,Z,P,br);

            VP[0] = br.VB[0][i+1];
            VP[1] = br.VB[1][i+1];

            out = br.URSOB(VP[0],VP[1]);
            P = out[0];
            //cout << P << " Hybrid2nd"<<endl;


            FP = FV(VP,Z,P,br);

            Vint[0][i] = br.VB[0][i] - (dt/(2*br.dx))*(FP[0] - FM[0]);
            Vint[1][i] = br.VB[1][i] - (dt/(2*br.dx))*(FP[1] - FM[1]);
        }

        // ---------------------- 2nd cycle -------------------

        /* //testing
        tie(V, B_om, D_om) = Test();

        cout << B_om(0,0) << "  " << B_om(0,1) << endl;
        cout << B_om(1,0) << "  " << B_om(1,1) << endl;

        cout << D_om(0,0) << "  " << D_om(0,1) << endl;
        cout << D_om(1,0) << "  " << D_om(1,1) << endl;


        cout << V[0] << endl;
        cout << V[1] << endl;
        cout << V[1] << endl; */


        for ( int i = 1; i < (Np - 1); i++ ){

            V[0] =  Vint[0][i];
            V[1] =  Vint[1][i];

            V_cur[0] = (br.VBO[0][i] + br.VBO[0][i+1])/2;
            V_cur[1] = (br.VBO[1][i] + br.VBO[1][i+1])/2;

            tie(B_om, D_om) = SH_Hybrid2nd(V_cur, Z,br);

            V_cur[0] = br.VBO[0][i+1] - br.VBO[0][i];
            V_cur[1] = br.VBO[1][i+1] - br.VBO[1][i];

            V = SumV(V,B_om.MulVr(V_cur));

            V_cur[0] = Vint[0][i+1] + Vint[0][i] - (br.VBO[0][i+1] + br.VBO[0][i]);
            V_cur[1] = Vint[1][i+1] + Vint[1][i] - (br.VBO[1][i+1] + br.VBO[1][i]);

            V = SumV(V,D_om.MulVr(V_cur));

            V_cur[0] = (br.VBO[0][i] + br.VBO[0][i-1])/2;
            V_cur[1] = (br.VBO[1][i] + br.VBO[1][i-1])/2;

            tie(B_om, D_om) = SH_Hybrid2nd(V_cur, Z,br);

            V_cur[0] = br.VBO[0][i-1] - br.VBO[0][i];   // take (-1) into account
            V_cur[1] = br.VBO[1][i-1] - br.VBO[1][i];   // take (-1) into account

            V = SumV(V,B_om.MulVr(V_cur));

            V_cur[0] = (Vint[0][i+1] + Vint[0][i] - (br.VBO[0][i+1] + br.VBO[0][i]))*(-1);
            V_cur[1] = (Vint[1][i+1] + Vint[1][i] - (br.VBO[1][i+1] + br.VBO[1][i]))*(-1);

            V = SumV(V,D_om.MulVr(V_cur));

            V_cur[0] = br.VBO[0][i];
            V_cur[1] = br.VBO[1][i];

            FPR = br.FPRCH(V_cur[0],V_cur[1],i);
            FPR[0] = dt*FPR[0];
            FPR[1] = dt*FPR[1];

            V = SumV(V,FPR);

            br.VB[0][i] = V[0];
            br.VB[1][i] = V[1];
        }
        /*if (br.ID ==  1){
            cout << "afterHybrid2nd"<<endl;
            cout<< br.VBO[1][br.pts-1] << " VBO "<< br.VB[1][br.pts-1] << " VB "<< br.ID <<endl;
            
            cout<< br.VBO[1][br.pts-1] << " VBO "<< br.VB[1][br.pts-1] << " VB "<< br.ID <<endl;
            cout <<endl;
        }*/

        // boundary, not necessary ?, existed in previous version
        br.VB[0][0] = br.VBO[0][0];
        br.VB[1][0] = br.VBO[1][0];

        br.VB[0][Np-1] = br.VBO[0][Np-1];
        br.VB[1][Np-1] = br.VBO[1][Np-1];
        
        //br.VB[1][Np-1] = br.VB[1][Np-2];
        //cout<<br.VB[1][Np-1]<<" Predictor2 "<< br.ID<<endl;
	}

    void Hybrid2ndpump (Zadacha& Z , Derevo& Tr , Vetv& br){
        // Hybrid 1-2 order scheme

        long Np;
        double zvuk, P, Pr, Eps,Tcur;

        vector< vector<double> >  Vint(Z.Cor,vector<double>(br.pts)), V1(Z.Cor,vector<double>(br.pts));

        vector<double> FM(Z.Cor),FP(Z.Cor),UM(Z.Cor),UP(Z.Cor),VM(Z.Cor),VP(Z.Cor),FPR(Z.Cor),V(Z.Cor),V_cur(Z.Cor);

        vector<double> out(3);

        Matrix SL(Z.Cor,Z.Cor), SP(Z.Cor,Z.Cor), OM(Z.Cor,Z.Cor), OMO(Z.Cor,Z.Cor), B_om(Z.Cor,Z.Cor), D_om(Z.Cor,Z.Cor);

        Tcur = T;
        Np = br.pts;
        br.VBO = br.VB;

        // ---------------------- 1st cycle -------------------

        Vint = br.VB;

        for ( int i = 1; i < (Np - 1); i++ ){

            VM[0] = br.VB[0][i-1]; //S
            VM[1] = br.VB[1][i-1];

            out = br.MusclePump(VM[0],VM[1],Tcur);
            P = out[0];

            FM = FV(VM,Z,P,br);

            VP[0] = br.VB[0][i+1];
            VP[1] = br.VB[1][i+1];

            out = br.MusclePump(VP[0],VP[1], Tcur);
            P = out[0];
            //cout << P << " Hybrid2ndpump"<<endl;


            FP = FV(VP,Z,P,br);

            Vint[0][i] = br.VB[0][i] - (dt/(2*br.dx))*(FP[0] - FM[0]);
            Vint[1][i] = br.VB[1][i] - (dt/(2*br.dx))*(FP[1] - FM[1]);
        }

        // ---------------------- 2nd cycle -------------------

        /* //testing
        tie(V, B_om, D_om) = Test();

        cout << B_om(0,0) << "  " << B_om(0,1) << endl;
        cout << B_om(1,0) << "  " << B_om(1,1) << endl;

        cout << D_om(0,0) << "  " << D_om(0,1) << endl;
        cout << D_om(1,0) << "  " << D_om(1,1) << endl;


        cout << V[0] << endl;
        cout << V[1] << endl;
        cout << V[1] << endl; */


        for ( int i = 1; i < (Np - 1); i++ ){

            V[0] =  Vint[0][i];
            V[1] =  Vint[1][i];

            V_cur[0] = (br.VBO[0][i] + br.VBO[0][i+1])/2;
            V_cur[1] = (br.VBO[1][i] + br.VBO[1][i+1])/2;

            tie(B_om, D_om) = SH_Hybrid2ndpump(V_cur, Z,br);

            V_cur[0] = br.VBO[0][i+1] - br.VBO[0][i];
            V_cur[1] = br.VBO[1][i+1] - br.VBO[1][i];

            V = SumV(V,B_om.MulVr(V_cur));

            V_cur[0] = Vint[0][i+1] + Vint[0][i] - (br.VBO[0][i+1] + br.VBO[0][i]);
            V_cur[1] = Vint[1][i+1] + Vint[1][i] - (br.VBO[1][i+1] + br.VBO[1][i]);

            V = SumV(V,D_om.MulVr(V_cur));

            V_cur[0] = (br.VBO[0][i] + br.VBO[0][i-1])/2;
            V_cur[1] = (br.VBO[1][i] + br.VBO[1][i-1])/2;

            tie(B_om, D_om) = SH_Hybrid2ndpump(V_cur, Z,br);

            V_cur[0] = br.VBO[0][i-1] - br.VBO[0][i];   // take (-1) into account
            V_cur[1] = br.VBO[1][i-1] - br.VBO[1][i];   // take (-1) into account

            V = SumV(V,B_om.MulVr(V_cur));

            V_cur[0] = (Vint[0][i+1] + Vint[0][i] - (br.VBO[0][i+1] + br.VBO[0][i]))*(-1);
            V_cur[1] = (Vint[1][i+1] + Vint[1][i] - (br.VBO[1][i+1] + br.VBO[1][i]))*(-1);

            V = SumV(V,D_om.MulVr(V_cur));

            V_cur[0] = br.VBO[0][i];
            V_cur[1] = br.VBO[1][i];

            FPR = br.FPRCH(V_cur[0],V_cur[1],i);
            FPR[0] = dt*FPR[0];
            FPR[1] = dt*FPR[1];

            V = SumV(V,FPR);

            br.VB[0][i] = V[0];
            br.VB[1][i] = V[1];
        }


        // boundary, not necessary ?, existed in previous version
        br.VB[0][0] = br.VBO[0][0];
        br.VB[1][0] = br.VBO[1][0];

        br.VB[0][Np-1] = br.VBO[0][Np-1];
        br.VB[1][Np-1] = br.VBO[1][Np-1];
    }
};

