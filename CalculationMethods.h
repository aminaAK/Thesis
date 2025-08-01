//#include "TreeInit.h"
#include "Predictor.h"
#include "Granuslov.h"
//#include "omp.h"

namespace CalculationMethods {
	using namespace Predictor;
	using namespace Granuslov;


	// this includes autoregulation. We do not use autoregulation for coronary, hence beta = 0
	void VesselFlowAveraging(Zadacha& Z, Vetv& br, long Nave){

	    vector<double> out(3);
	    double P, QQ0, PP0, alpha, beta;



        if (Nave < 0){
             //averaging over first cycle
            for (long i = 0; i < br.pts; i++ ){

                out = br.URSOB(br.VB[0][i],br.VB[1][i]);
                P = out[0];

                br.TD.Pave_prev = br.TD.Pave_prev + P*dt/(br.pts*Z.h_period_curr);
                br.TD.Qave_prev = br.TD.Qave_prev + br.VB[0][i]*br.VB[1][i]*dt/(br.pts*Z.h_period_curr);
            }

        }else{

            if (Nave > Z.Nave_cycles){

                if (Z.Nave_cycles > 1){
                    br.TD.Qave_prev = br.TD.Qave_next;
                    br.TD.Pave_prev = br.TD.Pave_next;
                }

                br.TD.Qave_next = br.TD.Qave_curr;
                br.TD.Pave_next = br.TD.Pave_curr;

                br.TD.Qave_curr = 0;
                br.TD.Pave_curr = 0;

                if (br.TD.Qave_prev > 0){

                    QQ0 = abs(br.TD.Qave_next/br.TD.Qave_prev);

                    if (QQ0 > 1) br.TD.cnewQ = br.TD.c*exp((3/8)*log(QQ0));
                    else br.TD.cnewQ = br.TD.c*exp((1/4)*log(QQ0));
                }

                if (br.TD.Pave_prev > 0){

                    PP0 = abs(br.TD.Pave_next/br.TD.Pave_prev);

                    if (PP0 > 1) br.TD.cnewP = br.TD.c*exp((1/2)*log(PP0));
                    else br.TD.cnewP = br.TD.c*exp((1/2)*log(PP0));         // same for now
                }
            }


            alpha  = 1.0;
            beta = 0;

            //br.TD.c = br.TD.c + beta*(br.TD.cnewP - br.TD.c)*dt/(Z.h_period_curr*alpha);

            for (long i = 0; i < br.pts; i++ ){
                out = br.URSOB(br.VB[0][i],br.VB[1][i]);
                P = out[0];
                //cout<<br.ID<< ' '<<br.VB[1][i]<<' '<< i<<endl;

                br.TD.Pave_curr = br.TD.Pave_curr + P*dt/(br.pts*Z.h_period_curr);
                br.TD.Qave_curr = br.TD.Qave_curr + br.VB[0][i]*br.VB[1][i]*dt/(br.pts*Z.h_period_curr);
            }
        }
    
	}
    void VesselFlowPumpAveraging(Zadacha& Z, Vetv& br, long Nave){

        vector<double> out(3);
        double P, QQ0, PP0, alpha, beta, Tcur;
        Tcur = T - Z.T_last_Hbeat;
        if (Nave < 0){
             //averaging over first cycle
            for (long i = 0; i < br.pts; i++ ){
                
                out = br.MusclePump(br.VB[0][i],br.VB[1][i], Tcur);
                P = out[0];

                br.TD.Pave_prev = br.TD.Pave_prev + P*dt/(br.pts*Z.h_period_curr);
                br.TD.Qave_prev = br.TD.Qave_prev + br.VB[0][i]*br.VB[1][i]*dt/(br.pts*Z.h_period_curr);
            }

        }else{

            if (Nave > Z.Nave_cycles){

                if (Z.Nave_cycles > 1){
                    br.TD.Qave_prev = br.TD.Qave_next;
                    br.TD.Pave_prev = br.TD.Pave_next;
                }

                br.TD.Qave_next = br.TD.Qave_curr;
                br.TD.Pave_next = br.TD.Pave_curr;

                br.TD.Qave_curr = 0;
                br.TD.Pave_curr = 0;

                if (br.TD.Qave_prev > 0){

                    QQ0 = abs(br.TD.Qave_next/br.TD.Qave_prev);

                    if (QQ0 > 1) br.TD.cnewQ = br.TD.c*exp((3/8)*log(QQ0));
                    else br.TD.cnewQ = br.TD.c*exp((1/4)*log(QQ0));
                }

                if (br.TD.Pave_prev > 0){

                    PP0 = abs(br.TD.Pave_next/br.TD.Pave_prev);

                    if (PP0 > 1) br.TD.cnewP = br.TD.c*exp((1/2)*log(PP0));
                    else br.TD.cnewP = br.TD.c*exp((1/2)*log(PP0));         // same for now
                }
            }


            alpha  = 1.0;
            beta = 0;

            //br.TD.c = br.TD.c + beta*(br.TD.cnewP - br.TD.c)*dt/(Z.h_period_curr*alpha);

            for (long i = 0; i < br.pts; i++ ){
                
                out = br.MusclePump(br.VB[0][i],br.VB[1][i], Tcur);
                P = out[0];

                br.TD.Pave_curr = br.TD.Pave_curr + P*dt/(br.pts*Z.h_period_curr);
                br.TD.Qave_curr = br.TD.Qave_curr + br.VB[0][i]*br.VB[1][i]*dt/(br.pts*Z.h_period_curr);
            }
        }
        //cout << br.ID <<"VVesselFlowPumpAveraging"<<endl;
    }



	inline void CalcVesselVertexes( Zadacha& Z ){
        
        //double a;
        for (long i = 0; i < Z.Ntr; i++ ){
            for  (long j = 0; j < TreeLst[i].Nkn; j++ ){
                 /*if (j == 5) {
                     cout << " CalcVesselVertexes "<<TreeLst[i].Nkn<<endl;
                     cout<<TreeLst[0].B[0].VB[1][TreeLst[0].B[0].pts-1]<<" before Grtoch  "<<endl;
                     cout<<endl;
                 }*/
                 Grtoch(Z,TreeLst[i],TreeLst[i].K[j]);
                /*if (j == 5) {
                     cout<<TreeLst[0].B[0].VB[1][TreeLst[0].B[0].pts-1]<<" after Grtoch  "<<endl;
                     cout<<endl;
                 }*/
                 //a = TreeLst[0].B[0].VB[1][TreeLst[0].B[0].pts-1];
                 //cout<<a<<" a "<< j<<endl;
             }
             //cout<<a<<" a"<<endl;
             //TreeLst[0].B[0].VB[1][TreeLst[0].B[0].pts-1] = a;
             // end of for loop  j < Tr.Nbr
             //cout<<TreeLst[0].B[0].VB[1][TreeLst[0].B[0].pts-1]<<" END of for loop  j < Tr.Nbr CalcVesselVertexes "<<endl;
	    } // end of for loop  i < Z.Nkn
        //cout<<TreeLst[0].B[0].VB[1][TreeLst[0].B[0].pts-1]<<" END  CalcVesselVertexes "<<endl;
	}

	inline void CalcEdges( Zadacha& Z ) {

	    long Nave, isCorPr, n;

	    double Tcur;

	    //calculating averages, all averages are linked to heart period
        if ((Z.useFlowAveraging == 1)&&(Z.N_heart_cycles >= 2))
            if (Z.N_heart_cycles == 2) Nave = -1;
            else Nave = Z.N_heart_cycles - 2;


	    for (long i = 0; i < Z.Ntr; i++ ){

            #pragma omp parallel for
             for  (long j = 0; j < TreeLst[i].Nbr; j++ ){
                 //cout<<j<<"for"<<endl;
                 if ((j == 11) or (j == 12)){
                 //if ((j == 16) or (j == 17) or (j == 18) or (j == 23)){
                     //cout<<j<<"pump"<<endl;
                     
                     Hybrid2ndpump(Z,TreeLst[i],TreeLst[i].B[j]);
                     
                 }
                 else{
                     /*if(j==0){
                         cout<<TreeLst[i].B[j].VB[1][TreeLst[i].B[j].pts-1]<<  " CalcEdges "<<TreeLst[i].B[j].ID<<endl;
                         cout<<TreeLst[0].B[0].VB[1][TreeLst[0].B[0].pts-1]<<  " CalcEdges "<<TreeLst[0].B[0].ID<<endl;
                         cout << endl;
                     }*/
                     Hybrid2nd(Z,TreeLst[i],TreeLst[i].B[j]);
                     //cout<<j<<"comm"<<endl;
                 }

                 if ((Z.useFlowAveraging == 1)&&(Z.N_heart_cycles >= 2)){
                     //if ((j == 16) or (j == 17) or (j == 18) or (j == 23)){
                     if ((j == 11) or (j == 12)){
                         //cout<<j<<"pump"<<endl;
                         VesselFlowPumpAveraging(Z,TreeLst[i].B[j], Nave);
                         
                     }
                     else {
                         VesselFlowAveraging(Z,TreeLst[i].B[j],Nave);
                         //cout<<j<<"comm"<<endl;
                     }
                 }
             }   // end of for loop  j < Tr.Nbr
	    } // end of for loop  i < Z.Ntr

        if ((Z.useFlowAveraging == 1)&&(Z.N_heart_cycles >= 2))
            if (Nave > Z.Nave_cycles)
                Z.Nave_cycles = Nave;

	}

	inline double GetDtAboutSMax(Zadacha& Z, double ASmax ){

        return Z.Kur/ASmax;

	}

};
