//#include "TreeInit.h"
#include "Predictor.h"
#include "Granuslov.h"


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

                out = br.URSOB(br.VB[0][i],br.VB[1][i],br.ID, T, i);
                P = out[0];

                br.TD.Pave_prev = br.TD.Pave_prev + P*dt/(br.pts);
                br.TD.Qave_prev = br.TD.Qave_prev + br.VB[0][i]*br.VB[1][i]*dt/(br.pts);
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

            br.TD.c = br.TD.c + beta*(br.TD.cnewP - br.TD.c)*dt/(Z.h_period_curr*alpha);

            for (long i = 0; i < br.pts; i++ ){
                out = br.URSOB(br.VB[0][i],br.VB[1][i], br.ID, T, i);
                P = out[0];

                br.TD.Pave_curr = br.TD.Pave_curr + P*dt/(br.pts);
                br.TD.Qave_curr = br.TD.Qave_curr + br.VB[0][i]*br.VB[1][i]*dt/(br.pts);
            }
        }

	}


	inline void CalcVesselVertexes( Zadacha& Z ){


         for (long i = 0; i < Z.Ntr; i++ ){

            #pragma omp parallel for
            for  (long j = 0; j < TreeLst[i].Nkn; j++ ){

                 Grtoch(Z,TreeLst[i],TreeLst[i].K[j]);

            }   // end of for loop  j < Tr.Nbr
	    } // end of for loop  i < Z.Nkn

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

                 Hybrid2nd(Z,TreeLst[i],TreeLst[i].B[j]);

                 if ((Z.useFlowAveraging == 1)&&(Z.N_heart_cycles >= 2))
                    VesselFlowAveraging(Z,TreeLst[i].B[j],Nave);

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
