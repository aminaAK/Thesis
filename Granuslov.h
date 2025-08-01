#include "TDGranuslov.h"

namespace Granuslov {

    using namespace CalcUtils;
    using namespace TaskData;
    using namespace TDGranuslov;


    /*
    Calculate function F in a X for Newton's method, Bernoulli
    */

    vector<double> FUZ_Ber (Zadacha& Z, vector<double> Xs, vector<double> Xu, vector<Vetv*> brin, vector<Vetv*> brou, long N, vector<double> alf,  vector<double> bet){

        long szin, szou, j;

        vector<double> e(N), F(N), P(N), out(3);

        double rho;

        rho = 1.0;

        szin = brin.size();
        szou = brou.size();

        for (long i = 0; i < szin; i++){

            e[i] = 1;
            F[i] = 0;
            //tie(alf[i],bet[i]) = IncomingCompatibilityCoeffs(Z, *(brin[i]));

            out = (*(brin[i])).URSOB(Xs[i],Xu[i]);
            P[i] = out[0];
        }

        for (long i = 0; i < szou; i++){

            j = i + szin;
            e[j] = - 1;
            F[j] = 0;
            //tie(alf[j],bet[j]) = OutgoingCompatibilityCoeffs(Z, *(brou[i]));

            out = (*(brou[i])).URSOB(Xs[j],Xu[j]);
            P[j] = out[0];
        }

        F[0] = e[0]*(alf[0]*Xs[0] + bet[0])*Xs[0];

        for (long i = 1; i < N; i++){

            F[i] = P[i-1]/rho - P[i]/rho + 0.5*pow(alf[i-1]*Xs[i-1] + bet[i-1],2) - 0.5*pow(alf[i]*Xs[i] + bet[i],2);
            F[0] = F[0] + e[i]*(alf[i]*Xs[i] + bet[i])*Xs[i];
        }

        return F;
    }


   
    vector<double> FUZ_Ber_pump (Zadacha& Z, vector<double> Xs, vector<double> Xu, vector<Vetv*> brin, vector<Vetv*> brou, long N, vector<double> alf,  vector<double> bet){

        long szin, szou, j;

        vector<double> e(N), F(N), P(N), out(3);

        double rho, Tcur;

        rho = 1.0;
        Tcur = T-Z.T_last_Hbeat;
        szin = brin.size();
        szou = brou.size();

        for (long i = 0; i < szin; i++){

            e[i] = 1;
            F[i] = 0;
            //tie(alf[i],bet[i]) = IncomingCompatibilityCoeffs(Z, *(brin[i]));

            out = (*(brin[i])).MusclePump(Xs[i],Xu[i], Tcur);
            P[i] = out[0];
        }

        for (long i = 0; i < szou; i++){

            j = i + szin;
            e[j] = - 1;
            F[j] = 0;
            //tie(alf[j],bet[j]) = OutgoingCompatibilityCoeffs(Z, *(brou[i]));

            out = (*(brou[i])).MusclePump(Xs[j],Xu[j], Tcur);
            P[j] = out[0];
        }

        F[0] = e[0]*(alf[0]*Xs[0] + bet[0])*Xs[0];

        for (long i = 1; i < N; i++){

            F[i] = P[i-1]/rho - P[i]/rho + 0.5*pow(alf[i-1]*Xs[i-1] + bet[i-1],2) - 0.5*pow(alf[i]*Xs[i] + bet[i],2);
            F[0] = F[0] + e[i]*(alf[i]*Xs[i] + bet[i])*Xs[i];
        }

        return F;
    }


/*
Calculate Jacobian J in a X for Newton's method, Bernoulli
*/

    Matrix TDYacobian_Ber(Zadacha& Z, vector<double> Xs, vector<double> Xu, vector<Vetv*> brin, vector<Vetv*> brou, long N, vector<double> alf,  vector<double> bet){

        long szin, szou, j;

        vector<double> e(N), F(N), P(N), Ps(N), out(3);

        Matrix YAC(N,N);

        double rho;

        rho = 1.0;

        szin = brin.size();
        szou = brou.size();

        for (long i = 0; i < szin; i++){

            e[i] = 1;
            F[i] = 0;
            //tie(alf[i],bet[i]) = IncomingCompatibilityCoeffs(Z, *(brin[i]));

            out = (*(brin[i])).URSOB(Xs[i],Xu[i]);
            P[i] = out[0];
            Ps[i] = out[1];
        }

        for (long i = 0; i < szou; i++){

            j = i + szin;
            e[j] = - 1;
            F[j] = 0;
            //tie(alf[j],bet[j]) = OutgoingCompatibilityCoeffs(Z, *(brou[i]));

            out = (*(brou[i])).URSOB(Xs[j],Xu[j]);
            P[j] = out[0];
            Ps[j] = out[1];
        }

        YAC.Zero();

        for (long i = 0; i < N; i++)
            YAC(0,i) = e[i]*(2*alf[i]*Xs[i]+ bet[i]);

        for (long i = 1; i < N; i++){

            YAC(i,i-1) = Ps[i-1]/rho + (alf[i-1]*Xs[i-1] + bet[i-1])*alf[i-1];
            YAC(i,i) = (-1)*(Ps[i]/rho + (alf[i]*Xs[i] + bet[i])*alf[i]);
        }

        return YAC;
    }
    Matrix TDYacobian_Ber_pump(Zadacha& Z, vector<double> Xs, vector<double> Xu, vector<Vetv*> brin, vector<Vetv*> brou, long N, vector<double> alf,  vector<double> bet){

        long szin, szou, j;

        vector<double> e(N), F(N), P(N), Ps(N), out(3);

        Matrix YAC(N,N);

        double rho, Tcur;
        Tcur = T - Z.T_last_Hbeat;
        rho = 1.0;

        szin = brin.size();
        szou = brou.size();

        for (long i = 0; i < szin; i++){

            e[i] = 1;
            F[i] = 0;
            //tie(alf[i],bet[i]) = IncomingCompatibilityCoeffs(Z, *(brin[i]));

            out = (*(brin[i])).MusclePump(Xs[i],Xu[i],Tcur);
            P[i] = out[0];
            Ps[i] = out[1];
        }

        for (long i = 0; i < szou; i++){

            j = i + szin;
            e[j] = - 1;
            F[j] = 0;
            //tie(alf[j],bet[j]) = OutgoingCompatibilityCoeffs(Z, *(brou[i]));

            out = (*(brou[i])).MusclePump(Xs[j],Xu[j], Tcur);
            P[j] = out[0];
            Ps[j] = out[1];
        }

        YAC.Zero();

        for (long i = 0; i < N; i++)
            YAC(0,i) = e[i]*(2*alf[i]*Xs[i]+ bet[i]);

        for (long i = 1; i < N; i++){

            YAC(i,i-1) = Ps[i-1]/rho + (alf[i-1]*Xs[i-1] + bet[i-1])*alf[i-1];
            YAC(i,i) = (-1)*(Ps[i]/rho + (alf[i]*Xs[i] + bet[i])*alf[i]);
        }

        return YAC;
    }


    vector<double> FUZ_Puaz(Zadacha& Z, vector<double> Xs, vector<double> Xu, vector<Vetv*> brin, vector<Vetv*> brou, long N,  vector<double> alf,  vector<double> bet, double R){

        long szin, szou, j;

        vector<double> e(N), F(N), P(N), out(3);

        double rho;

        rho = 1.0;

        szin = brin.size();
        szou = brou.size();

        for (long i = 0; i < szin; i++){

            e[i] = 1;
            F[i] = 0;
            //tie(alf[i],bet[i]) = IncomingCompatibilityCoeffs(Z, *(brin[i]));

            out = (*(brin[i])).URSOB(Xs[i],Xu[i]);
            P[i] = out[0];
        }

        for (long i = 0; i < szou; i++){

            j = i + szin;
            e[j] = - 1;
            F[j] = 0;
            //tie(alf[j],bet[j]) = OutgoingCompatibilityCoeffs(Z, *(brou[i]));

            out = (*(brou[i])).URSOB(Xs[j],Xu[j]);
            P[j] = out[0];
        }
        double vlv;
        vlv = VLVFunc(Xs[0], alf[0],bet[0]);
        F[0] = e[0]*(alf[0]*Xs[0] + bet[0])*Xs[0] + e[1]*(alf[1]*Xs[1] + bet[1])*Xs[1];
        //F[1] = (P[0] - P[1]) - vlv*R*Xs[0]*(alf[0]*Xs[0] + bet[0]);
        F[1] = vlv*(P[0] - P[1]) - R*Xs[0]*(alf[0]*Xs[0] + bet[0]);
        
        return F;
    }
    vector<double> FUZ_Puaz_Pump(Zadacha& Z, vector<double> Xs, vector<double> Xu, vector<Vetv*> brin, vector<Vetv*> brou, long N,  vector<double> alf,  vector<double> bet, double R){

        long szin, szou, j;

        vector<double> e(N), F(N), P(N), out(3);

        double rho, Tcur;
        Tcur = T - Z.T_last_Hbeat;
        rho = 1.0;

        szin = brin.size();
        szou = brou.size();

        for (long i = 0; i < szin; i++){

            e[i] = 1;
            F[i] = 0;
            //tie(alf[i],bet[i]) = IncomingCompatibilityCoeffs(Z, *(brin[i]));

            out = (*(brin[i])).MusclePump(Xs[i],Xu[i],Tcur);
            P[i] = out[0];
        }

        for (long i = 0; i < szou; i++){

            j = i + szin;
            e[j] = - 1;
            F[j] = 0;
            //tie(alf[j],bet[j]) = OutgoingCompatibilityCoeffs(Z, *(brou[i]));

            out = (*(brou[i])).MusclePump(Xs[j],Xu[j],Tcur);
            P[j] = out[0];
        }
        double vlv;
        vlv = VLVFunc(Xs[0], alf[0],bet[0]);
        F[0] = e[0]*(alf[0]*Xs[0] + bet[0])*Xs[0] + e[1]*(alf[1]*Xs[1] + bet[1])*Xs[1];
        //F[1] = (P[0] - P[1]) - vlv*R*Xs[0]*(alf[0]*Xs[0] + bet[0]);
        F[1] = vlv*(P[0] - P[1]) - R*Xs[0]*(alf[0]*Xs[0] + bet[0]);
        return F;
    }

    /*
    Calculate Jacobian J in a X for Newton's method, Puazel
    */

    Matrix TDYacobian_Puaz(Zadacha& Z, vector<double> Xs, vector<double> Xu, vector<Vetv*> brin, vector<Vetv*> brou, long N,  vector<double> alf,  vector<double> bet, double R){

        long szin, szou, j;

        vector<double> e(N), F(N), P(N), Ps(N), out(3);

        Matrix YAC(N,N);

        double rho;

        rho = 1.0;

        szin = brin.size();
        szou = brou.size();

        for (long i = 0; i < szin; i++){

            e[i] = 1;
            F[i] = 0;
            //tie(alf[i],bet[i]) = IncomingCompatibilityCoeffs(Z, *(brin[i]));

            out = (*(brin[i])).URSOB(Xs[i],Xu[i]);
            P[i] = out[0];
            Ps[i] = out[1];
        }

        for (long i = 0; i < szou; i++){

            j = i + szin;
            e[j] = - 1;
            F[j] = 0;
            //tie(alf[j],bet[j]) = OutgoingCompatibilityCoeffs(Z, *(brou[i]));

            out = (*(brou[i])).URSOB(Xs[j],Xu[j]);
            P[j] = out[0];
            Ps[j] = out[1];
        }

        YAC.Zero();
        double vlv, vlvd;

        vlv = VLVFunc(Xs[0], alf[0], bet[0]);
        //vlvd = VLVFuncDer(Xs[0], alf[0], bet[0]);
        for (long i = 0; i < N; i++)
            YAC(0,i) = e[i]*(2*alf[i]*Xs[i]+ bet[i]);

        for (long i = 1; i < N; i++){

            /*YAC(i,i-1) = Ps[i-1] - vlv*(R*(2*alf[i-1]*Xs[i-1]+ bet[i-1])) + (P[i-1] - P[i]) - vlvd * R*Xs[i-1]*(alf[i-1]*Xs[i-1] + bet[i-1]);
            YAC(i,i) = (-1)*Ps[i];*/
            YAC(i,i-1) = vlv*Ps[i-1] -(R*(2*alf[i-1]*Xs[i-1]+ bet[i-1])) + vlv*(1-vlv) * (P[i-1] - P[i]) -  R*Xs[i-1]*(alf[i-1]*Xs[i-1] + bet[i-1]);
            YAC(i,i) = (-1)*Ps[i]*vlv;
        }

        return YAC;
    }

    Matrix TDYacobian_Puaz_Pump(Zadacha& Z, vector<double> Xs, vector<double> Xu, vector<Vetv*> brin, vector<Vetv*> brou, long N,  vector<double> alf,  vector<double> bet, double R){

        long szin, szou, j;

        vector<double> e(N), F(N), P(N), Ps(N), out(3);

        Matrix YAC(N,N);

        double rho, Tcur;
        
        Tcur = T - Z.T_last_Hbeat;

        rho = 1.0;

        szin = brin.size();
        szou = brou.size();

        for (long i = 0; i < szin; i++){

            e[i] = 1;
            F[i] = 0;
            //tie(alf[i],bet[i]) = IncomingCompatibilityCoeffs(Z, *(brin[i]));

            out = (*(brin[i])).MusclePump(Xs[i],Xu[i],Tcur);
            P[i] = out[0];
            Ps[i] = out[1];
        }

        for (long i = 0; i < szou; i++){

            j = i + szin;
            e[j] = - 1;
            F[j] = 0;
            //tie(alf[j],bet[j]) = OutgoingCompatibilityCoeffs(Z, *(brou[i]));

            out = (*(brou[i])).MusclePump(Xs[j],Xu[j],Tcur);
            P[j] = out[0];
            Ps[j] = out[1];
        }

        YAC.Zero();
        double vlv, vlvd;

        vlv = VLVFunc(Xs[0], alf[0], bet[0]);
        //vlvd = VLVFuncDer(Xs[0], alf[0], bet[0]);
        for (long i = 0; i < N; i++)
            YAC(0,i) = e[i]*(2*alf[i]*Xs[i]+ bet[i]);

        for (long i = 1; i < N; i++){

            /*YAC(i,i-1) = Ps[i-1] - vlv*(R*(2*alf[i-1]*Xs[i-1]+ bet[i-1])) + (P[i-1] - P[i]) - vlvd * R*Xs[i-1]*(alf[i-1]*Xs[i-1] + bet[i-1]);
            YAC(i,i) = (-1)*Ps[i];*/
            YAC(i,i-1) = vlv*Ps[i-1] -(R*(2*alf[i-1]*Xs[i-1]+ bet[i-1])) + vlv*(1-vlv) * (P[i-1] - P[i]) -  R*Xs[i-1]*(alf[i-1]*Xs[i-1] + bet[i-1]);
            YAC(i,i) = (-1)*Ps[i]*vlv;
        }

        return YAC;
    }

    /*
    Calculate velocities at the end of Newton's method iteration
    */
    vector<double> TDNewtonPostproc (Zadacha& Z, vector<double> Xs, vector<Vetv*> brin, vector<Vetv*> brou, long N, vector<double> alf,  vector<double> bet){

        vector<double> U(N);
        long szin, szou, j;

        szin = brin.size();
        szou = brou.size();

        for (long i = 0; i < szin; i++){

            //tie(alfa,beta) = IncomingCompatibilityCoeffs(Z, *(brin[i]));

            U[i] = alf[i]*Xs[i] + bet[i];
           
        }

        for (long i = 0; i < szou; i++){

            j = i + szin;
            //tie(alfa,beta) = OutgoingCompatibilityCoeffs(Z, *(brou[i]));

            U[j] = alf[j]*Xs[j] + bet[j];
        }

        return U;

    }


    void CalculateCommonKnot(Zadacha& Z, long kID, vector<Vetv*> brin, vector<Vetv*> brou, long N){

        long szin, szou, Nmax, i, j;

        double Det, X0_norm, F_norm, F2_norm, F_2_norm, eps;

        Matrix YAC(N,N), YAC_inv(N,N);

        vector<double> F(N), X0(N), Xs(N), Xu(N), X0_cur(N), F2(N), F_2(N);
        vector<double> alf(N), bet(N);

        szin = brin.size();
        szou = brou.size();
        //cout<<kID<<' '<<szin<<' '<<szou<<endl;
        Nmax = 1000;
        
        for ( long i = 0; i < szin; i++ ){

            Xs[i] = (*(brin[i])).VBO[0][(*(brin[i])).pts - 1];
            Xu[i] = (*(brin[i])).VBO[1][(*(brin[i])).pts - 1];

            //Xs[i] = (*(brin[i])).VBO[0][(*(brin[i])).pts - 1];
            //Xu[i] = (*(brin[i])).VBO[1][(*(brin[i])).pts - 1];
            //cout<<Xu[i]<<' '<<kID<<" in"<<T<<endl;
            tie(alf[i],bet[i]) = IncomingCompatibilityCoeffs(Z, *(brin[i]));
        }

        for ( long i = 0; i < szou; i++ ){

            //Xs[szin + i] = (*(brou[i])).VBO[0][1];
            //Xu[szin + i] = (*(brou[i])).VBO[1][1];
            j = i + szin;
            Xs[j] = (*(brou[i])).VBO[0][0];
            Xu[j] = (*(brou[i])).VBO[1][0];
            //cout<<Xu[i]<<' '<<kID<<" out"<<endl;
            tie(alf[j],bet[j]) = OutgoingCompatibilityCoeffs(Z, *(brou[i]));
        }

        i = 1;
        X0_norm = 1.0;
        eps = 0.00001;
        
        //if ((kID == 9)or(kID == 11) or(kID ==15)) {
        
        /*if (kID == 15){
            
            while ((i < Nmax)&&(X0_norm > eps)){
                
                F = FUZ_Ber_pump(Z, Xs, Xu, brin, brou, N, alf, bet);
                
                YAC = TDYacobian_Ber_pump(Z, Xs, Xu, brin, brou, N, alf, bet);
                
                YAC_inv = YAC.InvMatrix();
                
                for (long j = 0; j < N; j++)
                    F[j] = (-1)*F[j];
                
                X0 = YAC_inv.MulVr(F);
                
                for (long j = 0; j < N; j++)
                    Xs[j] = Xs[j] + X0[j];
                
                X0_norm = 0;
                
                for (long j = 0; j < N; j++)
                    X0_norm = X0_norm + X0[j]*X0[j];
                
                X0_norm = sqrt(X0_norm);
                
                Xu = TDNewtonPostproc(Z,Xs, brin, brou, N, alf, bet);
                
                i++;
                
            }
        }
        else {*/
        while ((i < Nmax)&&(X0_norm > eps)){
            
            F = FUZ_Ber(Z, Xs, Xu, brin, brou, N, alf, bet);
            
            YAC = TDYacobian_Ber(Z, Xs, Xu, brin, brou, N, alf, bet);
            
            YAC_inv = YAC.InvMatrix();
            
            for (long j = 0; j < N; j++)
                F[j] = (-1)*F[j];
            
            X0 = YAC_inv.MulVr(F);
            
            for (long j = 0; j < N; j++)
                Xs[j] = Xs[j] + X0[j];
            
            X0_norm = 0;
            
            for (long j = 0; j < N; j++)
                X0_norm = X0_norm + X0[j]*X0[j];
            
            X0_norm = sqrt(X0_norm);
            
            Xu = TDNewtonPostproc(Z,Xs, brin, brou, N, alf, bet);
            
            
            i++;
            
        }
        //}
            

        if (i > Nmax)
            cout << "Warning [CalculateCommonKnot]: Iteration Limit ; kID " << kID << endl;

        /*if (T > 0.001){
        //if (kID == 2){
            for ( long i = 0; i < N; i++ ){

                cout << i << "   Xs = " << Xs[i] <<"   Xu = " << Xu[i] << endl;
            }
            cout << "  T =  "  << T << endl;
        }*/


        // save results
        for ( long i = 0; i < szin; i++ ){

            (*(brin[i])).VB[0][(*(brin[i])).pts - 1] = Xs[i];
            (*(brin[i])).VB[1][(*(brin[i])).pts - 1] = Xu[i];
            //if (kID ==5){ cout <<Xs[i]<<" Sin "<<Xu[i]<< " Uin "<< kID<<" kIDin"<<' '<<(*(brin[i])).ID <<endl;
            //    cout <<(*(brin[i])).VB[1][(*(brin[i])).pts -1] <<' '<<(*(brin[i])).ID <<endl;
            //}
        }

        for ( long i = 0; i < szou; i++ ){

            (*(brou[i])).VB[0][0] = Xs[szin + i];
            (*(brou[i])).VB[1][0] = Xu[szin + i];
            //if (kID ==2) cout <<Xs[i]<<" Sout "<<Xu[i]<< " Uout "<< kID<<" kIDout"<<endl;
        }
        

    }
/*void CalculateValve(Zadacha& Z, long kID, vector<Vetv*> brin, vector<Vetv*> brou, long N){
        long szin, szou, Nmax, i;
        vector<double> F(N), Xs(N), Xu(N), F_2(N), P(N), Ps(N), out(3), X0(N);
        vector<double> alf(N), bet(N), e(N);
        Matrix YAC(N,N), YAC_inv(N,N);
        double vlv, R, X0_norm, eps, rho;
        Nmax = 1000;
        
        szin = brin.size();
        szou = brou.size();
        
        Xs[0] = (*(brin[0])).VBO[0][(*(brin[0])).pts - 1];
        Xu[0] = (*(brin[0])).VBO[1][(*(brin[0])).pts - 1];
        Xs[1] = (*(brou[0])).VBO[0][0];
        Xu[1] = (*(brou[0])).VBO[1][0];
        out = (*(brin[0])).URSOB(Xs[0],Xu[0]);
        P[0] = out[0];
        Ps[0] = out[1];
        out = (*(brou[0])).URSOB(Xs[1],Xu[1]);
        P[1] = out[0];
        Ps[1] = out[1];
        vlv = VLVFunc(Xs[0], alf[0], bet[0]);
        
        
        tie(alf[0],bet[0]) = IncomingCompatibilityCoeffs(Z, *(brin[0]));
        tie(alf[1],bet[1]) = OutgoingCompatibilityCoeffs(Z, *(brou[0]));
        rho = 1.0;
        R = 1.0;

        i = 1;
        X0_norm = 1.0;
        eps = 0.00001;
        
        
        e[0] = 1;
        e[1] = -1;
           
        F[0] = (alf[0]*Xs[0] + bet[0])*Xs[0] - (alf[1]*Xs[1] + bet[1])*Xs[1];
        
        if (alf[0]*Xs[0] + bet[0] < 0)
            F[1] = 0*(P[0] - P[1]) - Xs[0]*(alf[0]*Xs[0] + bet[0]);
        else
            F[1] = (P[0] - P[1]) - Xs[0]*(alf[0]*Xs[0] + bet[0]);
        
      for (long i = 0; i < N; i++)
            YAC(0,i) = e[i]*(2*alf[i]*Xs[i]+ bet[i]);
        if (alf[0]*Xs[0] + bet[0] < 0){
            for (long i = 1; i < N; i++){
                YAC(i,i-1) = 0;
                YAC(i,i) = 0;
            }
        }
        else {
            for (long i = 1; i < N; i++){
                YAC(i,i-1) = Ps[i-1] + (alf[i-1]*Xs[i-1] + bet[i-1])*alf[i-1];;
                YAC(i,i) = (-1)*Ps[i];
            }
        }

        
                

        //YAC = TDYacobian_Puaz(Z, Xs, Xu, brin, brou, N, alf, bet, R);

        YAC_inv = YAC.InvMatrix();

        for (long j = 0; j < N; j++)
            F[j] = (-1)*F[j];

        X0 = YAC_inv.MulVr(F);

        for (long j = 0; j < N; j++)
            Xs[j] = Xs[j] + X0[j];

        

        Xu = TDNewtonPostproc(Z,Xs, brin, brou, N, alf, bet);
        

    
        
        if (i > Nmax)
            cout << "Warning [CalculateCommonKnot]: Iteration Limit ; kID " << kID << endl;
    
        // save results
        for ( long i = 0; i < szin; i++ ){

            (*(brin[i])).VB[0][(*(brin[i])).pts - 1] = Xs[i];
            (*(brin[i])).VB[1][(*(brin[i])).pts - 1] = Xu[i];
        }

        for ( long i = 0; i < szou; i++ ){
            (*(brou[i])).VB[0][0] = Xs[szin + i];
            (*(brou[i])).VB[1][0] = Xu[szin + i];
        }

    }
        
 */
        
        
        void CalculateValve(Zadacha& Z, long kID, vector<Vetv*> brin, vector<Vetv*> brou, long N){

        long szin, szou, Nmax, i;

        double Det, X0_norm, F_norm, F2_norm, F_2_norm, eps, R;

        Matrix YAC(N,N), YAC_inv(N,N);

        vector<double> F(N), X0(N), Xs(N), Xu(N), X0_cur(N), F2(N), F_2(N);
        vector<double> alf(N), bet(N);

        szin = brin.size();
        szou = brou.size();

        //Nmax = 1000;
        Nmax = 1000;

        for ( long i = 0; i < szin; i++ ){

           // Xs[i] = (*(brin[i])).VBO[0][(*(brin[i])).pts - 2];
            //Xu[i] = (*(brin[i])).VBO[1][(*(brin[i])).pts - 2];

            Xs[i] = (*(brin[i])).VBO[0][(*(brin[i])).pts - 1];
            Xu[i] = (*(brin[i])).VBO[1][(*(brin[i])).pts - 1];
        }

        for ( long i = 0; i < szou; i++ ){

            //Xs[szin + i] = (*(brou[i])).VBO[0][1];
            //Xu[szin + i] = (*(brou[i])).VBO[1][1];

            Xs[szin + i] = (*(brou[i])).VBO[0][0];
            Xu[szin + i] = (*(brou[i])).VBO[1][0];
        }

        i = 1;
        X0_norm = 1.0;
        //eps = 0.00001;
        eps = 0.00001;
        tie(alf[0],bet[0]) = IncomingCompatibilityCoeffs(Z, *(brin[0]));
        tie(alf[1],bet[1]) = OutgoingCompatibilityCoeffs(Z, *(brou[0]));
        //cout<<alf[0]<<' '<<bet[0]<< ' '<<alf[1]<<' '<<bet[1]<<endl;
        R = 0.1;
        //R = 0.1;
        while ((i < Nmax)&&(X0_norm > eps)){
                
                F = FUZ_Puaz(Z, Xs, Xu, brin, brou, N, alf, bet, R);
                //cout<<F[0]<<F[1]<<endl;
                
                YAC = TDYacobian_Puaz(Z, Xs, Xu, brin, brou, N, alf, bet, R);
                
                YAC_inv = YAC.InvMatrix();
                
                for (long j = 0; j < N; j++)
                    F[j] = (-1)*F[j];
                
                X0 = YAC_inv.MulVr(F);
                
                for (long j = 0; j < N; j++)
                    Xs[j] = Xs[j] + X0[j];
                
                X0_norm = 0;
                
                for (long j = 0; j < N; j++)
                    X0_norm = X0_norm + X0[j]*X0[j];
                
                X0_norm = sqrt(X0_norm);
                
                Xu = TDNewtonPostproc(Z,Xs, brin, brou, N, alf, bet);
                
                
                i++;
                
                
            }
        
                
            

        if (i > Nmax)
            cout << "Warning [CalculateCommonKnot]: Iteration Limit ; kID " << kID << endl;
     
        /*if (T > 0.001){
        //if (kID == 2){
            for ( long i = 0; i < N; i++ ){

                cout << i << "   Xs = " << Xs[i] <<"   Xu = " << Xu[i] << endl;
            }
            cout << "  T =  "  << T << endl;
        }*/

    
        // save results
        for ( long i = 0; i < szin; i++ ){

            (*(brin[i])).VB[0][(*(brin[i])).pts - 1] = Xs[i];
            (*(brin[i])).VB[1][(*(brin[i])).pts - 1] = Xu[i];
        }

        for ( long i = 0; i < szou; i++ ){

            (*(brou[i])).VB[0][0] = Xs[szin + i];
            (*(brou[i])).VB[1][0] = Xu[szin + i];
        }


    }
        


    inline void Grtoch(Zadacha& Z , Derevo& Tr , Uzel& kn){
        

        if ((kn.IG == FLOW)&&((kn.Nin + kn.Nou) == 1))
            TDGrtoch(Z, Tr, kn); // inner or outer knot
        else if (((kn.Nin + kn.Nou) > 1)&&(kn.IG == BIFURCATION)){
            /*if (kn.ID == 6){
                cout<<"before CalculateCommonKnot"<<endl;
                cout<<"br.ID 1"<<endl;
                cout<<TreeLst[0].B[0].VB[1][TreeLst[0].B[0].pts-1]<<" VB"<<endl;
                cout<<endl;
            }*/
        //CallID is removed
            CalculateCommonKnot(Z, kn.ID, kn.Bin, kn.Bou, kn.Nou + kn.Nin); // simplified, kn.ID - for error messages
            /*if (kn.ID == 6){
                cout<<"after CalculateCommonKnot"<<endl;
                cout<<"br.ID 1"<<endl;
                cout<<TreeLst[0].B[0].VB[1][TreeLst[0].B[0].pts-1]<<" VB"<<endl;
                cout<<endl;
            }*/
            
            
            /* Debug

            cout << "kn.ID  " << kn.ID << endl;

            for (long i = 0; i < kn.Bin.size(); i++ ){

                    cout << "IDin  " << (*kn.Bin[i]).ID << endl;
                    cout << "S  " << (*(kn.Bin[i])).VB[0][(*(kn.Bin[i])).pts - 1];
                    cout << "   U  " << (*(kn.Bin[i])).VB[1][(*(kn.Bin[i])).pts - 1] << endl;
            }

            for (long i = 0; i < kn.Bou.size(); i++ ){

                    cout << "IDou  " << (*kn.Bou[i]).ID << endl;
                    cout << "S  " << (*(kn.Bou[i])).VB[0][0];
                    cout << "   U  " << (*(kn.Bou[i])).VB[1][0] << endl;
            }

            cout << "kn.ID  " << kn.ID << endl;
            */


        }
        if (kn.IG == VALVE ){
            CalculateValve(Z, kn.ID, kn.Bin, kn.Bou, kn.Nou + kn.Nin);
            //Debug

            /*cout << "kn.ID  " << kn.ID << endl;

            for (long i = 0; i < kn.Bin.size(); i++ ){

                    cout << "IDin  " << (*kn.Bin[i]).ID << endl;
                    cout << "S  " << (*(kn.Bin[i])).VB[0][(*(kn.Bin[i])).pts - 1];
                    cout << "   U  " << (*(kn.Bin[i])).VB[1][(*(kn.Bin[i])).pts - 1] << endl;
            }

            for (long i = 0; i < kn.Bou.size(); i++ ){

                    cout << "IDou  " << (*kn.Bou[i]).ID << endl;
                    cout << "S  " << (*(kn.Bou[i])).VB[0][0];
                    cout << "   U  " << (*(kn.Bou[i])).VB[1][0] << endl;
            }

            cout << "kn.ID  " << kn.ID << endl;
            //cout << "kn.ID VALVE  " << kn.ID << endl;
             */
        }
    
         
            
    }
};
