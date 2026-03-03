

namespace TDGranuslov {

    using namespace CalcUtils;
    using namespace TaskData;
    using namespace Subroutines;

   /* double getInVelocity(Zadacha& Z) {
        const static long GranuslovPars = 6;
        double U0, tp, ts, T_, Vs, Ves, Vd, Tc, Tcur;
        double a1, b1, c1, a2, b2, c2;

        ifstream fin("paths.ini", ifstream::in);
        fin >> Globals::SharedDirectory;
        fin.close();

        string filename = SharedDirectory + "ini" + slash + "granuslov.ini";

        double tmp[GranuslovPars];
        ReadFile(filename, tmp, GranuslovPars);
        tp = tmp[0];
        ts = tmp[1];
        T_ = tmp[2];
        Vs = tmp[3];
        Ves = tmp[4];
        Vd = tmp[5];


        if (T > (Z.T_last_Hbeat + Z.h_period_curr)) {

            // new heart cycle

            Z.T_last_Hbeat = Z.T_last_Hbeat + Z.h_period_curr;
            Z.N_heart_cycles += 1;
        }

        double delV = Vs - Vd;
        double delV_D = Ves - Vd;

        a1 = (delV_D / ts - delV / tp) / (tp - ts);
        b1 = delV / tp + a1 * tp;
        c1 = Vd;

        a2 = ((2 * a1 * ts + b1) * (delV_D * ts + Vd * T_ - Ves * Vs) - Ves * delV_D) / ((T_ - ts) * (2 * a1 * ts + b1) - delV_D);
        c2 = (Vd * T_ - Ves * ts - a2 * (T_ - ts)) / delV_D;
        b2 = (Ves - a2) * (ts + c2);

        Tcur = T - Z.T_last_Hbeat;


        cout << "Tcur = " << Tcur << endl;
        if (Tcur < ts) {
            U0 = -a1 * (Tcur * Tcur) + b1 * Tcur + c1;
        }
        else {
            U0 = a2 + b2 / (Tcur + c2);
        }

        return U0;
    } */
    double getInFlow (Zadacha& Z){

        double Qin, Tc, Tcur, Ts, Tf, freq, Tstep;

//        if (T > (Z.T_last_Hbeat + Z.h_period_curr)){
//
//            // new heart cycle
//
//            Z.T_last_Hbeat = Z.T_last_Hbeat + Z.h_period_curr;
//            Z.N_heart_cycles += 1;
//        }

        Tc = 60/Z.HR;
        Tcur = T ;
        Ts = 0;
        //Ts = 0.2*Tc;
        Tf = 0.3;
        Tstep = 1.33; //sec
        freq = 1/Tstep*60; //min

        //if ((Tcur > Ts)&&(Tcur < Tf))
         //   Qin = (Z.SV*PI/(2*(Tf - Ts)))*sin(PI*(Tcur - Ts)/(Tf - Ts));
        //else
        //    Qin = 0;
        Qin = 0.3*exp(0.0328*freq)* (1 - exp(-3*T));
        //Qin = 20 * sin(2 * PI * T / 10);
        //Qin = 0.3*exp(0.0328*freq);

        return Qin;
    }

    vector<double> TDFlowToSU (Zadacha& Z , Derevo& Tr , Uzel& kn, double Qin){

        vector<double> V(Z.Cor);
        double alfa, beta;

        if (kn.Nou == 1){ // beginning of a branch

            tie(alfa, beta) = OutgoingCompatibilityCoeffs(Z,*(kn.Bou[0]),0);

            V[1] = (beta + sqrt(beta*beta + 4*alfa*Qin))/2;
            V[0] = (- beta + sqrt(beta*beta + 4*alfa*Qin))/(2*alfa);
        }

        if (kn.Nin == 1){

            tie(alfa, beta) = IncomingCompatibilityCoeffs(Z,*(kn.Bin[0]),9);

            V[1] = (beta - sqrt(beta*beta + 4*alfa*Qin))/2;
            V[0] = (- beta - sqrt(beta*beta + 4*alfa*Qin))/(2*alfa);

        }

        return V;
    }



    // Calculate function and derivative for Newton's method

    tuple<double, double> FuncNewtS(double Scur, Vetv& br, double alfa, double beta, double Pout){

        double Ps, Pcur, Func, FuncS;
        vector<double> out(3);

        out = br.URSOB(Scur, br.VBO[1][br.pts - 1], br.ID, T, br.pts - 1);
        Pcur = out[0];
        Ps = out[1];

        if (Pcur < 0)
            cout << "Warning [FuncNewtS]: P < 0 ; brID " << br.ID << endl;

        Func = alfa*Scur*Scur + beta*Scur - (Pcur - Pout)/br.TD.R;
        FuncS = 2*alfa*Scur + beta - Ps/br.TD.R;

        return {Func, FuncS};

    }

    // Newton's method for outlet Q = (P - Pveins)/R
    vector<double> CalcRtoSUnewt (Zadacha& Z , Derevo& Tr , Vetv& br){

        vector<double> V(Z.Cor);
        double Pout, S_norm, Scur, Snew, alfa, beta, Func, FuncS, Pcor, eps;
        long i, Nmax;

        i = 0;
        Nmax = 10000;
        eps = 0.001;

        Pout = Z.Pveins*1333.2;

        Pcor = 0; // external Pressure, we assume it to be 0

        tie(alfa,beta) = IncomingCompatibilityCoeffs(Z,br,br.pts-1);

        Scur = ( -beta - sqrt(beta * beta ) ) / (2. * alfa);
        Snew = br.TD.S0*(1 + log(1 + (Pout - Pcor)/(br.TD.c*br.TD.c)));

        /*cout << "BrID: " << br.ID << endl;
        cout << "Scur: " << Scur << "   Snew:  "<< Snew <<  endl;
        cout << "alfa: " << alfa << "   beta:  "<< beta <<  endl;
        cout << endl;*/


        if ((Scur <= 0)||(Snew <=0))
            cout << "Warning [CalcRtoSUnewt]: S <= 0 ; brID " << br.ID << endl;

        if (Scur < Snew)
            Scur = Snew;

        S_norm = 1;

        while ((i < Nmax)&&(S_norm > eps)){

            tie(Func,FuncS) = FuncNewtS(Scur,br,alfa,beta,Pout);

            Snew = Scur - (Func/FuncS);

            if (abs(Scur) > 0.00000001)
                S_norm = abs((Snew - Scur)/Scur);
            else{
                cout << "Warning [CalcRtoSUnewt]: Scur = 0 ; brID " << br.ID << endl;
                S_norm = 0;
            }

            i++;
            Scur = Snew;
        }

        if (i > Nmax)
            cout << "Warning [CalcRtoSUnewt]: Iteration Limit ; brID " << br.ID << endl;

        V[0] = Snew;
        V[1] = alfa*Snew + beta;

        return V;
    }
    double VLVFunc (double S,double alf,double bet) {
        
        double vlvwidth,vlv;
        vlvwidth = 100;
        //vlv = atan(vlvwidth*(alf*S+bet));
        //vlv = atan(vlvwidth*(alf*S+bet))*2;
        vlv = 1/(1+exp(-vlvwidth*(alf*S+bet-0.028)));
        //cout<< vlv<< ' '<<(alf*S+bet)<<endl;
        return vlv;
    }
    double VLVFuncDer (double S,double alf,double bet) {
        
        double vlvwidth,vlv,v;
        //vlvwidth = 100;
        //vlv = vlvwidth/(1+pow(vlvwidth*(alf*S+bet),2));
        //vlv = vlvwidth/(1+pow(vlvwidth*(alf*S+bet),2)/PI);
        v = VLVFunc(S,alf,bet);
        vlv = v*(1-v);
        return vlv;
    }


    vector<double> CalcEndFlow (Zadacha& Z , Derevo& Tr , Vetv& br){
        
        vector<double> V(Z.Cor);
        double Pout,Scur,  alfa, beta, Pcor, Tstep, Tc, Tcur,Pmax,a,b; //Tstep - period of human step
        
        Pout = Z.Pveins*1333.2;
        Pcor = 0;

        //Pcor = Pmax/2*(1+sin(2*PI*Tcur/Tstep)); // external Pressure, we assume it to be 0 - Amina added
        
        tie(alfa,beta) = IncomingCompatibilityCoeffs(Z, br, br.pts - 1);
        
        //Scur = br.TD.S0*(1 + log(1 + (Pout - Pcor)/(br.TD.c*br.TD.c)));
        //V[0] = Scur;
        //V[1] = alfa*Scur + beta;
        //V[0] = br.VB[0][8];
        //V[1] = br.VB[1][8];
        V[0]=br.VB[0][br.pts - 1];
        V[1]=br.VB[1][br.pts - 1];
        
//        if (beta*beta + 4*alfa*br.VB[0][br.pts - 1]*br.VB[1][br.pts - 1] >= 0){
//            V[0] = abs(-beta + sqrt( beta*beta + 4*alfa*br.VB[0][br.pts - 1]*br.VB[1][br.pts - 1]));
//            V[1] = alfa*V[0]+beta;
//        } else {
//            V[0] = 1e-6;
//            V[1] = alfa*V[0]+beta;
//            
//        }
//        a = alfa;
//        b = beta;
//        double Q_prev = br.VB[0][br.pts - 2] * br.VB[1][br.pts - 2];
//
//            //a*S^2 + b*S - Q_prev = 0
//        double D = b * b + 4 * a * Q_prev; 
//        double S = 0;
//        if (D < 0) {
//            //S = 0;
//            S = br.TD.S0;
//        } else {
//            double sqrtD = sqrt(D);
//            double S1 = (-b + sqrtD) / (2 * a);
//            double S2 = (-b - sqrtD) / (2 * a);
//   
//            if (S1 > 0 && S2 > 0) {
//                S = (S1 < S2) ? S1 : S2;
//                //S = S2;
//            } else if (S1 > 0) {
//                S = S1;
//            } else if (S2 > 0) {
//                S = S2;
//            } else {
//                
//                S = br.TD.S0;
//                }
//            }
//        V[0] = S;
//        V[1] = a * S + b;
        
         
        return V;
    }

    //calculate boundary knots
    void TDGrtoch(Zadacha& Z , Derevo& Tr , Uzel& kn){

        long idx1,idx2;
        double Qin, alf, bet, Tc, Qin_coeff, Pout, Rcor, Tcur, Ts, Tf, Pconst;
        vector<double> V(Z.Cor);

        Vetv *brp = &(Tr.B[0]);

        // arterial network

        if ((Tr.ID == PULMART)||(Tr.ID == SYSART)){

            if ((kn.IG == FLOW)&&(kn.Nou == 1)){
                /*if (kn.ID == 4){
                    brp = kn.Bou[0];
                    cout<<brp<<endl;
                    idx1 = 0;
                    idx2 = 1;
                    Qin = 0;
                    V[1] = 0;
                    V[0] = (*brp).TD.S0;
                    
                }else{*/
                    //aorta
                brp = kn.Bou[0];
                
                idx1 = 0;
                idx2 = 1;
                
                Tc = 60/Z.HR;
                Qin = getInFlow(Z);
                
                //V[1] = getInVelocity(Z)*0.01;
                //tie(alf, bet) = OutgoingCompatibilityCoeffs(Z,*brp);
                //V[0] = (V[1] - bet)/alf;
                
                V = TDFlowToSU(Z, Tr, kn, Qin);
                
            }
            if ((kn.IG == FLOW)&&(kn.Nin == 1)){
                
                brp = kn.Bin[0];

                idx1 = (*brp).pts - 1;
                idx2 = (*brp).pts - 2;
                //V = CalcRtoSUnewt(Z, Tr, *brp);

                V = CalcEndFlow(Z, Tr, *brp);
                

                //terminal artery

                /*brp = kn.Bin[0];

                idx1 = (*brp).pts - 1;
                idx2 = (*brp).pts - 2;

                Tc = 60/Z.HR;
                Tcur = (T - Z.T_last_Hbeat)/Tc;
                Ts = 0;
                Tf = 0.3;

                Rcor = 0;

                V = CalcRtoSUnewt(Z, Tr, *brp);*/

            }
            /*else {
                brp = kn.Bin[0];
            }*/
            (*brp).VB[0][idx1] = V[0];
            (*brp).VB[1][idx1] = V[1];

            cout << "BrID: " << (*brp).ID << endl;
            cout << "Brlen: " << (*brp).len << "   BrD:  "<< (*brp).width <<  endl;
            cout << "V: " << V[0] << "  " << V[1] << endl;
            cout << endl;

        }

    }
};
