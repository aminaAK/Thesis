

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
    double getInFlow (Zadacha& Z){//поток на входе

            double Qin, Tc, Tcur, Ts, Tf,freq,Tstep;

            if (T > (Z.T_last_Hbeat + Z.h_period_curr)){

                // new heart cycle

                Z.T_last_Hbeat = Z.T_last_Hbeat + Z.h_period_curr;
                Z.N_heart_cycles += 1;
            }

            Tc = 60/Z.HR;
            Tcur = T ;
            Ts = 0;
            //Ts = 0.2*Tc;
            Tf = 0.3;
            Tstep = 1.33; //sec
            freq = 1/Tstep*60; //min

            
            /*if ((Tcur > Ts)&&(Tcur < Tf))
                Qin = (Z.SV*PI/(2*(Tf - Ts)))*sin(PI*(Tcur - Ts)/(Tf - Ts));
            else
                Qin = 0;
             */
           //Qin = 20*(sin(2*PI*Tcur-PI*0.3)+0.8);
           Qin = 0.3*exp(0.0328*freq);
        /*if ((Tcur > 2)&&(Tcur < 4))
                Qin = 5*(sin(2*PI*Tcur-PI*0.3)+0.8);
            else
                Qin = 0;*/
            
            //Qin = Z.SV/Tc*(sin(2*PI*Tcur/Tc)+0.5);
            //cout<<"Qin "<<Qin<<endl;
            return Qin;
    }

    double VLVFunc (double S,double alf,double bet) {
        
        double vlvwidth,vlv;
        vlvwidth = 590;
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


    vector<double> TDFlowToSU (Zadacha& Z , Derevo& Tr , Uzel& kn, double Qin){

        vector<double> V(Z.Cor);
        double alfa, beta;

        if (kn.Nou == 1){ // beginning of a branch

            tie(alfa, beta) = OutgoingCompatibilityCoeffs(Z,*(kn.Bou[0]));
            //cout<<"ain bin"<<alfa<<' '<<beta<<endl;
            V[1] = (beta + sqrt(beta*beta + 4*alfa*Qin))/2;
            V[0] = (- beta + sqrt(beta*beta + 4*alfa*Qin))/(2*alfa);
            //cout<<"V "<<V[0]<<' '<<V[1]<<endl;
        }

        if (kn.Nin == 1){

            tie(alfa, beta) = IncomingCompatibilityCoeffs(Z,*(kn.Bin[0]));

            V[1] = (beta - sqrt(beta*beta + 4*alfa*Qin))/2;
            V[0] = (- beta - sqrt(beta*beta + 4*alfa*Qin))/(2*alfa);

        }

        return V;
    }



    // Calculate function and derivative for Newton's method

    tuple<double, double> FuncNewtS(double Scur, Vetv& br, double alfa, double beta, double Pout){

        double Ps, Pcur, Func, FuncS;
        vector<double> out(3);

        out = br.URSOB(Scur, br.VBO[1][br.pts - 1]);
        Pcur = out[0];
        Ps = out[1];

        if (Pcur < 0)
            cout << "Warning [FuncNewtS]: P < 0 ; brID " << br.ID << endl;

        Func = alfa*Scur*Scur + beta*Scur - (Pcur - Pout)/br.TD.R;
        FuncS = 2*alfa*Scur + beta - Ps/br.TD.R;

        return {Func, FuncS};

    }
    vector<double> CalcEndFlow (Zadacha& Z , Derevo& Tr , Vetv& br){
        
        vector<double> V(Z.Cor);
        double Pout,Scur,  alfa, beta, Pcor, Tstep, Tc, Tcur,Pmax; //Tstep - period of human step
        
        Pout = Z.Pveins*1333.2;
        Pcor = 0;
    
        //Pcor = Pmax/2*(1+sin(2*PI*Tcur/Tstep)); // external Pressure, we assume it to be 0 - Amina added
        
        tie(alfa,beta) = IncomingCompatibilityCoeffs(Z,br);
        
        Scur = br.TD.S0*(1 + log(1 + (Pout - Pcor)/(br.TD.c*br.TD.c)));
        V[0] = Scur;
        V[1] = alfa*Scur + beta;
        return V;
    }
    // Newton's method for outlet Q = (P - Pveins)/R

    vector<double> CalcRtoSUnewt (Zadacha& Z , Derevo& Tr , Vetv& br){

        vector<double> V(Z.Cor);
        double Pout, S_norm, Scur, Snew, alfa, beta, Func, FuncS, Pcor, eps, Tstep, Tc, Tcur,Pmax; //Tstep - period of human step
        long i, Nmax;
        Tstep = 0.5;
        i = 0;
        Nmax = 10000;
        eps = 0.001;
        Tc = 60/Z.HR;
        Tcur = (T - Z.T_last_Hbeat)/Tc;
        
        //Pout = Z.Pveins*1333.2;
        Pout = Z.Pveins*1333.2*0.9899;
        Pcor = 0;
        Pmax = 0.5*1333.2; // Pmax = 10 kPa 1mmHg = 133.32 Pa
        //Pcor = Pmax/2*(1+sin(2*PI*Tcur/Tstep)); // external Pressure, we assume it to be 0 - Amina added

        tie(alfa,beta) = IncomingCompatibilityCoeffs(Z,br);

        Scur = ( -beta - sqrt(beta * beta ) ) / (2. * alfa);
        Snew = br.TD.S0*(1 + log(1 + (Pout - Pcor)/(br.TD.c*br.TD.c)));
        //cout << "Scur= " << Scur <<"b="<<beta<<"a="<<alfa<<"br="<<br.ID<< endl;

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

    //calculate boundary knots
    void TDGrtoch(Zadacha& Z , Derevo& Tr , Uzel& kn){

        long idx1,idx2;
        double Qin, alf, bet, Tc, Qin_coeff, Pout, Rcor, Tcur, Ts, Tf, Pconst;
        vector<double> V(Z.Cor);

        Vetv *brp = &(Tr.B[0]);

        // arterial network

        if ((Tr.ID == PULMART)||(Tr.ID == SYSART)){

            if ((kn.IG == FLOW)&&(kn.Nou == 1)&&(kn.ID != 16)){

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
            //if ((kn.IG == FLOW)&&(kn.Nin == 1)&&(kn.ID==16)){
            if (kn.ID==16){

                brp = kn.Bou[0];

                idx1 = 0;
                idx2 = 1;

                Tc = 60/Z.HR;
                Qin = 0;

                //V[1] = getInVelocity(Z)*0.01;
                //tie(alf, bet) = OutgoingCompatibilityCoeffs(Z,*brp);
                //V[0] = (V[1] - bet)/alf;

                V = TDFlowToSU(Z, Tr, kn, Qin);
                //cout<<V[0]<<' '<<V[1]<< " S U"<< endl;
                
            }
            if ((kn.IG == FLOW)&&(kn.Nin == 1)&&(kn.ID != 16)){

                //terminal artery

                brp = kn.Bin[0];

                idx1 = (*brp).pts - 1;
                idx2 = (*brp).pts - 2;

                Tc = 60/Z.HR;
                Tcur = (T - Z.T_last_Hbeat)/Tc;
                Ts = 0;
                Tf = 0.3;

                Rcor = 0;
                V = CalcEndFlow(Z, Tr, *brp);
                

                //V = CalcRtoSUnewt(Z, Tr, *brp);

            }
            /*else {
                brp = kn.Bin[0];
            }*/
            
            //cout<<idx1<<endl;
            (*brp).VB[0][idx1] = V[0];
            (*brp).VB[1][idx1] = V[1];

           /* cout << "BrID: " << (*brp).ID << endl;
            cout << "Brlen: " << (*brp).len << "   BrD:  "<< (*brp).width <<  endl;
            cout << "V: " << V[0] << "  " << V[1] << endl;
            cout << endl;*/

        }

    }
};
