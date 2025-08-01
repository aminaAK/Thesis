

namespace Subroutines {

    using namespace CalcUtils;
    using namespace TaskData;

  /*
	Функция берёт значения в ближних и в дальних от узла точках
        для входящей ветви
*/
  inline tuple<vector<double>,vector<double>,vector<double>,vector<double>> OCNFI(Zadacha& Z, Vetv& br){

      vector<double> OC(Z.Cor),OF(Z.Cor),NC(Z.Cor),NF(Z.Cor);

      OC[0] = br.VBO[0][br.pts - 1];
      OC[1] = br.VBO[1][br.pts - 1];

      OF[0] = br.VBO[0][br.pts - 2];
      OF[1] = br.VBO[1][br.pts - 2];

      NC[0] = br.VB[0][br.pts - 1];
      NC[1] = br.VB[1][br.pts - 1];

      NF[0] = br.VB[0][br.pts - 2];
      NF[1] = br.VB[1][br.pts - 2];

      return {OC,OF,NC,NF};
  }

  /*
	Функция берёт значения в ближних и в дальних от узла точках
        для исходящей ветви
*/
  inline tuple<vector<double>,vector<double>,vector<double>,vector<double>> OCNFO(Zadacha& Z, Vetv& br){

      vector<double> OC(Z.Cor),OF(Z.Cor),NC(Z.Cor),NF(Z.Cor);

      OC[0] = br.VBO[0][0];
      OC[1] = br.VBO[1][0];

      OF[0] = br.VBO[0][1];
      OF[1] = br.VBO[1][1];

      NC[0] = br.VB[0][0];
      NC[1] = br.VB[1][0];

      NF[0] = br.VB[0][1];
      NF[1] = br.VB[1][1];

      return {OC,OF,NC,NF};
  }




/* Sobstvennye vektora i matritsy OM, OMO
 Stroki matritsy OM - levye sobstvennye vektora
 Stolbtsy matritsy OMO - pravyye sobstvennye vektora
 Ne pereputayte gde levye, a gde pravyye*/

    inline tuple<Matrix,Matrix> OMEGAB(vector<double> V, Zadacha& Z, double zvuk, double P, Vetv& br){

            Matrix OMO(Z.Cor,Z.Cor), OM(Z.Cor,Z.Cor);

            OMO(0,0) = V[0]/(2*zvuk);
            OMO(0,1) = V[0]/(2*zvuk);
            OMO(1,0) = -0.5;
            OMO(1,1) = 0.5;

            OM(0,0) = zvuk/V[0];
            OM(0,1) = -1;
            OM(1,0) = zvuk/V[0];
            OM(1,1) = 1;

            return {OMO,OM};

	}

	inline tuple<vector<double>,Matrix,Matrix> SZH_Hybr (vector<double> V, Zadacha& Z, double zvuk, double P, double gamma, Vetv& br){

            Matrix B_h(Z.Cor,Z.Cor), D_h(Z.Cor,Z.Cor);
            vector<double> S(Z.Cor),Sigm(Z.Cor);

            S[0] = V[1] - zvuk;
            S[1] = V[1] + zvuk;

            Sigm[0] = S[0]*dt/br.dx;
            Sigm[1] = S[1]*dt/br.dx;

            B_h(0,0) = abs(Sigm[0])*(1+5*(1 - gamma)*(1 - abs(Sigm[0]))/19)/2;
            B_h(0,1) = 0;
            B_h(1,1) = abs(Sigm[1])*(1+5*(1 - gamma)*(1 - abs(Sigm[1]))/19)/2;
            B_h(1,0) = 0;

            D_h(0,0) = 6*(1-gamma)*Sigm[0]*(1/abs(Sigm[0]) - 1)/19;
            D_h(0,1) = 0;
            D_h(1,1) = 6*(1-gamma)*Sigm[1]*(1/abs(Sigm[1]) - 1)/19;
            D_h(1,0) = 0;

            return {S,B_h,D_h};

	}


    /* Compatibility conditions coefficients U = a*S + b for the right side of the branch (incoming)  */

    tuple<double, double> IncomingCompatibilityCoeffs(Zadacha& Z, Vetv& br)
    {

        double P, Ps, zvuk, sigm, a, b, gamma, Tcur;

        vector<double> out(3);

        vector<double> S(Z.Cor), OC(Z.Cor),OF(Z.Cor),NC(Z.Cor),NF(Z.Cor), FPR(Z.Cor), V(Z.Cor);

        Matrix OM(Z.Cor,Z.Cor), OMO(Z.Cor,Z.Cor), B_h(Z.Cor,Z.Cor), D_h(Z.Cor,Z.Cor);

        tie(OC,OF,NC,NF) =  OCNFI(Z, br);
        Tcur = T - Z.T_last_Hbeat;
        //cout<< Tcur << endl;
        V = OC;
        //if ((br.ID == 16) or (br.ID == 17) or (br.ID == 18) or (br.ID== 23)){
        if ((br.ID == 11) or (br.ID == 12)){
            out = br.MusclePump(V[0],V[1],Tcur);
            P = out[0];
            Ps = out[1];
            zvuk = out[2];
            //cout << br.ID<<"incoming" <<endl;
        }
        else {
            out = br.URSOB(V[0],V[1]);
            P = out[0];
            Ps = out[1];
            zvuk = out[2];
        }

        gamma = 1; // not used here, we only need S

        tie(S, B_h, D_h) = SZH_Hybr(V, Z, zvuk, P, gamma, br);
        tie(OMO,OM) = OMEGAB(V, Z, zvuk, P, br);
        FPR = br.FPRCH(V[0],V[1],br.pts - 1);

        sigm = S[1]*dt/br.dx;
        //cout<<"sigm="<<sigm<<endl;

        a = (-1)*OM(1,0)/OM(1,1);

        if (a > 0)
        {
         cout << "Error, a > 0 for incoming BR.ID " << br.ID<< endl;
        }

        b = (OC[1] + sigm*NF[1] - a*(OC[0] + sigm*NF[0]) + FPR[1]*dt - a*FPR[0]*dt)/(1 + sigm);

        return {a,b};

    }

    /* Compatibility conditions coefficients U = a*S + b for the left side of the branch (outgoing)  */

    tuple<double, double> OutgoingCompatibilityCoeffs (Zadacha& Z, Vetv& br)
    {

        double P, Ps, zvuk, sigm, a, b, gamma, Tcur;

        vector<double> out(3);

        vector<double> S(Z.Cor), OC(Z.Cor),OF(Z.Cor),NC(Z.Cor),NF(Z.Cor), FPR(Z.Cor), V(Z.Cor);

        Matrix OM(Z.Cor,Z.Cor), OMO(Z.Cor,Z.Cor), B_h(Z.Cor,Z.Cor), D_h(Z.Cor,Z.Cor);

        tie(OC,OF,NC,NF) =  OCNFO(Z, br);
        Tcur = T - Z.T_last_Hbeat;
        //cout<< Tcur << endl;
        V = OC;
        //if ((br.ID == 16) or (br.ID == 17) or (br.ID == 18) or (br.ID== 23)){
        if ((br.ID == 11) or (br.ID == 12)){
            out = br.MusclePump(V[0],V[1], Tcur);
            P = out[0];
            Ps = out[1];
            zvuk = out[2];
            //cout << br.ID<<"outgoing" <<endl;
        }
        else {
            out = br.URSOB(V[0],V[1]);
            P = out[0];
            Ps = out[1];
            zvuk = out[2];
            
        }


        gamma = 1; // not used here, we only need S

        tie(S, B_h, D_h) = SZH_Hybr(V, Z, zvuk, P, gamma, br);

        //cout << "S: " << S[0] << "  " << S[1] << endl;
        //cout << "zvuk: " << zvuk << endl;
        tie(OMO,OM) = OMEGAB(V, Z, zvuk, P, br);
        FPR = br.FPRCH(V[0],V[1],0);

        sigm = S[0]*dt/br.dx;

        a = (-1)*OM(0,0)/OM(0,1);

         if (a < 0)
        {
            cout << "Error, a < 0 for outgoing BR.ID " << br.ID<< endl;
            
        }

        b = (OC[1] - sigm*NF[1] + a*(sigm*NF[0] - OC[0]) + FPR[1]*dt - a*FPR[0]*dt)/(1 - sigm);

        return {a,b};

    }



};
