#include "TreeInitTD.h"

namespace TreeInit {
	using namespace TreeInitTD;

	inline void TreeInitialization( Derevo& Tr , Zadacha& Z ) {
		Uzel *knt1 , *knt2;
		string NID;
		long tmp_i, addBr, addKn;
		double tmp_d;

		long ID , KnID1 , KnID2;

		NID = '1';
        Tr.dirname = "tree" + trim( Adjustl( NID ) );
		//"C:\Users\Dmitr\source\repos\1d_model\t_retinal_structure"

		/* not used in coronary
		{
			char wrotmnenogi [ 10 ];
			sprintf( wrotmnenogi , "%ld" , Tr.ID );
			NID = string( wrotmnenogi );
		}
		if ( Tr.ID == BRONCHIAL ) {
			NID = "bronch";
			//InitAlveolar();
		}
		// trim == long to string //
        Tr.dirname = "tree" + trim( Adjustl( NID ) );
		Tr.treefilename = trim( SharedDirectory ) + "tree" + trim( Adjustl( NID ) ) + slash + "tree.tre";
		Tr.knotfilename = trim( SharedDirectory ) + "tree" + trim( Adjustl( NID ) ) + slash + "knot.tre";
		Tr.branchfilename = trim( SharedDirectory ) + "tree" + trim( Adjustl( NID ) ) + slash + "branch.tre";
		Tr.TDknotExternalFilename = trim( SharedDirectory ) + "tree" + trim( Adjustl( NID ) ) + slash + "tdknotexternal.tre";
		Tr.TDbranchfilename = trim( SharedDirectory ) + "tree" + trim( Adjustl( NID ) ) + slash + "tdbranch.tre";
		Tr.TDImpactfilename = trim( SharedDirectory ) + "tree" + trim( Adjustl( NID ) ) + slash + "tdimpact.tre";*/
		Tr.treefilename = trim( SharedDirectory ) + "tree" + trim( Adjustl( NID ) ) + slash + "tree.tre";
		Tr.knotfilename = trim( SharedDirectory ) + "tree" + trim( Adjustl( NID ) ) + slash + "knot.tre";
		Tr.branchfilename = trim( SharedDirectory ) + "tree" + trim( Adjustl( NID ) ) + slash + "branch.tre";
		Tr.inputdata =  trim( SharedDirectory ) + "input_data.tre";
		Tr.FFRfile =  trim( SharedDirectory ) + "FFR.tre";

        ofstream fou;
        fou.open(Tr.FFRfile);
        fou << (-1) << endl;
        fou.close();

		cout << "Intialization of tree " << endl;
		ifstream fin( Tr.treefilename , ifstream::in );
		fin >> tmp_i;
		Tr.Nkn = tmp_i ;

		fin >> tmp_i;
		Tr.Nbr = tmp_i ;

		//cout << Tr.Nkn << endl;

		//cout << Tr.treefilename << endl;

		fin.close( );

		Tr.K = new  Uzel [ Tr.Nkn ];
		Tr.B = new  Vetv [ Tr.Nbr ];


		//----------------------------------- Reading Knots ----------------------------------

		ifstream fin1( Tr.knotfilename , ifstream::in );


		for ( long i = 0; i < Tr.Nkn; i++ ) {
			fin1 >> ID;
			Tr.K [ ID - 1 ].ID = ID;
			fin1 >> Tr.K [ ID - 1 ].C.X;
			fin1 >> Tr.K [ ID - 1 ].C.Y;	//
			fin1 >> Tr.K [ ID - 1 ].C.Z;	//
			fin1 >> Tr.K [ ID - 1 ].IG; //
            //cout << "Tr.K  " << Tr.K [ID-1].IG << endl;
		}

		fin1.close( );



		//----------------------------------- Reading Branches  --------------------------------


		//Aorta


		ifstream fin2( Tr.branchfilename , ifstream::in );


		for ( long i = 0; i < Tr.Nbr; i++ ) {
			fin2 >> ID;
			ID = ID;
			Tr.B [ ID - 1 ].ID = ID;
			Tr.B [ ID - 1 ].myTreeID = Tr.ID;
			Tr.B [ ID - 1 ].InvertPoints = 0;

			fin2 >> KnID1 >> KnID2;		//
			fin2 >> Tr.B [ ID - 1 ].len;		//
			fin2 >> Tr.B [ ID - 1 ].width;	//
			fin2 >> Tr.B [ ID - 1 ].pts;		//

			// check if enough points
			if ( ( Tr.B [ ID - 1].len > 15.0 ) && ( Tr.B [ ID - 1 ].pts < 15 ) ) {
				Tr.B [ ID - 1 ].pts = 15;
			}

			Tr.B [ ID - 1 ].Kn1 = &(Tr.K [ KnID1 - 1  ]);

			/*cout << Tr.B[ID - 1].Kn1 <<endl;
			cout << (*(Tr.B[ID - 1].Kn1)).ID <<endl;
			cout << Tr.B[ID - 1].Kn1 ->ID << endl;
			cout << Tr.B[ID - 1].Kn1->C.X << endl;
			cout <<Tr.K[0].C.X << endl;
			cout << Tr.K[0].ID << endl;*/
			Tr.B [ ID - 1 ].Kn2 = &(Tr.K [ KnID2 - 1  ]);

			if ( Tr.B [ ID - 1 ].pts > 1 ) {
				Tr.B [ ID - 1 ].dx = Tr.B [ ID - 1 ].len / double( Tr.B [ ID - 1 ].pts - 1 );
			}
			else {
				cout << "Error! Wrong points number at (tree, branch): " << Tr.ID << Tr.B [ ID - 1 ].ID << endl;
				Tr.B [ ID - 1 ].dx = Tr.B [ ID - 1 ].len;
			}





			Tr.B [ ID - 1 ].VB.resize( Z.Cor );
			Tr.B [ ID - 1 ].VBO.resize( Z.Cor );
			for ( long j = 0; j < Z.Cor; j++ ) {
				Tr.B [ ID - 1 ].VB [ j ].resize( Tr.B [ ID - 1 ].pts );
				Tr.B [ ID - 1 ].VBO [ j ].resize( Tr.B [ ID - 1 ].pts );
			}
            //cout << Tr.B [ ID - 1 ].VBO [ 1 ][Tr.B [ ID - 1 ].pts-1]<< endl;
            //cout << "Tr.B  " << Tr.B [ID-1].ID << endl;

		}

		fin2.close( );





		//----- Making list of branches incoming and outgoing to every knot -----
		for ( long i = 0; i < Tr.Nkn; i++ ) {
			Tr.K [ i ].Nou = 0;
			Tr.K [ i ].Nin = 0;
		}
		long Iin, Iou;
		if ( ( Tr.ID == PULMVEN ) || ( Tr.ID == SYSVEN ) ) {
			//  Venous tree edges direction e1 <-- e2
			for ( long i = 0; i < Tr.Nbr; i++ ) {
				Tr.B [ i ].InvertPoints = 1;
				Iin = Tr.B [ i ].Kn1 -> ID;
				Iou = Tr.B [ i ].Kn2 -> ID;
                cout << "Iin " << Iin << endl;
                cout << "Iou " << Iou << endl;
				Tr.K [ Iin - 1 ].Nin += 1;
				Tr.K [ Iou - 1 ].Nou += 1;
				Tr.K [ Iin - 1 ].Bin.push_back( &Tr.B [ i ] );
				Tr.K [ Iou - 1 ].Bou.push_back( &Tr.B [ i ] );
			}
		}
		else {
			//  Arterial and Lymphatic trees edges direction e1 --> e2
			for ( long i = 0; i < Tr.Nbr; i++ ) {
				Tr.B [ i ].InvertPoints = 0;
				Iou = Tr.B [ i ].Kn1 -> ID;
				Iin = Tr.B [ i ].Kn2 -> ID;
				Tr.K [ Iin - 1 ].Nin += 1;
				Tr.K [ Iou - 1 ].Nou += 1;
				Tr.K [ Iin - 1 ].Bin.push_back( &Tr.B [ i ] );
				Tr.K [ Iou - 1 ].Bou.push_back( &Tr.B [ i ] );
			}
		}

		// Marking input and output as BIFURCATION
		for ( long i = 0; i < Tr.Nkn; i++ ) {
			if ( ( Tr.K [ i ].Nin + Tr.K [ i ].Nou ) == 1 ) {
				Tr.K [ i ].IG = FLOW;
            } else if(( Tr.K [ i ].Nin==1)&&( Tr.K [ i ].Nou==1)){
                Tr.K [ i ].IG = VALVE;
            } else {
				Tr.K [ i ].IG = BIFURCATION;
			}
		}

        if (Z.Debug == 2)
        {
            cout << "Total knots loaded: " << Tr.Nkn << endl << "--------------------------------------" << endl;
            for ( long i = 0; i < Tr.Nkn; i++ ) printUzel ( Tr.K[i] );
            cout << "Total branches loaded: " << Tr.Nbr << endl << "--------------------------------------" << endl;
            for ( long i = 0; i < Tr.Nbr; i++ ) printVetv ( Tr.B[i] );
        }

	}

	inline void LoadMultiKnots( ) {  // not used
		long  KntID , szin , szou , LTrID;
		Derevo LTree;

		Globals::MKnots.multiKnotsfilename = trim( SharedDirectory ) + "multiknots.tre";
		ifstream fin( Globals::MKnots.multiKnotsfilename , ifstream::in );
		fin >> Globals::MKnots.sz;

		if ( Globals::MKnots.sz > 0 ) {
			Globals::MKnots.Lst.resize( Globals::MKnots.sz );
			for ( long i = 0; i < Globals::MKnots.sz; i++ ) {
				fin >> Globals::MKnots.Lst [ i ].sz;
				Globals::MKnots.Lst [ i ].GrpLst.resize( Globals::MKnots.Lst [ i ].sz );
				szin = 0;
				szou = 0;
				for ( long j = 0; j < Globals::MKnots.Lst [ i ].sz; j++ ) {
					fin >> MKnots.Lst [ i ].GrpLst [ j ].TrID;
					fin >> MKnots.Lst [ i ].GrpLst [ j ].sz;
					Globals::MKnots.Lst [ i ].GrpLst [ j ].KntLst.resize( Globals::MKnots.Lst [ i ].GrpLst [ j ].sz );
					for ( long k = 0; k < Globals::MKnots.Lst [ i ].GrpLst [ j ].sz; k++ ) {
						fin >> KntID;
						LTrID = MKnots.Lst [ i ].GrpLst [ j ].TrID;
						LTree = TreeLst [ LTrID ];
						Globals::MKnots.Lst [ i ].GrpLst [ j ].KntLst [ k ] = LTree.K [ KntID ];
						szin = szin + Globals::MKnots.Lst [ i ].GrpLst [ j ].KntLst [ k ].Nin;
						szou = szou + Globals::MKnots.Lst [ i ].GrpLst [ j ].KntLst [ k ].Nou;
					}
				}
				Globals::MKnots.Lst [ i ].szin = szin;
				Globals::MKnots.Lst [ i ].szou = szou;
			}
		}
		fin.close( );

		cout << "Total multiknots loaded: " << MKnots.sz << endl;
	}
};
