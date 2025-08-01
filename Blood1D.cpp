//
//  Blood1D.cpp
//  thesis
//
//  Created by  Amina on 11.03.2025.
//

#include "Blood1D.hpp"
// 1DLym.cpp: определяет точку входа для консольного приложения.

//#pragma once
//#include "TreeInit.h"
#include "CalculationMethods.h"
#include <chrono>
//using namespace TreeInit;
using namespace CalculationMethods;

int main( ) {
    long steps_skipped;
    double lastWriteT;
    double stepWriteT;
    //string protocol_str , SChar , STemp , LString;

    Zadacha Z;
    MultiKnotTreeGroup LKnotTreeGroup;
    MultiKnot MKnot;
    Vetv Branch;
    Uzel Knot;

    cout << "Start intialization: " << endl;

    LoadGlobals( );
    timeForTrack = 0.;    //для отслеживания времени работы
    dt = 0.00001;                // Начальный шаг по времени
    T = 0.;                        // Начальное время
    steps_skipped = 0;    // кол-во шагов по времени без записи
    Nz = 0;                        // кол-во записанных шагов по времени
    last_write = 0.;
    last_write_dump = 0.; // время последнего сохранения промежуточных расчетов

    lastWriteT = 0.0;
    stepWriteT = 0.5;

    //------------------------------ Инициализация задачи --------------------------

    Z.IdentifyTask( );
    //Z.PrintTask();


    //------------------------ Инициализация сосудистых сетей ----------------------
    if ( Z.Ntr > 0 ) {
        TreeLst.resize( Z.Ntr );
        for ( long i = 0; i < Z.Ntr; i++ ) {
            TreeLst [ i ].ID = i+1;
            TreeInitialization( TreeLst [ i ] , Z );
            TreeInitializationTD( TreeLst [ i ] , Z );
            InitDataTD( TreeLst [ i ] , Z );
        }
        //LoadMultiKnots( );
    }
    //Call WriteTreeParams(Z)

    //------------------------- Инициализация модели сердца ------------------------
    // !!! InitHeart( );
    //--------------------------- Цикл по времени ----------------------------------
    cout  << "Start calculations" << endl;

    //std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    //std::chrono::steady_clock::time_point curtime;
    //std::chrono::steady_clock::time_point mark1;
    //std::chrono::steady_clock::time_point mark2;
    //std::chrono::steady_clock::time_point mark3;
    //std::chrono::steady_clock::duration dur_nodes, dur_inn;

    //dur_inn = std::chrono::milliseconds(0);
    //dur_nodes = std::chrono::milliseconds(0);
//#define SOL
//#ifdef SOL
    cout << "T = " << T << " , Time_calc = " << time_calc << " , dt = " << dt << endl;

    
    while ( ( T < time_calc ) && ( dt > 0. ) ) {
        //cout<<TreeLst[0].B[0].VB[1][TreeLst[0].B[0].pts-1]<<" T1 "<<T<<endl;
        Smax = 0.; //Globals::Smax = 0.0;
        T += dt;
        if (T > last_write){
            //cout << last_write << endl;
            last_write += write_delta;

            TDWrite(Z);
        }
        //mark1 = std::chrono::steady_clock::now();
        //cout<<TreeLst[0].B[0].VB[1][TreeLst[0].B[0].pts-1]<<" T2 "<<T<<endl;

        CalcEdges(Z);
        //cout<<"CalcEdges"<<endl;
        //cout<<TreeLst[0].B[0].VBO[1][TreeLst[0].B[0].pts-1]<<" T3 "<<T<<endl;

        //mark2 = std::chrono::steady_clock::now();

        CalcVesselVertexes(Z);
        //cout<<"CalcVesselVertexes"<<endl;
        //cout<<TreeLst[0].B[0].VB[1][TreeLst[0].B[0].pts-1]<<" T4 "<<T<<endl;

        //mark3 = std::chrono::steady_clock::now();


        //CalcBronchTree( Z );
        //CalcMatterTransfer( Z );
        dt = GetDtAboutSMax( Z , Smax ); //Находится рекомендованный правилом Куранта шаг

        //------------------------------ Запись результатов ----------------------------
        if (Z.Debug == 2)
        {
            steps_skipped = steps_skipped + 1;  // Количество циклов с момента прошлой записи
            if ( WriteCalcData( Z ) == 1 )
            {
                steps_skipped = 0;
            }
        }

         //dur_inn = dur_inn + (mark2 - mark1);
         //dur_nodes = dur_nodes + (mark3 - mark2);

     /*   if (T > lastWriteT){

           // curtime = std::chrono::steady_clock::now();

            cout << "T = " << T << "  Calc. time = " << std::chrono::duration_cast<std::chrono::seconds>(curtime - start).count()  << " sec"<<endl;
           // cout << "Inner points  " <<   std::chrono::duration_cast<std::chrono::milliseconds>(mark2 - mark1).count()  << " ms"<<endl;
            // cout << "Boundaries  " <<   std::chrono::duration_cast<std::chrono::milliseconds>(mark3 - mark2).count()  << " ms"<<endl;

            if (Z.Debug > 0)
            {
                cout << "Inner points  " <<   std::chrono::duration_cast<std::chrono::seconds>(dur_inn).count()  << " sec"<<endl;
                cout << "Boundaries  " <<   std::chrono::duration_cast<std::chrono::seconds>(dur_nodes).count()  << " sec"<<endl;
            }




            lastWriteT = lastWriteT + stepWriteT;
        }*/
        //if ( TestToBreak( Z ) == 0 ) break;
        if ( ( T > last_write_dump + write_dump ) && ( Z.useSaveDump == 1 ) ) {
            for ( long i = 0; i < Z.Ntr; i++ ) {
                SaveDump( TreeLst [ i ] , Z );
            }
            last_write_dump += write_dump;
        }
        //cout << "T = " << T << " , Time_calc = " << time_calc << " , dt = " << dt << endl;;


    }  // end of main time loop
//#endif
    //----------converting results for JCNetwork-------------------
/*#ifdef JCNetwork
    for ( long i = 0; i < Z.Ntr; i++ ) {
        ConvertResults( TreeLst [ i ] , Z );
    }
#endif
 */

    //WriteFFR(TreeLst[0] , Z);

    cout << "WORKS!!!" << endl;
    //getchar( );
}

