//########################################################################
//  Модуль определения задача-зависимых параметров 
//  Задача-зависимыми называются параметры, наличие или отсутствие 
//  которых определяется типом задачи
//########################################################################
#include <string>
#include <vector>
#include "CalcUtils.h"

namespace TreeDataDepend {
	using namespace std;
	using namespace CalcUtils;

	const static long IT_BLOOD = 1;
	const static long IT_SUBSTANCE = 2;
	const static long GT_LUMPED = 1;
	const static long GT_DISTRIBUTED = 2;

	//________________________________________________________________________
	//________________________________________________________________________
	//     Задача-зависимые параметры ветви
	//________________________________________________________________________
	class TDVetvParameter {
	public:
		double R , Rinit;   // сопротивление (сконцентрировано в нижней части сосуда)
		double delta;    // коэффициент кровопотери в сосуде
		double c , cnewQ , cnewP;      // скорость пульсовой волны в линейном приближении
		double Ps;      // коэф. для расчета функции состояния
		double Smin , Smax;  // максимальное и минимальное значения поперечного сечения сосуда
		double S0 , S0new; // поперечное сечение в покое
		double Par , Pgr;   // Альвеолярное (Par) и плевральное (Pgr) давление, в большом кргуе должно быть 0
		double d_start;
		long gravity , pmaxflag; // гравитация
		double Qave_curr , Rave_curr , Pave_curr , RadiusAve_curr;
		double Qave_prev , Rave_prev , Pave_prev , RadiusAve_prev;
		double Qave0 , Qave_next , Pave0 , Pave_next;
		double Save_curr , Save_next , Save_prev , Save0;
		double *Pmax , *Pmaxcurr;
		double MDP , PD_aor_0 , PD_aor_EECP , MSP , PS_aor_0 , PS_aor_EECP , EE , V_0 , V_comp; //// EECP
		long EEflag;
		vector< vector<double> > conc , conc_old; // концентрации веществ
		vector<string> fnameC; // концентрации веществ записывается в эти файлы
		vector<string> fnameVar; //string fnameS , fnameU , fnameP;
		vector<string> fnameTD; //*fnamePmax , *fnamePave , *fnameQave , *fnameMDP , *fnameMSP , *fnameEE;
		string strID;
	};

	//________________________________________________________________________
	//     Компонента альвеолярного объема
	//________________________________________________________________________
	class TDAlveolar {
	public:
		double A , I , Ig , R , Rg , C , Y [ 4 ];
	};

	//________________________________________________________________________
	//     Задача-зависимые параметры узла
	//________________________________________________________________________
	class TDUzelParameter {
	public:
		long TextID;         // номер дерева содержащего внешние узлы
		long Kext_sz;         // количество внешних узлов
		vector<long> KextID; // номера внешних узлов (по версии внешнего дерева)
		TDAlveolar *alv;
		long NgasExch;
		vector<long> gasExchKnt;
		// LYMPHATIC 
		double LymJVdst;
	};

	//________________________________________________________________________
	//     Внешнее воздействие на ветвь
	//________________________________________________________________________
	class TDExternalImpact {
	public:
		long BrID;     // номер сосуда
		long IType;     // тип воздействия (сток / приток крови / концентрации веществ)
		long GType;     // тип геометрии воздействия: 1 - точечное, 2 - распределено по всей длине сосуда
		double len;       // если тип геометрии точечный, то место воздействия (в % длины сосуда, от его начала)
		long SubstID;    // номер вещества (если сток / приток концентрации)
		double alfa;       // коэффициент определяющий воздействие (скорость стока / притока и т.п.)
		double Tstart , Tend;   // моменты начала и окончания воздействия (абсолютное)
	};
}
