#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define Debug_Mode 0
#define MAX(a,b) ( (a) > (b) ? (a) : (b) )
#define MIN(a,b) ( (a) < (b) ? (a) : (b) )
#define ROUND(a) ( (a) > (0) ? (int)(a+0.5) : (int)(a-0.5) )
#define ABS(a) ( (a) > (0) ? (a) : -(a) )

#define Fs 50
#define Bin_size 1 // (1) 0.5 sec; (2) 1 sec
#define FilterType 1 // (0) no filter(minus the average) (1) IIR
#define NoDownSample 1 // (1) turn off down sample autocorr

#if (Bin_size == 1)
#define Di 25
#elif(Bin_size == 2)
#define Di 50
#endif

// step_counting_stage
#define Amin 200 // 200 / FacAmp;
#define Phigh 110 // acf based only
#define Plow 8
#define MinWalkTime 4
#define MinWalkStep 8
#define Fac_acf 10000 // normalize acf values
#define Fac_Step 1000
#define Phys_win_size 1000
#define Num_wave_in_bin 6

#define Dinc 2
#define Win_acf 150 // autocorr window
#define Max_lag 75 // tau
#define Acf_len 76

#if (Bin_size == 1)
#define Win_acfDi 6
#define Ovlp 125 // Win_acf - Di
#define OvlpDi 5 // Ovlp / Di;
#define OvlpBin 10 // 2 * Di; 5 sec
#define MaxCompBin1 40
#define MaxCompBin2 50
#define Win_acc 1125
#define Win_IIR 1175
#define Sec1 2
#define Sec0p5 1
#define Sec2 4
#define Sec3 6
#define Sec4 8
#define Sec5 10
#elif(Bin_size == 2)
#define Win_acfDi 3
#define Ovlp 100 // Win_acf - Di
#define OvlpDi 3 // Ovlp / Di + 1;
#define OvlpBin 5 // 2 * OvlpDi - 1; 5 sec
#define MaxCompBin1 20
#define MaxCompBin2 25
#define Win_acc 1100
#define Win_IIR 1150
#define Sec1 1
#define Sec0p5 1/2
#define Sec2 2
#define Sec3 3
#define Sec4 4
#define Sec5 5
#endif

typedef double RealType;
typedef int intType;

typedef struct _Eul {
	intType nbin;
	RealType *max_acc;
	RealType *min_acc;
	RealType *acc_mean;
	intType *acf_val;
	intType *pitch1;
	intType *pitch2;
	intType *p2ave;
	intType *acf_temp; // 2023_0209
	intType *wd;
	intType *wdWrk1;
	intType *wdWrk2;
	intType *LeftBnd;
	intType *RightBnd;
}EUL;

typedef struct _Sstamp {
	intType *bgn;
	intType *end;
	intType *step;
	intType Glb_cur_tbgn;
	intType Glb_cur_tend;
}SSTAMP;