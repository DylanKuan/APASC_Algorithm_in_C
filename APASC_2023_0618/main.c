#include "APASC.h"
#include <time.h>
// »P step_count_20230310_before ¤@¼Ë

// input
static float ax[1200000] = { 0 };
static float ay[1200000] = { 0 };
static float az[1200000] = { 0 };
static float x[180000] = { 0 }; // one hour

// allocate_Euler_struct
static RealType max_acc[MaxCompBin2] = { 0 };
static RealType min_acc[MaxCompBin2] = { 0 };
static RealType acc_mean[MaxCompBin2] = { 0 };
static intType acf_val[MaxCompBin2] = { 0 };
static intType pitch1[MaxCompBin2] = { 0 };
static intType pitch2[MaxCompBin2] = { 0 };
static intType acf_temp[MaxCompBin2] = { 0 }; // 2023_0209
static intType wd[MaxCompBin2] = { 0 };
static intType wdWrk1[MaxCompBin2] = { 0 };
static intType wdWrk2[MaxCompBin2] = { 0 };
static intType LeftBnd[MaxCompBin2] = { 0 };
static intType RightBnd[MaxCompBin2] = { 0 };

// allocate_Sstamp_struct
static intType bgn[2] = { 0 };
static intType end[2] = { 0 };
static intType step[2] = { 0 };

// IIR filter
static RealType IIR_wrk[Win_IIR] = { 0 };
static RealType CA[5] = { 1.0, -0.044230632246867074, -1.6901413811887278, 0.03228430657945458, 0.73272603037181583 };
static RealType CB[5] = { 0.8559306018148144, 0.0, -1.7118612036296288, 0.0, 0.8559306018148144 };

// acc
static  RealType buffer1[Win_acc] = { 0 };
static  RealType buffer2[Win_acc] = { 0 };
static  RealType acc_temp[Win_acf] = { 0 }; // autocorr // size : Max_Win_acf
static  RealType acf[Acf_len] = { 0 }; // autocorr // size : Max_lag

static intType GminTable[62] = { 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1366, 1234, 1120, 1020, 895, 726, 590, 480, 389, 314, 253, 202, 160, 142, 126, 112, 99, 88, 79, 70, 62, 55, 48, 43, 37, 33, 28, 25, 21, 18, 15, 12, 9, 7, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

void print_script(intType *, RealType *, RealType *);
int get_data(float *, float *, float *);
int get_exp_data(float *);

void main(int argc, char *argv[]) {

	double t1, t2, dt1;
	int nn;
	intType tt, ib = 0, ic = 0, nbin = 0, cwin = 1, firstCal = 1, lasrCal = 0, Win_step, File_step = 0;
	RealType *acc_buf, *acc_comp, *temp;
	intType maxCompBin = MaxCompBin1;
	EUL Eul;
	EUL *Eul_ptr = &Eul;
	SSTAMP Sstamp;
	SSTAMP *Sstamp_ptr = &Sstamp;

	// Euler_struct
	Eul_ptr->nbin = 0;
	Eul_ptr->max_acc = max_acc;
	Eul_ptr->min_acc = min_acc;
	Eul_ptr->acc_mean = acc_mean;
	Eul_ptr->acf_val = acf_val;
	Eul_ptr->pitch1 = pitch1;
	Eul_ptr->pitch2 = pitch2;
	Eul_ptr->acf_temp = acf_temp;
	Eul_ptr->wd = wd;
	Eul_ptr->wdWrk1 = wdWrk1;
	Eul_ptr->wdWrk2 = wdWrk2;
	Eul_ptr->LeftBnd = LeftBnd;
	Eul_ptr->RightBnd = RightBnd;
	// allocate_Sstamp_struct
	Sstamp_ptr->bgn = bgn;
	Sstamp_ptr->end = end;
	Sstamp_ptr->step = step;
	Sstamp_ptr->Glb_cur_tbgn = 0;
	Sstamp_ptr->Glb_cur_tend = 0;
	// buffer
	acc_buf = &buffer1[0];
	acc_comp = &buffer2[0];

	print_script(GminTable, CA, CB); // specify_global_constants_2016_0729

	//nn = get_data(ax, ay, az); // input_accelerometer_data4_new_format
	nn = get_exp_data(x);

	t1 = clock();
	for (tt = 0; tt < nn; tt++) { // global
		ib++; // relative to bin
		ic++; // point relative to comp domain

		//acc_buf[ic - 1] = sqrt(ax[tt] * ax[tt] + ay[tt] * ay[tt] + az[tt] * az[tt]);
		acc_buf[ic - 1] = sqrt(x[tt]);

		if (ib == Di || tt + 1 == nn) {
			nbin++;
			ib = 0;
		}
		if (nbin == maxCompBin || tt + 1 == nn) {
			if (tt == nn) 
				lasrCal = 1;
			memcpy(&acc_comp[0], &acc_buf[ic - Ovlp], Ovlp * sizeof(RealType));
			temp = &acc_comp[0];
			acc_comp = &acc_buf[0];
			acc_buf = &temp[0];

			Win_step = APASC(ic, nbin, acc_comp, IIR_wrk, CA, CB, acf, Eul_ptr, Sstamp_ptr, GminTable, acc_temp, cwin, firstCal, lasrCal);

			//printf("cwin = %d  nbin = %d  maxCompBin=%d  winStep=%d\n", cwin, nbin, maxCompBin, Win_step);
			File_step += Win_step;
			if (tt + 1 == nn) 
				break;
			cwin++;
			firstCal = 0;
			nbin = OvlpBin;
			maxCompBin = MaxCompBin2;
			ic = Ovlp;
		}
	}
	t2 = clock();
	printf("File_step = %d\n", File_step);
	dt1 = (t2 - t1) / CLOCKS_PER_SEC;
	printf("t1 = %6f\tt2 = %6f\tdt1 = %6f\n", t1, t2, dt1);
	return;
}