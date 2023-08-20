#include "XXAPASC.h"
#pragma warning( disable : 4996)/*no fopen, fscanf warning*/

void print_script(intType *GminTable, RealType *CA, RealType *CB) {

	int ntab, i, fac_filter, nCA, nCB;
	printf("C code Defintion Scripy Begins\n\n");
	printf("#define Fs %d\n", Fs);
	printf("#define Dinc %d\n", Dinc);
	printf("#define Fac_acf %d\n", Fac_acf);
	printf("#define Ovlp %d\n", Ovlp);
	printf("#define Phys_win_size %d\n", Phys_win_size);
	printf("#define Di %d\n", Di);
	printf("#define Win_acf %d\n", Win_acf);
	printf("#define OvlpDi %d\n", OvlpDi);
	printf("#define Win_acfDi %d\n", Win_acfDi);
	printf("#define MaxCompBin1 %d\n", MaxCompBin1);
	printf("#define MaxCompBin2 %d\n", MaxCompBin2);
	printf("#define Max_lag %d\n", Max_lag);
	printf("#define Phigh %d\n", Phigh);
	printf("#define Plow %d\n", Plow);
	printf("#define Num_wave_in_bin %d\n", Num_wave_in_bin);
	printf("#define Amin %d\n", Amin);
	printf("#define MinWalkTime %d\n", MinWalkTime);
	printf("#define MinWalkStep %d\n", MinWalkStep);
	printf("#define Sec1 %d\n", Sec1);
	printf("#define Sec0p5 %d\n", Sec0p5);
	printf("#define Sec3 %d\n", Sec3);
	printf("#define Sec4 %d\n", Sec4);
	printf("#define Sec5 %d\n", Sec5);

	ntab = 62;
	printf("GminTable[%d] = {0,", ntab + 1);
	for (i = 0; i < ntab; i++) {
		if (i < ntab - 1) 
			printf("%d, ", GminTable[i]);
		else 
			printf("%d}\n", GminTable[i]);
	}
	//********************** Print filter coefficient **********************
	fac_filter = 32768;
	nCA = 5;
	printf("int CA[6] = {0, ");
	for (i = 0; i < nCA; i++) {
		if (i < nCA - 1) 
			printf("%d, ", ROUND(fac_filter*CA[i]));
		else 
			printf("%d}\n", ROUND(fac_filter*CA[i]));
	}
	nCB = 5;
	printf("int CB[6] = {0, ");
	for (i = 0; i < nCB; i++) {
		if (i < nCA - 1) 
			printf("%d, ", ROUND(fac_filter*CB[i]));
		else 
			printf("%d}\n", ROUND(fac_filter*CB[i]));
	}
	printf("coefficient have been normalized by fac_filter true for  fac_filter=%d, fs=%d, passband = 0.75~24Hz", fac_filter, Fs);
	printf("\n\n");
	printf("C code Defintion ends\n");
}

int get_data(float *ax, float *ay, float *az) {

	int n = 0, i = 0, j = 0, k = 0;
	FILE *fpRead_a;
	char path[] = "..//..//acc_data//";
	char fileType[] = ".txt";
	char fileName[] = "1_exp-walk-1230-2014//";
	//char fileName[] = "2_exp-wrist-0115-2015//";
	//char fileName[] = "3_exp_positive_0312//";
	//char fileName[] = "4_exp_positive_0318//"; // 0.5noDown_wrong; 0.5Down_wrong; 1noDown_wrong; 1Down_wrong;
	//char fileName[] = "5_exp_0320//"; // 0.5noDown_wrong; 0.5Down_wrong;
	//char fileName[] = "6_exp_0325//";
	//char fileName[] = "7_exp_0331//"; // 0.5Down_wrong; 1Down_wrong;
	//char fileName[] = "8_exp1201-2014//";
	//char fileName[] = "9_paper-data//"; // 1noDown_wrong;
	//char fileName[] = "11_total-ieee//";
	//char fileName[] = "101_exp-non-walk-1230-2014//";
	//char fileName[] = "102_exp_negative_0313//";
	char filePath[200];
	strcpy(filePath, path);
	strcat(filePath, fileName);
	char file_a[300];

	// ax
	strcpy(file_a, filePath);
	strcat(file_a, "x");
	strcat(file_a, fileType);
	fpRead_a = fopen(file_a, "r");
	while (fscanf(fpRead_a, "%f", &ax[i]) != EOF) {
		n++;
		i++;
	}
	fclose(fpRead_a);
	// ay
	strcpy(file_a, filePath);
	strcat(file_a, "y");
	strcat(file_a, fileType);
	fpRead_a = fopen(file_a, "r");
	while (fscanf(fpRead_a, "%f", &ay[j]) != EOF) {
		j++;
	}
	fclose(fpRead_a);
	// az
	strcpy(file_a, filePath);
	strcat(file_a, "z");
	strcat(file_a, fileType);
	fpRead_a = fopen(file_a, "r");
	while (fscanf(fpRead_a, "%f", &az[k]) != EOF) {
		k++;
	}
	fclose(fpRead_a);
	return(n);
}

int get_exp_data(float *x) {

	int n = 0, i = 0;
	FILE *fpRead_a;
	char path[] = "..//..//acc_data//one_hour_sig_2023_0621//";
	//char fileName[] = "cosw.txt";
	//char fileName[] = "ftone.txt";
	char fileName[] = "walk.txt";
	//char fileName[] = "moto.txt";

	char filePath[200];
	strcpy(filePath, path);
	strcat(filePath, fileName);

	fpRead_a = fopen(filePath, "r");
	while (fscanf(fpRead_a, "%f", &x[i]) != EOF) {
		n++;
		i++;
	}
	fclose(fpRead_a);
	return(n);
}