#include "APASC.h"

// ref : Z.-E Xu, “Accelerometer-based accurate step counting algorithm,” Master, Department of Mechanical Engineering, National Central University, Taoyuan, 2023. https://etd.lib.nctu.edu.tw/cgi-bin/gs32/ncugsweb.cgi?o=dncucdr&s=id=%22GC110323073%22.&searchmode=basic
intType APASC(intType n, intType nbin, RealType *acc_comp, RealType *IIR_wrk, RealType *CA, RealType *CB, RealType *acf, EUL *Eul_ptr, SSTAMP *Sstamp_ptr, intType *GminTable, RealType *acc_temp, intType cwin, intType firstCal, intType lastCal) {

	intType bgnPt_acc, bgnBin_acc, bgnBin_acf, phyBgnBin, phyEndBin, compBgnPt, bgnBin, endBin, Win_step, del_step;
	if (firstCal) {
		bgnPt_acc = 0;
		bgnBin_acc = 0;
		bgnBin_acf = 0;
	}
	else {
		bgnPt_acc = Ovlp;
		bgnBin_acf = OvlpDi;
		if (Di == 25) 
			bgnBin_acc = 2 * OvlpDi;
		else if (Di == 50) 
			bgnBin_acc = 2 * OvlpDi - 1;
	}

	filter_signal(n, acc_comp, IIR_wrk, CA, CB, 4);
	acc_bin_extremum(n, nbin, bgnPt_acc, bgnBin_acc, acc_comp, Eul_ptr);
	acc_bin_mean(nbin, Eul_ptr);
	process_similarity_and_pitch(n, nbin, bgnBin_acf, acc_comp, acf, acc_temp, Eul_ptr);

	Eul_ptr->nbin = nbin;
	phyBgnBin = OvlpDi;
	switch (Di) {
	case 25:
		phyEndBin = Eul_ptr->nbin - OvlpDi - 1;
		compBgnPt = (cwin - 1) * Phys_win_size + 1 - 2 * Ovlp; // in global coord
		break;
	case 50:
		phyEndBin = Eul_ptr->nbin - OvlpDi;
		compBgnPt = (cwin - 1) * Phys_win_size + 1 - (Ovlp + Win_acf); // in global coord
		break;
	}
	bgnBin = 1;
	endBin = Eul_ptr->nbin;
	if (cwin == 1) {
		phyBgnBin = 0;
		compBgnPt = 1;
	}
	if (lastCal) 
		phyEndBin = Eul_ptr->nbin - 1;

	Win_step = eff_euler_online_step_count_1(nbin, phyBgnBin, phyEndBin, bgnBin, endBin, compBgnPt, Eul_ptr, Sstamp_ptr, GminTable);

	return(Win_step);
}

void filter_signal(intType n, RealType *acc_comp, RealType *IIR_wrk, RealType *CA, RealType *CB, intType zzL) {

	intType i;
	RealType sum = 0, ave;
	for (i = 0; i < n; i++)
		sum += acc_comp[i];
	ave = sum / n;
	switch (FilterType) {
	case 0:
		for (i = 0; i < n; i++) 
			acc_comp[i] -= ave;
	case 1:
		IIR_one_Fs_padding(CA, CB, IIR_wrk, acc_comp, n + Fs, 4, ave);
	}
}

void IIR_one_Fs_padding(RealType *CA, RealType *CB, RealType *IIR_wrk, RealType *sig, intType n, intType zzL, RealType ave) { // 20230209

	intType i, j, naxpy, y_tmp, fs = Fs;
	RealType a1 = CA[0], sigi;
	memset(&IIR_wrk[0], 0, n * sizeof(RealType));
	for (i = 0; i < n; i++) {
		if (n - i < zzL + 1) 
			naxpy = (n - 1) - i;
		else 
			naxpy = zzL;

		if (i < fs) 
			sigi = ave;
		else 
			sigi = sig[i - fs];

		for (j = 0; j <= naxpy; j++) {
			y_tmp = i + j;
			IIR_wrk[y_tmp] += sigi * CB[j];
		}
		if ((n - 1) - i < zzL) 
			naxpy = (n - 2) - i;
		else 
			naxpy = zzL - 1;

		a1 = -IIR_wrk[i];
		for (j = 0; j <= naxpy; j++) {
			y_tmp = (i + j) + 1;
			IIR_wrk[y_tmp] += a1 * CA[j + 1];
		}
	}
	memcpy(&sig[0], &IIR_wrk[Fs], (n - Fs) * sizeof(RealType));
}

void acc_bin_extremum(intType n, intType nbin, intType bgnPt_acc, intType bgnBin_acc, RealType *acc, EUL *Eul) {

	intType bin, endPt_acc, i;
	for (bin = bgnBin_acc; bin < nbin; bin++) {
		endPt_acc = MIN(bgnPt_acc + Di, n);
		Eul->max_acc[bin] = acc[bgnPt_acc];
		Eul->min_acc[bin] = acc[bgnPt_acc];
		for (i = bgnPt_acc + 1; i < endPt_acc; i++) {
			if (acc[i] > Eul->max_acc[bin]) 
				Eul->max_acc[bin] = acc[i];
			if (acc[i] < Eul->min_acc[bin]) 
				Eul->min_acc[bin] = acc[i];
		}
		bgnPt_acc += Di;
	}
}

void acc_bin_mean(intType nbin, EUL *Eul) {

	intType jj, jj2, i, ng2;
	RealType local_max, local_min;
	ng2 = MAX(Win_acfDi - 2, 0);
	for (jj = 0; jj < nbin; jj++) {
		jj2 = MIN(jj + ng2, nbin - 1);
		local_max = Eul->max_acc[jj];
		local_min = Eul->min_acc[jj];
		for (i = jj + 1; i <= jj2; i++) {
			if (Eul->max_acc[i] > local_max) 
				local_max = Eul->max_acc[i];
			if (Eul->min_acc[i] < local_min) 
				local_min = Eul->min_acc[i];
		}
		Eul->acc_mean[jj] = local_max - local_min;
	}
}

void process_similarity_and_pitch(intType n, intType nbin, intType bgnBin_acf, RealType *acc, RealType *acf, RealType *acc_temp, EUL *Eul) {

	intType bin, bgn_pt = 0, data_len, max_lag1;
	RealType *acc_ptr;
	data_len = Win_acf;
	for (bin = bgnBin_acf; bin < nbin; bin++) {
		if (Eul->acc_mean[bin] < Amin) {
			Eul->acf_val[bin] = 0;
			Eul->pitch1[bin] = (intType)(4 * Fs);
			Eul->pitch2[bin] = Eul->pitch1[bin];
		}
		else {
			acc_ptr = &acc[bgn_pt]; // the beginning point of the current acc window
			max_lag1 = MIN(data_len - 1, Max_lag);
			fast_cal_similarity_0717(bin, max_lag1, data_len, acf, acc_ptr, acc_temp, Eul);
		}
		bgn_pt += Di;
		data_len = MIN(n - bgn_pt, Win_acf);
	}
}

void fast_cal_similarity_0717(intType bin, intType max_lag0, intType acc0_len, RealType *acf, RealType *acc0, RealType *acc_temp, EUL *Eul) {

	intType i, max_lag, acf_val, pitch1, pitch2, idx, idx_slave, idx_max, idx_nbr1, idx_nbr2, idx_nbr3, dinc;
	RealType ave0 = 0, acf_nbr1 = 0, acf_nbr2 = 0, acf_nbr3 = 0, acf_max;
	dinc = Dinc;
	if (NoDownSample) 
		dinc = 1;
	max_lag = max_lag0 / dinc; // floor(max_lag0/Dinc)
	for (i = 0; i < acc0_len; i += dinc) 
		ave0 += acc0[i];
	ave0 /= (acc0_len/dinc);
	memcpy(&acc_temp[0], &acc0[0], acc0_len * sizeof(RealType)); // 2023_0209
	for (i = 0; i < acc0_len; i++) 
		acc_temp[i] -= ave0;

	cal_similarity0602_2016(bin, max_lag, acc0_len, dinc, acf, acc0, acc_temp, Eul, ave0);

	acf_val = Eul->acf_val[bin];
	pitch1 = Eul->pitch1[bin];
	pitch2 = Eul->pitch2[bin];
	if (acf_val < 2) 
		return;
	if (NoDownSample) 
		return;
	if (pitch1 < 13 && (pitch1 - pitch2) < 3) {
		ave0 = 0;
		for (i = 0; i < acc0_len; i++)
			ave0 += acc0[i];
		ave0 /= acc0_len;
		cal_similarity0602_2016(bin, max_lag0, acc0_len, 1, acf, acc0, acc_temp, Eul, ave0);
		return;
	}
	idx = pitch1 + 1;
	if (pitch2 != pitch1) 
		idx_slave = pitch2 + 1;
	else 
		idx_slave = -1;
	idx_nbr2 = Dinc * idx - 1;
	idx_nbr1 = idx_nbr2 - 1;
	idx_nbr3 = idx_nbr2 + 1;
	for (i = 0; i < acc0_len - idx_nbr1 + 1; i++) 
		acf_nbr1 += acc_temp[i] * acc_temp[idx_nbr1 - 1 + i];
	for (i = 0; i < acc0_len - idx_nbr2 + 1; i++) 
		acf_nbr2 += acc_temp[i] * acc_temp[idx_nbr2 - 1 + i];
	for (i = 0; i < acc0_len - idx_nbr3 + 1; i++) 
		acf_nbr3 += acc_temp[i] * acc_temp[idx_nbr3 - 1 + i];

	acf_max = acf_nbr1;
	idx_max = 1;
	if (acf_nbr2 > acf_max) {
		idx_max = 2;
		acf_max = acf_nbr2;
	}
	if (acf_nbr3 > acf_max) {
		idx_max = 3;
		acf_max = acf_nbr3;
	}
	pitch1 = idx_nbr2 + idx_max - 2;
	if (idx_slave == -1) 
		pitch2 = pitch1;
	else 
		pitch2 = Dinc * idx_slave - 1;
	pitch1 = pitch1 - 1;
	pitch2 = pitch2 - 1;

	Eul->acf_val[bin] = acf_val;
	Eul->pitch1[bin] = pitch1;
	Eul->pitch2[bin] = pitch2;
}

void cal_similarity0602_2016(intType bin, intType max_lag, intType acc0_len, intType dinc, RealType *acf, RealType *acc0, RealType *acc_temp, EUL *Eul, RealType ave) {

	intType i, nacf, half_n, idx, pitch1, pitch2, jdx1, jdx2, idx_slave, acf_val;
	RealType acf00, acf_val_max;
	nacf = max_lag + 1;
	autocorr_brute_force_CY(max_lag, acc0_len, acc_temp, acf, dinc, ave);
	idx = step_local_max(0, nacf, acf);

	if (idx == -1) 
		acf_val_max = 0;
	else 
		acf_val_max = acf[idx];
	if (idx > -1) {
		if (10 * (idx + 1) < 35 * Fs && acf[MAX(idx - 3, 0)] > acf[idx]) 
			idx = step_local_max(MIN(idx + 1, max_lag - 1), nacf, acf);
	}
	pitch1 = idx;
	pitch2 = idx;
	if (idx <= 0 || nacf < 10) {
		Eul->acf_val[bin] = 0;
		Eul->pitch1[bin] = 4 * Fs;
		Eul->pitch2[bin] = Eul->pitch1[bin];
		return;
	}
	acf_val_max = acf[idx];
	acf00 = acf[0];

	acf_val = normalized_acf_2023_0225(acf00, acf_val_max, idx, acc0_len, dinc, acc0);

	jdx1 = MAX(Num_wave_in_bin + 1, (idx + 1) * 3 / 8 - 1);
	jdx2 = (idx + 1) * 5 / 8;
	jdx2 = MIN(jdx2 + 1, max_lag);
	idx_slave = step_local_max(jdx1, jdx2, acf);
	if (idx_slave > jdx1) 
		pitch2 = idx_slave;
	Eul->acf_val[bin] = acf_val;
	Eul->pitch1[bin] = pitch1;
	Eul->pitch2[bin] = pitch2;
}

void autocorr_brute_force_CY(intType tau, intType n, RealType *sig, RealType *acf, intType dinc, RealType ave) {

	intType i, j, k;
	RealType sum;
	for (i = 0, j = 0; i <= tau; i++, j += dinc) {
		sum = 0;
		for (k = 0; j + k < n; k += dinc) {
			sum += sig[k] * sig[j + k];
		}
		acf[i] = sum;
	}
}

intType step_local_max(intType i1, intType i2, RealType *x) {

	intType i, idx = -1;
	RealType init_max = 0;
	for (i = i1 + 1; i < i2 - 1; i++) {
		if (x[i] >= x[i + 1] && x[i] >= x[i - 1] && (2 * x[i] > x[i + 1] + x[i - 1])) {
			if (x[i] > init_max) {
				init_max = x[i];
				idx = i;
			}
		}
	}
	return(idx);
}

intType normalized_acf_2023_0225(RealType acf00, RealType acf_val_max, intType idx, intType acc0_len, intType dinc, RealType *acc0) {

	intType i, j, k, n, len, v0, acf_val;
	RealType nA2 = 0, nB2 = 0, sumAB = 0, nA, nB, max2, min2;
	if (10 * acf_val_max < (2 * acf00)) {
		acf_val = (intType)(Fac_acf * acf_val_max / acf00);
		return(acf_val);
	}

	j = 0;
	k = idx * dinc;
	n = ROUND(1.0 * acc0_len / dinc); //********************************************************************************* if acc0_len = 25 dinc = 2 bug, add 1.0
	len = n - idx;

	for (i = 0; i < len; i++) {
		nA2 += acc0[j] * acc0[j];
		nB2 += acc0[k] * acc0[k];
		sumAB += acc0[j] * acc0[k];
		j += dinc;
		k += dinc;
	}
	nA = sqrt(nA2);
	nB = sqrt(nB2);
	max2 = nA;
	min2 = nB;
	if (nB > nA) {
		max2 = nB;
		min2 = nA;
	}
	if (1) {
		v0 = (intType)(Fac_acf * sumAB / (nA * nB));
		if (max2 >= 3 * min2)
			acf_val = (intType)((147 * v0) / 256);
		else if (max2 >= 2 * min2)
			acf_val = (intType)((-(intType)(33 * max2 / min2) + 247) * v0 / 256);
		else if (max2 >= min2)
			acf_val = (intType)((-(intType)(75 * max2) / min2 + 331) * v0 / 256);
		else
			acf_val = (intType)(v0);
	}
	else {
		// coe_max = 1.5; coe_min = 0.5;
		acf_val = (intType)(Fac_acf * sumAB * sqrt(min2) / (nA * nB * sqrt(max2)));
	}
	return(acf_val);
}

intType eff_euler_online_step_count_1(intType nbin, intType phyBgnBin, intType phyEndBin, intType bgnBin, intType endBin, intType compBgnPt, EUL *Eul0, SSTAMP *Sstamp_ptr, intType *GminTable) {

	intType numBnd, winStep = 0;
	memcpy(&Eul0->acf_temp[0], &Eul0->acf_val[nbin - OvlpBin], OvlpDi * sizeof(intType));
	memcpy(&Eul0->acf_temp[OvlpDi], &Eul0->pitch1[nbin - OvlpBin], OvlpDi * sizeof(intType));
	memcpy(&Eul0->acf_temp[2 * OvlpDi], &Eul0->pitch2[nbin - OvlpBin], OvlpDi * sizeof(intType));

	modify_sim_and_pitch_12(nbin, Eul0);
	numBnd = set_motion_status_6_2(nbin, Eul0, GminTable);
	winStep = process_each_window_Euler_only_0625(numBnd, phyBgnBin, phyEndBin, compBgnPt, winStep, Eul0, Sstamp_ptr);
	copy_ghost_bins_2023_0208(nbin, Eul0);
	return(winStep);
}

void modify_sim_and_pitch_12(intType nbin, EUL *Eul0) {

	intType i, jj, jj0, jj1, jj2, thresh2, jn1;
	RealType local_max, local_min, amp0, amp1;
	for (jj = 0; jj < nbin; jj++) {
		jj0 = MAX((jj - Sec0p5), 0);
		local_max = Eul0->max_acc[jj0];
		local_min = Eul0->min_acc[jj0];
		for (i = jj0 + 1; i <= jj; i++) {
			if (local_max < Eul0->max_acc[i]) 
				local_max = Eul0->max_acc[i];
			if (local_min > Eul0->min_acc[i]) 
				local_min = Eul0->min_acc[i];
		}
		amp0 = local_max - local_min;

		jj1 = MIN(jj + 1, nbin - 1);
		jj2 = MIN(jj + 1 + Sec1, nbin - 1);
		local_max = Eul0->max_acc[jj1];
		local_min = Eul0->min_acc[jj1];
		for (i = jj1 + 1; i <= jj2; i++) {
			if (local_max < Eul0->max_acc[i]) 
				local_max = Eul0->max_acc[i];
			if (local_min > Eul0->min_acc[i]) 
				local_min = Eul0->min_acc[i];
		}
		amp1 = local_max - local_min;

		if (amp0 < 300) 
			thresh2 = 5;
		else 
			thresh2 = 2;

		if (10 * amp0 >= thresh2 * amp1) {
			if (Eul0->acc_mean[jj] > 1800) {
				if (Eul0->pitch2[jj] < 16) 
					Eul0->acf_val[jj] = MIN(Eul0->acf_val[jj] + (Fac_acf * 0.0938), 10000); // 16 = floor(Fs / 3); Fs = 50; // 0.0938 = 6/64;
				else 
					Eul0->acf_val[jj] = MIN(Eul0->acf_val[jj] + (Fac_acf * 0.0625), 10000); // 0.0625 = 1/16;
			}
		}
		else {
			if (amp0 < 300) {
				Eul0->acf_val[jj] = 0;
				Eul0->pitch1[jj] = 4 * Fs;
				Eul0->pitch2[jj] = Eul0->pitch1[jj];
			}
			else {
				jn1 = MAX(jj - 1, 0);
				Eul0->acf_val[jj] = Eul0->acf_val[jn1];
				Eul0->pitch1[jj] = Eul0->pitch1[jn1];
				Eul0->pitch2[jj] = Eul0->pitch2[jn1];
			}
		}
	}
}

intType set_motion_status_6_2(intType nbin, EUL *Eul0, intType *GminTable) {

	intType numBnd;
	modify_irratic_pitch_7(nbin, Eul0);
	get_initial_walk_status_3(nbin, Eul0, GminTable);
	numBnd = trim_walk_status_8(nbin, Sec5, Eul0, GminTable);
	return(numBnd);
}

void modify_irratic_pitch_7(intType nbin, EUL *Eul0) {

	intType i, j, maxBin;
	for (i = 0; i < nbin; i++) {
		if (Eul0->pitch1[i] >= 4 * Fs) {
			maxBin = MIN(nbin, i + Sec4 + 1);
			for (j = i + 1; j < maxBin; j++) {
				if (Eul0->pitch1[j] < 4 * Fs) {
					Eul0->pitch1[i] = Eul0->pitch1[j];
					Eul0->pitch2[i] = Eul0->pitch2[j];
					break;
				}
			}
		}
	}
}

void get_initial_walk_status_3(intType nbin, EUL *Eul0, intType *GminTable) {

	intType idx, i, j, ip1, i1, i2, ip2 = -1, span, pp, c1, c2, thresh1, cond1, local_max;
	for (i = 0; i < nbin; i++) {
		ip1 = 0;
		if (i >= nbin - 3) 
			ip1 = MAX(nbin - 3, 1);
		i1 = MAX(0, i - 1);
		if (10 * Eul0->acc_mean[i] < 5 * Eul0->acc_mean[i1]) 
			ip1 = i1 + 1;
		i2 = MAX(0, i - 2);
		if (10 * Eul0->acc_mean[i] < 3 * Eul0->acc_mean[i2] && Eul0->acc_mean[i1] < Eul0->acc_mean[i2]) 
			ip1 = i2 + 1;
		if (ip1 > 0) {
			ip2 = MAX(1, i - Win_acfDi + 2);
			if (ip2 > ip1) 
				continue;
			local_max = Eul0->acf_val[ip2 - 1];
			for (j = ip2; j < ip1; j++) {
				if (local_max < Eul0->acf_val[j]) {
					local_max = Eul0->acf_val[j];
					ip2 = j + 1;
				}
			}
			for (j = i; j > ip2 - 1; j--) {
				Eul0->pitch1[j] = Eul0->pitch1[ip2 - 1];
				Eul0->pitch2[j] = Eul0->pitch2[ip2 - 1];
				Eul0->acf_val[j] = Eul0->acf_val[ip2 - 1];
			}
		}
	}

	get_initial_walk_status_part_II(0, nbin, 1, 1, Eul0, GminTable, Eul0->wd);

	span = 2;
	for (i = 0; i < nbin; i++) {
		pp = Eul0->pitch2[i];
		i1 = MAX(i - span, 0);
		i2 = MIN(i + span, nbin - 1);
		c1 = 0;
		c2 = 0;
		thresh1 = 2;
		for (j = i1; j <= i; j++) {
			if (10 * ABS(Eul0->pitch2[j] - pp) < thresh1 * pp) 
				c1++;
		}
		for (j = i; j <= i2; j++) {
			if (10 * ABS(Eul0->pitch2[j] - pp) < 2 * thresh1 * pp) 
				c2++;
		}
		cond1 = ((c1 + c2) >= 3);
		if (cond1 == 0 && Eul0->wd[i] > 1) Eul0->wd[i] = 0;
	}
}

void get_initial_walk_status_part_II(intType j1, intType j2, intType icase, intType floating_value, EUL *Eul0, intType *GminTable, intType *wd) {

	intType jj, pp, st_pitch, st_acc, st_simi;
	RealType gmin;  //********************************************************************************* intType?
	for (jj = j1; jj < j2; jj++) {
		st_pitch = 2 * (Eul0->pitch2[jj] < Phigh && Eul0->pitch2[jj] > Plow);
		pp = Eul0->pitch2[jj];
		pp = MAX(pp, 1);
		pp = MIN(pp, 60) - 1; // 60 = ceil(1.2 * Fs)
		if (icase == 1) 
			gmin = GminTable[pp] + Amin; // gmin = GlbCnt.GminTable(pp)/GlbCnt.FacAmp + GlbCnt.Amin;
		else 
			gmin = GminTable[pp] / 2 + Amin; // gmin = GlbCnt.GminTable(pp)/GlbCnt.FacAmp + GlbCnt.Amin;

		st_acc = 2 * (Eul0->acc_mean[jj] > gmin);   //********************************************************************************* intType?

		if (2 * Eul0->acf_val[jj] > Fac_acf) {
			st_simi = 2;
		}
		else if (10 * Eul0->acf_val[jj] < Fac_acf * 4) {
			st_simi = 0;
		}
		else {
			if (Eul0->acc_mean[jj] > 2000) 
				st_simi = 2;
			else 
				st_simi = 1;
		}
		if (st_acc == 0 || st_pitch == 0 || st_simi == 0) 
			wd[jj] = 0;
		else if (st_acc == 2 && st_pitch == 2 && st_simi == 2) 
			wd[jj] = 2;
		else 
			wd[jj] = floating_value;
	}
}

intType trim_walk_status_8(intType nbin, intType winSize, EUL *Eul0, intType *GminTable) {

	intType i, j, k, iloop, numBnd, jbgn, jend, majority, span, prob1, j1, j2, cc0, cc1, cc2, c1, c2, i1, i2, maxNum, num, jlast, acceptableNum;
	for (i = 0; i < nbin; i += winSize) {
		jbgn = i;
		jend = MIN(jbgn + winSize, nbin);
		cc0 = 0;
		cc1 = 0;
		cc2 = 0;
		for (j = jbgn; j < jend; j++) {
			if (Eul0->wd[j] == 0) 
				cc0++;
			else if (Eul0->wd[j] == 1) 
				cc1++;
			else 
				cc2++;
		}
		if (cc1 == 0) 
			continue;
		majority = 2 * (cc0 <= cc2);
		for (j = jbgn; j < jend; j++) {
			if (Eul0->wd[j] == 1) 
				Eul0->wd[j] = majority;
		}
	}
	for (iloop = 1; iloop <= 2; iloop++) {
		memcpy(&Eul0->wdWrk1[0], &Eul0->wd[0], nbin * sizeof(intType));
		if (Di == 25) {
			if (iloop == 1) {
				span = 10;
				prob1 = 41;
			}
			else {
				span = 2;
				prob1 = 51;
			}
		}
		else if (Di == 50) {
			if (iloop == 1) {
				span = 5;
				prob1 = 21;
			}
			else {
				span = 1;
				prob1 = 41;
			}
		}
		for (i = 0; i < nbin; i++) {
			if (Eul0->wd[i] == 0) 
				continue;
			j1 = MAX(i - span, 0);
			j2 = MIN(i + span, nbin - 1);
			c1 = 0; c2 = 0;
			for (j = j1; j <= i; j++) 
				c1 += (Eul0->wd[j] > 0);
			for (j = i; j <= j2; j++) 
				c2 += (Eul0->wd[j] > 0);
			if ((100 * c1) <= (prob1 * span) && (100 * c2) <= (prob1 * span)) 
				Eul0->wdWrk1[i] = 0;
		}
		memcpy(&Eul0->wd[0], &Eul0->wdWrk1[0], nbin * sizeof(intType));
	}

	numBnd = find_walk_boundary(nbin, Eul0);

	for (i = 0; i < numBnd; i++) {
		i1 = MIN(Eul0->RightBnd[i] + 1, nbin - 1);
		if (Di == 25) 
			i2 = MIN(Eul0->RightBnd[i] + 9, nbin);
		else if (Di == 50) 
			i2 = MIN(Eul0->RightBnd[i] + 5, nbin);
		maxNum = i2 - i1;
		memcpy(&Eul0->wdWrk1[i1], &Eul0->wd[i1], maxNum * sizeof(intType)); // wdWrk1(i1:i2) = wd(i1:i2);

		get_initial_walk_status_part_II(i1, i2, 2, 0, Eul0, GminTable, Eul0->wdWrk1);

		num = 0;
		jlast = -1;
		for (j = i1; j < i2; j++) {
			if (Eul0->wdWrk1[j] == 2) {
				num++;
				jlast = j;
			}
		}
		if (Di == 25) 
			acceptableNum = 2;
		else if (Di == 50) 
			acceptableNum = 1;
		if (num >= maxNum - acceptableNum) {
			for (j = i1; j <= jlast; j++) {
				Eul0->wd[j] = 2;
			}
		}
	}
	numBnd = find_walk_boundary(nbin, Eul0);
	kill_short_segment4(nbin, Sec1, Sec1, Eul0);
	kill_short_segment4(nbin, Sec4, Sec3, Eul0);
	numBnd = find_walk_boundary(nbin, Eul0);
	return(numBnd);
}

intType find_walk_boundary(intType nbin, EUL *Eul0) {

	intType i, c1 = -1, numBnd = 0;
	if (Eul0->wd[0] > 0) {
		Eul0->LeftBnd[0] = 0;
		c1 = 0;
	}
	for (i = 1; i < nbin; i++) {
		if (Eul0->wd[i] > 0 && Eul0->wd[i - 1] == 0) {
			c1++;
			Eul0->LeftBnd[c1] = i;
		}
	}

	for (i = 0; i < nbin - 1; i++) {
		if (Eul0->wd[i] > 0 && Eul0->wd[i + 1] == 0) {
			Eul0->RightBnd[numBnd] = i;
			numBnd++;
		}
	}
	if (Eul0->wd[nbin - 1] > 0) {
		Eul0->RightBnd[numBnd] = nbin - 1;
		numBnd++;
	}
	return(numBnd);
}

void kill_short_segment4(intType nbin, intType thresh1, intType thresh2, EUL *Eul0) {

	intType i, r, nr = 0, prevStatus, curStatus;
	memcpy(&Eul0->wdWrk1[0], &Eul0->LeftBnd[0], (10) * sizeof(intType));
	memcpy(&Eul0->wdWrk2[0], &Eul0->RightBnd[0], (10) * sizeof(intType));

	Eul0->wdWrk1[0] = 1;
	for (i = 1; i < nbin; i++) {
		if (Eul0->wd[i] != Eul0->wd[i - 1]) {
			Eul0->wdWrk2[nr] = i;
			nr++;
			Eul0->wdWrk1[nr] = i + 1;
		}
	}
	Eul0->wdWrk2[nr] = nbin;
	prevStatus = Eul0->wd[Eul0->wdWrk1[0] - 1];
	for (r = 1; r < nr; r++) {
		curStatus = Eul0->wd[Eul0->wdWrk1[r] - 1];
		if (prevStatus != curStatus && (Eul0->wdWrk2[r + 1] - Eul0->wdWrk1[r + 1]) >= thresh2 && (Eul0->wdWrk2[r - 1] - Eul0->wdWrk1[r - 1]) >= thresh2 && (Eul0->wdWrk2[r] - Eul0->wdWrk1[r]) <= thresh1) {
			for (i = Eul0->wdWrk1[r]; i <= Eul0->wdWrk2[r]; i++) 
				Eul0->wd[i - 1] = prevStatus;
		}
		prevStatus = curStatus;
	}
}

intType process_each_window_Euler_only_0625(intType numSeg, intType phyBgnBin, intType phyEndBin, intType compBgnPt, intType winStep, EUL *Eul0, SSTAMP *Sstamp) {

	intType s, jbgn, jend, binBgn, binEnd, segStep, del_step;
	if (numSeg == 0) 
		return(0);
	for (s = 0; s < numSeg; s++) {
		jbgn = Eul0->LeftBnd[s];
		jend = Eul0->RightBnd[s];
		binBgn = MAX(phyBgnBin, jbgn);
		binEnd = MIN(phyEndBin, jend);

		segStep = segment_step_by_integration_0609_2016(binBgn, binEnd, Eul0);
		winStep += segStep;

		Sstamp->Glb_cur_tbgn = compBgnPt + binBgn * Di;
		Sstamp->Glb_cur_tend = compBgnPt + (binEnd + 1) * Di - 1;
		Sstamp->bgn[1] = Sstamp->Glb_cur_tbgn;
		Sstamp->end[1] = Sstamp->Glb_cur_tend;
		Sstamp->step[1] = segStep;

		del_step = dump_small_step_0713(Sstamp);
		winStep -= del_step;
	}
	return(winStep);
}

intType segment_step_by_integration_0609_2016(intType binBgn, intType binEnd, EUL *Eul0) {

	intType i, freq, sumStep = 0;
	for (i = binBgn; i <= binEnd; i++) {
		freq = Fac_Step * Fs / Eul0->pitch1[i];
		if (Eul0->pitch1[i] != Eul0->pitch2[i]) 
			freq = 2 * freq;
		sumStep += freq;
	}
	sumStep = ROUND(1.0 * sumStep * Di / (Fac_Step * Fs)); // 0.0005 = Di / (Fac_Step * Fs); //********************************************************************************* intType?
	return(sumStep);
}

intType dump_small_step_0713(SSTAMP *Sstamp) {

	intType is_small_step, del_step = 0;
	if ((Sstamp->bgn[1] - Sstamp->end[0]) > Fs) {
		is_small_step = (Sstamp->end[0] - Sstamp->bgn[0]) < Fs * MinWalkTime || Sstamp->step[0] < MinWalkStep;
		if (is_small_step) 
			del_step = Sstamp->step[0];
		Sstamp->bgn[0] = Sstamp->bgn[1];
		Sstamp->end[0] = Sstamp->end[1];
		Sstamp->step[0] = Sstamp->step[1];
	}
	else {
		Sstamp->end[0] = Sstamp->end[1];
		Sstamp->step[0] += Sstamp->step[1];
	}
	return(del_step);
}

void copy_ghost_bins_2023_0208(intType nbin, EUL *Eul) {

	memcpy(&Eul->acf_val[0], &Eul->acf_temp[0], OvlpDi * sizeof(intType));
	memcpy(&Eul->pitch1[0], &Eul->acf_temp[OvlpDi], OvlpDi * sizeof(intType));
	memcpy(&Eul->pitch2[0], &Eul->acf_temp[2 * OvlpDi], OvlpDi * sizeof(intType));
	memcpy(&Eul->max_acc[0], &Eul->max_acc[nbin - OvlpBin], OvlpBin * sizeof(RealType));
	memcpy(&Eul->min_acc[0], &Eul->min_acc[nbin - OvlpBin], OvlpBin * sizeof(RealType));
}