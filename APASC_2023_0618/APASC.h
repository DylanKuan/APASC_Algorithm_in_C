#include "XXAPASC.h"

intType APASC(intType, intType, RealType *, RealType *, RealType *, RealType *, RealType *, EUL *, SSTAMP *, intType *, RealType *, intType, intType, intType);

void acc_bin_extremum(intType, intType, intType, intType, RealType *, EUL *);

void filter_signal(intType, RealType *, RealType *, RealType *, RealType *, intType);
void IIR_one_Fs_padding(RealType *, RealType *, RealType *, RealType *, intType, intType, RealType);

void acc_bin_mean(intType, EUL *);

void process_similarity_and_pitch(intType, intType, intType, RealType *, RealType *, RealType *, EUL *);
void fast_cal_similarity_0717(intType, intType, intType, RealType *, RealType *, RealType *, EUL *);
void cal_similarity0602_2016(intType, intType, intType, intType, RealType *, RealType *, RealType *, EUL *, RealType);
void autocorr_brute_force_CY(intType, intType, RealType *, RealType *, intType, RealType);
intType step_local_max(intType, intType, RealType *);
intType normalized_acf_2023_0225(RealType, RealType, intType, intType, intType, RealType *);

intType eff_euler_online_step_count_1(intType, intType, intType, intType, intType, intType, EUL *, SSTAMP *, intType *);
void modify_sim_and_pitch_12(intType, EUL *);
intType set_motion_status_6_2(intType, EUL *, intType *);
void modify_irratic_pitch_7(intType, EUL *);
void get_initial_walk_status_3(intType, EUL *, intType *);
void get_initial_walk_status_part_II(intType, intType, intType, intType, EUL *, intType *, intType *);
intType trim_walk_status_8(intType, intType, EUL *, intType *);
intType find_walk_boundary(intType, EUL *);
void kill_short_segment4(intType, intType, intType, EUL *);
intType process_each_window_Euler_only_0625(intType, intType, intType, intType, intType, EUL *, SSTAMP *);
intType segment_step_by_integration_0609_2016(intType, intType, EUL *);
intType dump_small_step_0713(SSTAMP *);
void copy_ghost_bins_2023_0208(intType, EUL *);