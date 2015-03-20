function ang_hypo = mid_pt_chg(ang, mid_start, flg, var)
% 
% 
% 
% 
% default var = 0.2;

ang_hypo = ang;
ang_hypo(mid_start) = ang(mid_start)+flg*var*0.1;
ang_hypo(mid_start-1) = ang(mid_start-1)-flg*var*0.1;

