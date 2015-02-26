function ang_hypo = mid_pt_chg(ang, mid_start, flg)
% 
% 
% 
% 
% 

ang_hypo = ang;
ang_hypo(mid_start-1) = ang(mid_start-1)+flg*0.2;
ang_hypo(mid_start-1) = ang(mid_start-1)-flg*0.2;

