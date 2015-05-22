function seq_cor = crossCheck(seq1, seq2, num)
% 
% 
% 
% 
% 
% 
% 
% 
%
% num should be a odd number eg. default num = 5;


shift = round((num-1)/2);
shift_vec = [-shift:shift];

len_seq = length(seq1);

seq_cor = zeros(num,1);

for ii = 1:num;
    if shift_vec(ii) < 0
        seq_cor(ii) = sum(abs(seq1(1:len_seq+shift_vec(ii))- seq2(-shift_vec(ii)+1:end)))/(len_seq-abs(shift_vec(ii)));
    elseif shift_vec(ii) > 0
        seq_cor(ii) = sum(abs(seq2(1:len_seq-shift_vec(ii))- seq1(shift_vec(ii)+1:end)))/(len_seq-abs(shift_vec(ii)));
    else 
        seq_cor(ii) = sum(abs(seq1- seq2))/len_seq;
    end
end
