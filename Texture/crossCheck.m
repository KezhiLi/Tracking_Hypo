function seq_cor = crossCheck(seq1_ori, seq2_ori, num)
% cross check two vectors, find its best match and how many shifts are
% needed. Here 'match' means that the absolute difference (after adjustment by minus its mean) is minimized
% within the corresponding locations. Do not consider outside values during
% comparision.
% 
% Input: seq1_ori : a vector, 1st sequence added into comparision
%        seq2_ori : a vector, 2nd sequence added into comparision
%        num: how many comparisions are in consideration
%             num should be a odd number eg. default num = 5/7/9;
% Output: seq_cor: a vector of length 'num', save the absolute difference between these two
%                  vectors. The middle value represents the difference without shift.  
%
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 16/06/2015
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary
% notices on any copies of the Software.

% maximum shift to left and right
shift = round((num-1)/2);
% shift vector that controls the shift process
shift_vec = [-shift:shift];

% adjust mean, to focus on 'trend' of the sequence
seq1 = seq1_ori - mean(seq1_ori);
seq2 = seq2_ori - mean(seq2_ori);

% sequence length
len_seq = length(seq1);

seq_cor = zeros(num,1);

% calculate absolute difference on both sides, do not consider values not
% overlap
for ii = 1:num;
    if shift_vec(ii) < 0
        seq_cor(ii) = sum(abs(seq1(1:len_seq+shift_vec(ii))- seq2(-shift_vec(ii)+1:end)))/(len_seq-abs(shift_vec(ii)));
    elseif shift_vec(ii) > 0
        seq_cor(ii) = sum(abs(seq2(1:len_seq-shift_vec(ii))- seq1(shift_vec(ii)+1:end)))/(len_seq-abs(shift_vec(ii)));
    else 
        seq_cor(ii) = sum(abs(seq1- seq2))/len_seq;
    end
end
