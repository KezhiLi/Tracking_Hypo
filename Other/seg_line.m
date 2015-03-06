function seg_pt = seg_line(len, seg_num);

seg_len = round(len/seg_num);

seg_pt = [1: seg_len :len,len]; 

reduc = length(seg_pt)-1-seg_num; 
if reduc > 0
    seg_pt = [seg_pt(1:end-reduc-1),seg_pt(end)];
end



