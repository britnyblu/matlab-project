%% this function is used to create an RT profile from one sample

function data_table=createRT(G1,S,interp,smooth,min_len,chroms,gap_file,cur_dir)
%import data
G1_cell=importdata(G1);
S_cell=importdata(S);
gaps=readtable(strcat(cur_dir,"/gaps/",gap_file));
data_table=array2table(zeros(0,3),'VariableNames',{'loc','RT','chr'});

for i=chroms
    %for this chromosome isolate relevant data
    value=strcat("chr",i);
    rel_G1_indeces=find(G1_cell.textdata(:,1)==value);
    cur_chrom_G1=G1_cell.data(rel_G1_indeces,1);
    rel_S_indeces=find(S_cell.textdata(:,1)==value);
    cur_chrom_S=S_cell.data(rel_S_indeces,1);
    windows_ends=cur_chrom_G1(1:200:end,1);
    
    %loop through S reads to fill out window counts
    if length(windows_ends)<3
       continue 
    end
    max=windows_ends(1,1);
    G1_counter=1;
    S_counter=0;
    for item=1:length(rel_S_indeces)
        if cur_chrom_S(item,1)<max
            S_counter=S_counter+1;
        else
            windows_ends(G1_counter,2)=S_counter;
            S_counter=1;
            if G1_counter>length(windows_ends)-1
               break 
            end
            G1_counter=G1_counter+1;
            max=windows_ends(G1_counter,1);
        end
    end
    
    %make ratios and z scores
    windows_ends(:,3)=windows_ends(:,2)/200;
    windows_ends(:,4)=(windows_ends(:,3)-mean(windows_ends(:,3),1))/std(windows_ends(:,3),1);

    %for each intergap - csaps to smooth
    cur_gaps=gaps(ismember(gaps.Var2,value),:);
    all_smooth_locs=[];
    all_smooth_zscores=[];
    for ind=1:height(cur_gaps)
        inter_gap_start=cur_gaps(ind,4).Var4;
        if ind==height(cur_gaps)
            break
            inter_gap_end=inf;
        else
            inter_gap_end=cur_gaps(ind+1,3).Var3;
        end
        cur_locs=windows_ends(windows_ends(:,1)>=inter_gap_start&windows_ends(:,1)<inter_gap_end,1);
        cur_zscores=windows_ends(windows_ends(:,1)>=inter_gap_start&windows_ends(:,1)<inter_gap_end,4);
        if length(cur_locs)<min_len
           continue 
        end
        F_func = csaps(cur_locs,cur_zscores,smooth);
        %rounding=log10(interp*1000);
        %interp_locs=roundn(inter_gap_start,rounding):interp*1000:round(inter_gap_end,rounding);
        interp_locs=round(inter_gap_start/interp)*interp:interp*1000:round(inter_gap_end/interp)*interp;
        interp_locs=transpose(interp_locs);
        smooth_zscores = fnval(interp_locs,F_func);
        all_smooth_locs=[all_smooth_locs; interp_locs];
        all_smooth_zscores=[all_smooth_zscores; smooth_zscores];
    end
    chro_table=array2table([all_smooth_locs all_smooth_zscores],'VariableNames',{'loc','RT'});
    chro_table.chr(:)=value;
    data_table=[data_table; chro_table];
    
end
data_table = data_table(:,{'chr' 'loc' 'RT'});
data_table=rmoutliers(data_table,'mean','DataVariables','RT');

end
