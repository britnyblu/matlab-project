%this wrapper function sorts the parameters to create RT and merges the
%results
function []=run_RT(folder,interp,smooth,min_len,genome,run_name,cur_dir)
if startsWith(genome,'h')
    chroms=strsplit(num2str(1:22));
else
    chroms=strsplit(num2str(1:19));
end
temp=strings(size(chroms));
[temp{:}]=chroms{:};
chroms=[temp 'X' 'Y'];
gap=strcat(genome,"_gaps.txt");
%first using the folder get each samples data and run it in the rt function
%in G1 samples get list of file names
G1_files=dir(strcat(folder,"\G1"));
full_table=array2table(zeros(0,3),'VariableNames',{'loc','RT','chr'});
for i=3:length(G1_files)
   sample_name=strsplit(G1_files(i).name,'_G1');
   sample_name=sample_name{1};
   G1_cur_file=strcat(folder,"\G1\",G1_files(i).name);
   S_cur_file=strcat(folder,"\S\",strrep(G1_files(i).name,"G1",'S'));
   table=createRT(G1_cur_file,S_cur_file,interp,smooth,min_len,chroms,gap,cur_dir);
   rt_col=3;
   table.Properties.VariableNames{rt_col}=strcat('RT_',sample_name);
   if height(full_table)<5
       full_table=table;
   else
       full_table=innerjoin(full_table,table);
   end
end

save(strcat(cur_dir,"/profiles/rt_profile_exp_",run_name),'full_table');

end