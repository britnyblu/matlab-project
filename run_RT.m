%this wrapper function sorts the parameters to create RT and merges the
%results
function []=run_RT(folder,interp,smooth,min_len,genome,run_name,cur_dir)

%distinguish between human genome and mouse genomes to choose chrom list
if startsWith(genome,'h')
    chroms=strsplit(num2str(1:22));
else
    chroms=strsplit(num2str(1:19));
end

%add X and Y chroms
temp=strings(size(chroms));
[temp{:}]=chroms{:};
chroms=[temp 'X' 'Y'];

%choose correct gap file
gap=strcat(genome,"_gaps.txt");

%first using the folder get each sample's data and run it in the rt function
%in G1 samples get list of file names
G1_files=dir(strcat(folder,"\G1")); %get file names
full_table=array2table(zeros(0,3),'VariableNames',{'loc','RT','chr'}); %prepare table
for i=3:length(G1_files) %for each file
   sample_name=strsplit(G1_files(i).name,'_G1');
   sample_name=sample_name{1};
   G1_cur_file=strcat(folder,"\G1\",G1_files(i).name);
   S_cur_file=strcat(folder,"\S\",strrep(G1_files(i).name,"G1",'S'));
   
   %now that we have all the parameter we need, run this sample in the
   %helper function
   table=createRT(G1_cur_file,S_cur_file,interp,smooth,min_len,chroms,gap,cur_dir);
   rt_col=3;
   table.Properties.VariableNames{rt_col}=strcat('RT_',strrep(sample_name,'.bed','')); %rename the column with the sample name
   
   %only add the data if it worked, less than 5 means not enough data in
   %the results to include this sample
   if height(full_table)<5
       full_table=table;
   else
       full_table=innerjoin(full_table,table);
   end
end

%save the data
save(strcat(cur_dir,"/profiles/rt_profile_exp_",run_name),'full_table');

end