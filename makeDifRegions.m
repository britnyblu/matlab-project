%%code to make a table listing all available exps and samples
function []=makeDifRegions(group1,group2,table,min_cutoff,output_name,pv_cutoff,run_len,cur_dir)
    %find rel indices of the groups
    group1_inds=[];
    for i=group1
        id=find(strcmpi(table.Properties.VariableNames,i));
        group1_inds=[group1_inds id];
    end
    group2_inds=[];
    for i=group2
        id=find(strcmpi(table.Properties.VariableNames,i));
        group2_inds=[group2_inds id];
    end    
    
    %prepare each group's sd of differences between value to window mean    
    sds=[];
    table.group1_means=mean(table{:,group1_inds},2);
    table.group2_means=mean(table{:,group2_inds},2);
    dif_to_mean=[];
    for i=1:length(group1_inds)
        dif_to_mean=[dif_to_mean; table.(group1_inds(i))-table.group1_means];
    end
    sds=std(dif_to_mean);
    dif_to_mean=[];
    for i=1:length(group2_inds)
        dif_to_mean=[dif_to_mean; table.(group2_inds(i))-table.group2_means];
    end
    sds=[sds std(dif_to_mean)];
    %iterate to get p value for each window
    df=1;
    Nrows = size(table,1);
    h = zeros(Nrows,1);
    for k = 1:Nrows
        h(k) = likelihoodtest(table{k,group1_inds},table{k,group2_inds},sds,df);
    end
    
    
    %add p values and fdr to table
    fdr=bonf_holm(h);
    [h_fdr, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(h);
    table.pv=h;
    table.fdr=adj_p;
    
    %label passing fdr cutoff
    pass=adj_p<pv_cutoff; 
    table.pass_fdr=pass;

    %label passing minimum delta cutoff
    table.delta=mean(table{:,group1_inds},2)-mean(table{:,group2_inds},2);
    table.pass_delta=abs(table.delta)>min_cutoff;
    table.pass_fdr_delta=(table.pass_fdr & abs(table.delta)>min_cutoff);
    
    %table(1:10,:)
    %sum(table.pass_fdr==true)
    
    %label consecutive regions and cutoff
    d = [true; diff(table.pass_fdr_delta) ~= 0; true];  % TRUE if values change
    n = diff(find(d));               % Number of repetitions
    Y = repelem(n, n);
    table.run_length=Y; %this shows how many windows in a row behave the same
    table.pass_all=(table.pass_fdr_delta & abs(table.run_length)>=run_len); %this takes all windows that pass both cutoffs (fdr and delta) for enough consecutive windows 
    
    %all_inds=[[1 2] group1_inds group2_inds];
    %table=table(:,all_inds);

    
    %save
    save(strcat(cur_dir,"/dif region files/dif_",output_name),'table');


end



%%This is the likelihood test
function pv=likelihoodtest(group1,group2,sds,df)
    df=1; %set degrees freedom
    
    %%This calculates the probability of these values belonging to a group
    %%with the given mean and sd
    function prob=probability(sample,mean,sd)
       prob=1/(sqrt(2*pi)*sd)*exp(-(power((sample-mean),2)/(2*sd^2)));
    end
    
    %likelihood we have two different groups
    two_group_likeli=1;
    for i=1:length(group1)
        two_group_likeli=two_group_likeli*probability(group1(i),mean(group1),sds(1));
    end
    for i=1:length(group2)
        two_group_likeli=two_group_likeli*probability(group2(i),mean(group2),sds(2));
    end    
  
    all_group=[group1 group2];
    
    %likelihood all the points actually belong to the same group
    one_group_likeli=1;
    for i=1:length(group1)
        one_group_likeli=one_group_likeli*probability(group1(i),mean(all_group),sds(1));
    end
    for i=1:length(group2)
        one_group_likeli=one_group_likeli*probability(group2(i),mean(all_group),sds(2));
    end    
    two=two_group_likeli;
    one=one_group_likeli;
    
    %likelihood ratio test and resulting p value
    pv=1-chi2cdf(-2*log(one/two),df);
    end
    
   