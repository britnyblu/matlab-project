%%code to make a table listing all available exps and samples
function []=makeDifRegions(group1,group2,table,min_cutoff,output_name,pv_cutoff,run_len,cur_dir)
    %find rel indices
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
    
    
    %iterate to get p values
    Nrows = size(table,1);
    h = zeros(Nrows,1);
    sds=[std(table{:,group1_inds}),std(table{:,group2_inds})];
    df=1;
    for k = 1:Nrows
        h(k) = likelihoodtest(table{k,group1_inds},table{k,group2_inds},sds,df);
    end
    
    
    %add p values and fdr to rel cols
    fdr=bonf_holm(h);
    all_inds=[[1 2] group1_inds group2_inds];
    table=table(:,all_inds);
    table.pv=h;
    table.fdr=fdr;
   % table=[table h fdr];

    
    %label passing cutoff
    pass=fdr<pv_cutoff;
    table.pass_fdr=pass;

    %check that delta is high enough - get inds
    table.delta=mean(table{:,group1_inds},2)-mean(table{:,group2_inds},2);
    table.pass_fdr_delta=(table.pass_fdr & abs(table.delta)>min_cutoff);
    
    %consecutive
    d = [true; diff(table.pass_fdr_delta) ~= 0; true];  % TRUE if values change
    n = diff(find(d));               % Number of repetitions
    Y = repelem(n, n);
    table.run_length=Y;
    table.pass_all=(table.pass_fdr_delta & abs(table.run_length)>=run_len);
    
    %save
    save(strcat(cur_dir,"/dif region files/dif_",output_name),'table');


end




function pv=likelihoodtest(group1,group2,sds,df)
    df=1;
    function prob=probability(sample,mean,sd)
       prob=1/(sqrt(2*pi)*sd)*exp(-(power((sample-mean),2)/(2*sd^2)));
    end
    two_group_likeli=probability(group1,mean(group1),sds(1))+probability(group2,mean(group2),sds(2));
    all_group=[group1 group2];
    one_group_likeli=probability(group1,mean(all_group),sds(1))+probability(group2,mean(all_group),sds(1));
    two=prod(two_group_likeli);
    one=prod(one_group_likeli);
    pv=1-chi2cdf(-2*log(one/two),df);
    end
    
   