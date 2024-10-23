%% Now...is there an age effect?
idx_pre = fMRI_hmm4_music.pre_dex==1;
idx_post = fMRI_hmm4_music.pre_dex==0;
idx_old = fMRI_hmm4_music.Ages>40;
idx_young = fMRI_hmm4_music.Ages<40;
%% Slice it up
Pre_FO = fMRI_hmm4_music.block_FO(idx_pre==1);
Post_FO = fMRI_hmm4_music.block_FO(idx_post==1);
Pre_Stim = fMRI_hmm4_music.stim_FO(idx_pre==1,:);
Post_Stim = fMRI_hmm4_music.stim_FO(idx_post==1,:);
Pre_FO_young = fMRI_hmm4_music.block_FO(idx_pre==1&idx_young==1);
Post_FO_young = fMRI_hmm4_music.block_FO(idx_post==1&idx_young==1);
Pre_FO_old = fMRI_hmm4_music.block_FO(idx_pre==1&idx_old==1);
Post_FO_old = fMRI_hmm4_music.block_FO(idx_post==1&idx_old==1);
Pre_Stim_young = fMRI_hmm4_music.stim_FO(idx_pre==1&idx_young==1,:);
Post_Stim_young = fMRI_hmm4_music.stim_FO(idx_post==1&idx_young==1,:);
Pre_Stim_old = fMRI_hmm4_music.stim_FO(idx_pre==1&idx_old==1,:);
Post_Stim_old = fMRI_hmm4_music.stim_FO(idx_post==1&idx_old==1,:);
%% Appetizer PLS
addpath(genpath('Pls'))
FO_data = fMRI_hmm4_music.Metrics.FO;
indata = {FO_data(idx_pre==1&idx_young==1);FO_data(idx_pre==1&idx_old==1)};
numsubs = nan(1,length(indata));
for i = 1:length(indata)
    tempdata = indata{i};
    tempdata(isnan(tempdata)) = [];
    indata{i} = tempdata;
    numsubs(i) = size(tempdata,1);
end
clear i option
option.method = 1;
option.num_perm = 500;
option.num_boot = 100;

app_res = pls_analysis(indata,numsubs,1,option);
%% Average it down for PLS:
Big_Mat = {Pre_Stim_young;Pre_Stim_old;Post_Stim_old};
Age_Task_Sess = cell(3,3);
%Pre_YMat = nan(size(Pre_Stim_young));
%Pre_OMat = nan(size(Pre_Stim_old));
%Post_OMat = nan(size(Post_Stim_old));
K = 4;
for mat = 1:length(Big_Mat)
    tempdata = Big_Mat{mat};
    for cat = 1:size(tempdata,2) %category (ss,popular, novel)
        submat = tempdata(:,cat);% get the relevant row
        outmat = nan(length(tempdata),K);%pre-allocate the output matrix
        for part = 1:length(outmat) %rows, or participants
            tempmat = submat{part}; %get the matrix
            if isempty(tempmat)
                continue
            else
                outmat(part,:) = mean(tempmat);%take the average
            end
        end
        Age_Task_Sess{cat,mat} = outmat;%save the average data
    end
end
clear Big_Mat K mat tempdata cat submat outmat tempmat part 
%% And now some PLSes!
addpath(genpath('Pls'))
indata = {cell2mat(Age_Task_Sess(:,1));cell2mat(Age_Task_Sess(:,2));...
    cell2mat(Age_Task_Sess(:,3))};
numcond = 3;
numsubs = zeros(1,length(indata));

for i = 1:length(indata)
    tempdata = indata{i};
    tempdata(any(isnan(tempdata),2),:)=[];
    indata{i} = tempdata;
    numsubs(i) = length(tempdata)/numcond;
end
clear i option

option.method = 1;
option.num_perm = 500;
option.num_boot = 100;

%presess_res = pls_analysis(indata(1:2),numsubs(1:2),numcond,option);
%oldsess_res = pls_analysis(indata(2:end),numsubs(2:end),numcond,option);
session_res = pls_analysis(indata,numsubs,numcond,option);
%% And now some PLSes!
addpath(genpath('Pls'))
%indata = fMRI_hmm4_music.PrePost_PLS.input(:,1);%young
indata = fMRI_hmm4_music.PrePost_PLS.input(:,2);%old

numcond = 3;

for i = 1:length(indata)
    tempdata = indata{i};
    tempdata(any(isnan(tempdata),2),:)=[];
    indata{i} = tempdata;
end
clear i option

indata = cell2mat(indata);
numsubs = length(indata)/numcond;

option.method = 1;
option.num_perm = 500;
option.num_boot = 100;

%young_res = pls_analysis({indata},numsubs,numcond,option);
old_res = pls_analysis({indata},numsubs,numcond,option);
%%

%% Plotting
res = young_res;
LV = 1;
K = 4;
p = res.perm_result.sprob(LV);
long_labs = {'SS Y Pre';'Pop Y Pre';'Nov Y Pre';...
    'SS O Pre';'Pop O Pre';'Nov O Pre';...
    'SS O Post';'Pop O Post';'Nov O Post'};

figure
subplot(1,2,1)
z = res.boot_result.orig_usc(:,LV);
bar(z)
hold on
yneg = res.boot_result.llusc(:,LV);
ypos = res.boot_result.ulusc(:,LV);
errorbar(1:length(z),z,yneg-z,ypos-z,'.')
colorbar off
xticklabels(long_labs)
xtickangle(45)
grid on
title(sprintf('PLS Stimuli: Age and Stimuli, p = %f',p),'FontSize',12)

x = res.boot_result.compare_u(:,LV);
x = reshape(x,K,[]);
%x(abs(x)<6) = 0;
subplot(1,2,2)
pbaspect([1 1 1])
shower_tile_plot(x);
xticks([])
yticks([])
ylabel('States')
colorbar
title(sprintf('LV %d',LV),'FontSize',12)

clear x y z p yneg ypos LV res ans option numsubs numcond

%% Save stuff into the array
Session_Slices.Pre_FO = Pre_FO;
Session_Slices.Post_FO = Post_FO;
Session_Slices.Pre_FO_young=Pre_FO_young;
Session_Slices.Pre_FO_old=Pre_FO_old;
Session_Slices.Post_FO_old=Post_FO_old;
Session_Slices.Pre_Stim=Pre_Stim; 
Session_Slices.Post_Stim=Post_Stim;
Session_Slices.Pre_Stim_young=Pre_Stim_young;
Session_Slices.Pre_Stim_old=Pre_Stim_old;
Session_Slices.Post_Stim_old=Post_Stim_old;

fMRI_hmm4_music.PrePost_PLS.input = Age_Task_Sess;
fMRI_hmm4_music.PrePost_PLS.results.general = session_res;
fMRI_hmm4_music.PrePost_PLS.results.Pre = presess_res;
fMRI_hmm4_music.PrePost_PLS.results.Old = oldsess_res;
fMRI_hmm4_music.Session_Slices = Session_Slices;
%% Plot some stuff!
figure
Age_stim = Age_Task_Sess(:,1:2);
stim_labs = {'Self-Select';'Popular';'Novel'};
for i = 1:3
    subplot(2,3,i)
    bar(nanmean(Age_stim{i,1}))
    grid on
    ylim([0 0.4])
    %imagesc(Age_stim{i,1})
    %caxis([0 50])
    ylabel('Fractional Occupancy')
    xticks(1:K)
    xlabel('State')
    %colorbar
    title(strcat(stim_labs{i},' Young'),'FontSize',12)

    subplot(2,3,i+3)
    bar(nanmean(Age_stim{i,2}))
    grid on
    %imagesc(Age_stim{i,2})
    %caxis([0 50])
    ylabel('Fractional Occupancy')
    xticks(1:K)
    xlabel('State')
    %colorbar
    title(strcat(stim_labs{i},' Old'),'FontSize',12)
end
clear i Age_stim stim_labs
%% Match Group
match_pre = [1;3;5;8;10;9;11;14;15;21;26;27;29;33;37];
match_post = [2;4;7;12;13;16;20;23;24;25;28;30;32;35;38];
tempdata = fMRI_hmm4_music.Metrics.FO;
indata = {[tempdata(match_pre,:);tempdata(match_post,:)]};

clear option

option.method = 1;
option.num_perm = 500;
option.num_boot = 100;

match_FO = pls_analysis(indata,15,2,option);
%% Slice it up
Pre_MBI_Stim = fMRI_hmm4_music.stim_FO(match_pre,:);
Post_MBI_Stim = fMRI_hmm4_music.stim_FO(match_post,:);
%% Average it down for PLS:
Big_Mat = {Pre_MBI_Stim;Post_MBI_Stim};
Match_Sess = cell(3,2);
K = 4;
for mat = 1:length(Big_Mat)
    tempdata = Big_Mat{mat};
    for cat = 1:size(tempdata,2) %category (ss,popular, novel)
        submat = tempdata(:,cat);% get the relevant row
        outmat = nan(length(tempdata),K);%pre-allocate the output matrix
        for part = 1:length(outmat) %rows, or participants
            tempmat = submat{part}; %get the matrix
            if isempty(tempmat)
                continue
            else
                outmat(part,:) = nanmean(tempmat);%take the average
            end
        end
        Match_Sess{cat,mat} = outmat;%save the average data
    end
end
clear Big_Mat K mat tempdata cat submat outmat tempmat part 
%% And now some PLSes!
addpath(genpath('Pls'))
indata = {cell2mat(Match_Sess(:,1));cell2mat(Match_Sess(:,2))};
numcond = 3;
numsubs = zeros(1,length(indata));

for i = 1:length(indata)
    tempdata = indata{i};
    tempdata(any(isnan(tempdata),2),:)=[];
    indata{i} = tempdata;
    numsubs(i) = length(tempdata)/numcond;
end
clear i option

option.method = 1;
option.num_perm = 500;
option.num_boot = 100;

matchsess_res = pls_analysis({cell2mat(indata(1:2))},15,numcond*2,option);
%% Plotting
res = matchsess_res;
LV = 1;
K = 4;
p = res.perm_result.sprob(LV);
long_labs = {'SS O Pre';'Pop O Pre';'Nov O Pre';...
    'SS O Post';'Pop O Post';'Nov O Post'};

figure
subplot(1,2,1)
z = res.boot_result.orig_usc(:,LV);
bar(z)
hold on
yneg = res.boot_result.llusc(:,LV);
ypos = res.boot_result.ulusc(:,LV);
errorbar(1:length(z),z,yneg-z,ypos-z,'.')
colorbar off
xticklabels(long_labs)
xtickangle(45)
grid on
title(sprintf('PLS Pre- & Post-MBI, p = %f',p),'FontSize',12)

x = res.boot_result.compare_u(:,LV);
x = reshape(x,K,[]);
%x(abs(x)<6) = 0;
subplot(1,2,2)
pbaspect([1 1 1])
shower_tile_plot(x);
xticks([])
yticks([])
ylabel('States')
colorbar
title(sprintf('LV %d',LV),'FontSize',12)

clear x y z p yneg ypos LV res ans option numsubs numcond
%% Plot some stuff!
figure
Age_stim = Match_Sess(:,1:2);
stim_labs = {'Self-Select';'Popular';'Novel'};
for i = 1:3
    subplot(2,3,i)
    imagesc(Age_stim{i,1})
    caxis([.1 .8])
    ylabel('Participants')
    xticks(1:5)
    xlabel('State')
    colorbar
    title(strcat('FO ',stim_labs{i},' Pre-MBI'),'FontSize',12)

    subplot(2,3,i+3)
    imagesc(Age_stim{i,2})
    caxis([.1 .8])
    ylabel('Participants')
    xticks(1:5)
    xlabel('State')
    colorbar
    title(strcat('FO ',stim_labs{i},' Post-MBI'),'FontSize',12)
end
clear i Age_stim stim_labs
%%
fMRI_hmm4_music.PrePost_PLS.input_match = Match_Sess;
fMRI_hmm4_music.PrePost_PLS.results.match = matchsess_res;