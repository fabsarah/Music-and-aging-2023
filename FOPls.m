%% Now, let's get weird!
%here, we're going to use the block_dex table (tells us which block is
%which for stimuli, ratings, etc) to mke some new matrices
tempdata = fMRI_hmm4_music;
stim_FO = cell(80,3);
for part = 1:80
    index = tempdata.block_dex{part};
    block = table2array(index(:,2));
    block(block==4) = 3;%combining novel and BP
    tempFO = tempdata.block_FO{part};
    for stim = 1:3
        newmat = tempFO(block==stim,:);
        stim_FO{part,stim} = newmat;
    end
end
fMRI_hmm4_music.stim_FO = stim_FO;
clear part index block tempFO stim newmat
% NB cell 65 has no data (NANs in the input)
%% Now, average and assemble into three more matrices
% still needs to be checked for indexing errors!
K = 4;
tempdata = fMRI_hmm4_music.stim_FO;
stim_FO_all = cell(4,1);
for stimcat = 1:3
    stimmat = nan(80,K);
    for part = 1:80
        partv = tempdata{part,stimcat};
        if isempty(partv)
            continue
        else
                tempFO = nanmean(partv);
                stimmat(part,:) = tempFO;
        end
    end
    stim_FO_all{stimcat} = stimmat;
end
fMRI_hmm4_music.stim_FO_all = stim_FO_all;
clear stimcat part stimmat part partv state tempFO
%% Here goes nothing!
addpath(genpath('Pls'))
clear option
indata_stim = cell2mat(stim_FO_all);
indata_stim(isnan(indata_stim))=0;
option.method = 1;
option.num_perm = 500;
option.num_boot = 100;

stim_res = pls_analysis({indata_stim},80,3,option);
%% Plot with Bars
% There is one significant LV here! Self-selected vs all other musics.
% Amazing!
res = stim_res;
LV = 1;
p = res.perm_result.sprob(LV);

figure
subplot(1,2,1)
z = res.boot_result.orig_usc(:,LV);
bar(z)
hold on
yneg = res.boot_result.llusc(:,LV);
ypos = res.boot_result.ulusc(:,LV);
errorbar(1:length(z),z,yneg-z,ypos-z,'.')
colorbar off
xticklabels({'Self-Select';'Popular';'Novel'})
xtickangle(45)
grid on
title(sprintf('PLS Stimuli: All Ages, p = %d',p),'FontSize',12)

x = res.boot_result.compare_u(:,LV);
x = reshape(x,4,[]);
%x(abs(x)<6) = 0;
subplot(1,2,2)
pbaspect([1 1 1])
shower_tile_plot(x)
xticks([])
yticks([])
ylabel('States')
colorbar
title(sprintf('LV %d',LV),'FontSize',12)
%% Now...is there an age effect?
Age_stim = cell(3,2);
idx_old = fMRI_hmm4_music.Ages>40;
idx_young = fMRI_hmm4_music.Ages<40;
for i = 1:3
    tempmat = fMRI_hmm4_music.stim_FO_all{i};
    Age_stim{i,1} = tempmat(idx_young==1,:);
    Age_stim{i,2} = tempmat(idx_old==1,:);
end
clear i
%%
%addpath(genpath('Pls'));
clear option
indata{1,1} = cell2mat(Age_stim(:,1));
indata{2,1} = cell2mat(Age_stim(:,2));
for i = 1:length(indata)
    tempdata = indata{i,:};
    tempdata(isnan(tempdata)) = 0;
    indata{i,:} = tempdata;
end
%indata(isnan(indata))=0;
option.method = 1;
option.num_perm = 500;
option.num_boot = 100;

%all_res = pls_analysis(indata,[41;39],4,option);
%young_res = pls_analysis({indata{1}},41,4,option);
old_res = pls_analysis({indata{2}},39,4,option);
%% testing PLS
indata{1,1} = Age_stim{1,1};
indata{2,1} = Age_stim{1,2};

for i = 1:length(indata)
    tempdata = indata{i,:};
    tempdata(isnan(tempdata)) = 0;
    indata{i,:} = tempdata;
end
clear i tempdata option
num_subjs = [41;39];
option.method = 1;
option.num_perm = 500;
option.num_boot = 100;

ss_res = pls_analysis(indata,num_subjs,1,option);


%% Go up a block to plot young and old. Stay here to plot everyone together
res = ss_res;
LV = 1;
p = res.perm_result.sprob(LV);

figure
subplot(1,2,1)
z = res.boot_result.orig_usc(:,LV);
bar(z)
hold on
yneg = res.boot_result.llusc(:,LV);
ypos = res.boot_result.ulusc(:,LV);
errorbar(1:length(z),z,yneg-z,ypos-z,'.')
colorbar off
xticklabels({'Self-Select Y';'Popular Y';'Novel Y';...
    'Self-Select O';'Popular O';'Novel O'})
xtickangle(45)
grid on
title(sprintf('PLS Stimuli: Age and Stimuli, p = %d',p),'FontSize',12)

x = res.boot_result.compare_u(:,LV);
x = reshape(x,5,[]);
%x(abs(x)<6) = 0;
subplot(1,2,2)
pbaspect([1 1 1])
shower_tile_plot(x)
xticks([])
yticks([])
ylabel('States')
colorbar
title(sprintf('LV %d',LV),'FontSize',12)
%% Save the variables!
fMRI_hmm4_music.FOPLS.input = Age_stim;
fMRI_hmm4_music.FOPLS.results.general = stim_res;
fMRI_hmm4_music.FOPLS.results.Young = young_res;
fMRI_hmm4_music.FOPLS.results.Old = old_res;
fMRI_hmm4_music.FOPLS.results.all = all_res;
fMRI_hmm4_music.FOPLS.results.ss = ss_res;
%%
figure
stim_labs = {'Self-Select';'Popular';'Novel'};
for i = 1:3
    subplot(1,3,i)
    imagesc(stim_FO_all{i})
    %caxis([0 50])
    ylabel('Participants')
    xticks(1:5)
    xlabel('State')
    colorbar
    title(stim_labs{i},'FontSize',12)
end
clear i
%%
figure
stim_labs = {'Self-Select';'Popular';'Novel'};
for i = 1:3
    subplot(1,3,i)
    imagesc(Age_stim{i,1})
    %caxis([0 50])
    ylabel('Participants')
    xticks(1:5)
    xlabel('State')
    colorbar
    title(strcat(stim_labs{i},' Young'),'FontSize',12)

    subplot(1,3,i+3)
    imagesc(Age_stim{i,2})
    %caxis([0 50])
    ylabel('Participants')
    xticks(1:5)
    xlabel('State')
    colorbar
    title(strcat(stim_labs{i},' Old'),'FontSize',12)
end
clear i