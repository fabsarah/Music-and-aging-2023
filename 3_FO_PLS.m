%% Now, let's start the analysis!
% This picks up where step 1 (HMM_Estimation) left off. Most of the data here is from the 
% structure array created in that step. If you're using this from "scratch", I've tried to 
% give descriptions throughout of what the matrices are/their dimensionality
% 
%
% First up, we're going to use the block_dex variable (tells us which block is
% which for stimuli, ratings, etc) to mke some new matrices.

tempages = Age; % this is a matrix with participants' ages 
fMRI_hmm4_music.Ages = tempages; %save it to the array

%%
tempdata = presentation_table; % this is a table of stimulus triggers for all participants
block_dex = cell(80,1);% 80 participants
start = 1;% initialize a counter, stimulus-wise
stop = start+23;% 24 stimuli for each
for part = 1:80
    tempdex = tempdata(start:stop,:);% go through the data 
    block_dex{part} = tempdex;% get the values, add them to the new array
    start = stop+1;% get the new start position
    stop = start+23;% and the new stop position
end
clear start stop part tempdex
fMRI_hmm4_music.block_dex = block_dex; %add it to the array

%% And now use this index to organize the FO matrices:
tempdata = fMRI_hmm4_music;
stim_FO = cell(80,3);% 80 participants, three stimulus types
for part = 1:80
    index = tempdata.block_dex{part};% get the participant's stimulus trigger order
    block = index(:,2);% for me, this information was stored in the second column
    tempFO = tempdata.block_FO{part};% get the participant's FO data
    for stim = 1:3 % 3 stimulus types
        newmat = tempFO(block==stim,:);% get the relevant FO
        stim_FO{part,stim} = newmat;%put it into the new matrix
    end
end
fMRI_hmm4_music.stim_FO = stim_FO;% save it to the array!
clear part index block tempFO stim newmat %make look nice

%% Now, average and assemble into three more matrices
% stim_FO gives us a participant*stimulus category cell array, so now let's do some 
% averaging so we can run PLS:

K = fMRI_hmm4_music.options.K;
tempdata = fMRI_hmm4_music.stim_FO;
stim_FO_all = cell(3,1);
for stimcat = 1:3 % stimulus category
    stimmat = nan(80,K);
    for part = 1:80
        partv = tempdata{part,stimcat}; %get the data
        if isempty(partv)
            continue
        else
                tempFO = nanmean(partv);% average it down
                stimmat(part,:) = tempFO; %put it in the temporary matrix
        end
    end
    stim_FO_all{stimcat} = stimmat;% add the matrix to the new cell array
end
fMRI_hmm4_music.stim_FO_all = stim_FO_all;% save it to the array
clear stimcat part stimmat part partv state tempFO % make look nice

%% PLS time!
addpath(genpath('Pls'))
clear option % make sure there isn't an old option file kicking around
indata_stim = cell2mat(stim_FO_all);% set our input data. We should have a task*1 cell array where each cell is a participant*K matrix of FO values
indata_stim(isnan(indata_stim))=0;% make sure there are no nans
option.method = 1;% mean-centred PLS
option.num_perm = 500;
option.num_boot = 100;
nparts = 80;% 80 participants
ntask = 3;% 3 tasks/stimulus categories

stim_res = pls_analysis({indata_stim},nparts,ntask,option); % run the analysis! 

%% Plot the LVs!
res = stim_res;
LV = 1;
p = res.perm_result.sprob(LV);% can check this whole array to see how many significant LVs you have. It works like a traditional p value
K = fMRI_hmm4_music.options.K;

figure
subplot(1,2,1)
z = res.boot_result.orig_usc(:,LV);% get the usc matrix. This shows us the "contrast" between tasks/stimuli types
bar(z)
hold on
yneg = res.boot_result.llusc(:,LV);
ypos = res.boot_result.ulusc(:,LV);
errorbar(1:length(z),z,yneg-z,ypos-z,'.')
xticklabels({'Self-Select';'Popular';'Novel'})
xtickangle(45)
grid on
title(sprintf('PLS Stimuli: All Ages, p = %d',p),'FontSize',12)

x = res.boot_result.compare_u(:,LV);
x = reshape(x,K,[]);
%x(abs(x)<6) = 0;% can threshold this for readability 
subplot(1,2,2)
pbaspect([1 1 1])
shower_tile_plot(x);%imagesc will also work, just not as "crisp"
xticks([])
yticks([])
ylabel('States')
colorbar
title(sprintf('LV %d',LV),'FontSize',12)

clear res K LV p x y* z
%% Now...is there an age effect?
% Here, we'll divide up our matrix into age groups:
Age_stim = cell(3,2);% pre-allocate a new cell
idx_old = fMRI_hmm4_music.Ages>40;% sorry :(
idx_young = fMRI_hmm4_music.Ages<40;
for i = 1:3
    tempmat = fMRI_hmm4_music.stim_FO_all{i};% get the FO data
    Age_stim{i,1} = tempmat(idx_young==1,:);separate into separate cells
    Age_stim{i,2} = tempmat(idx_old==1,:);
end
clear i
%%
%addpath(genpath('Pls'));
clear option
indata{1,1} = cell2mat(Age_stim(:,1));
indata{2,1} = cell2mat(Age_stim(:,2));
for i = 1:length(indata) % make sure there are no nans
    tempdata = indata{i,:};
    tempdata(isnan(tempdata)) = 0;
    indata{i,:} = tempdata;
end

option.method = 1;% mean-centred PLS
option.num_perm = 500;
option.num_boot = 100;
nparts = [length(indata{1});length(indata{2})];% we now have different group lengths
ntask = 3;

all_res = pls_analysis(indata,{nparts},ntask,option);
%young_res = pls_analysis({indata{1}},nparts(1),ntask,option);
%old_res = pls_analysis({indata{2}},nparts(2),ntask,option);

%% Go up a block for the syntax to plot young and old separately. Stay here to plot everyone together:
res = all_res;
LV = 1;
p = res.perm_result.sprob(LV);
K = fMRI_hmm4_music.options.K;

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
clear z y* x K res LV p

%% Save the results!
fMRI_hmm4_music.FOPLS.input = Age_stim;
fMRI_hmm4_music.FOPLS.results.general = stim_res;
fMRI_hmm4_music.FOPLS.results.Young = young_res;
fMRI_hmm4_music.FOPLS.results.Old = old_res;
fMRI_hmm4_music.FOPLS.results.all = all_res;

%% Plot the input data all together:
K = fMRI_hmm4_music.options.K;
ntask = 3;

figure
stim_labs = {'Self-Select';'Popular';'Novel'};
for i = 1:3
    subplot(1,ntask,i)
    imagesc(stim_FO_all{i})
    %caxis([0 50])
    ylabel('Participants')
    xticks(1:K)
    xlabel('State')
    colorbar
    title(stim_labs{i},'FontSize',12)
end
clear i

%% Plot the input data by age:
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
    title(strcat(stim_labs{i},' Younger Adults'),'FontSize',12)

    subplot(1,3,i+3)
    imagesc(Age_stim{i,2})
    %caxis([0 50])
    ylabel('Participants')
    xticks(1:5)
    xlabel('State')
    colorbar
    title(strcat(stim_labs{i},' Older Adults'),'FontSize',12)
end
clear i
