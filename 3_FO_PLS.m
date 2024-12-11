%% Now, let's start the analysis!
% This picks up where step 2 (testing the Ks) left off, and most of the data here is from the 
% structure array created in step 1 (HMM Estimation) using our "winning" K value. 
% 
% We could use the PLS results from our K tests, but I've shown the non-looped PLS
% here for posterity. Delicious, delicious posterity!
%
% If you're using this from "scratch", I've tried to 
% give descriptions throughout of what the matrices are/their dimensionality
%
% First up, we're going to use the block_dex variable (tells us which block is
% which for stimuli, ratings, etc) to mke some new matrices.

%% Compare folks who have clean pre- and post-intervention data
% these are indices matching up the participants' data pre- and post-intervention
match_pre = [1;3;5;8;10;9;11;14;15;21;26;27;29;33;37];
match_post = [2;4;7;12;13;16;20;23;24;25;28;30;32;35;38];

tempdata = fMRI_hmm4;
music_pre = tempdata.music_FO(match_pre,:);
music_post = tempdata.music_FO(match_post,:);

%% PLS Time!
addpath(genpath('Pls'))
clear option % make sure there isn't an old options file hanging around


indata = {[match_pre;match_post]};

option.method = 1;
option.num_perm = 500;
option.num_boot = 100;
nparts = length(match_pre);% how many participants in the analysis?
nconds = 2;%2 conditions, pre- and post-intervention

match_res = pls_analysis(indata,nparts,nconds,option);% run the analysis

%% Plotting
res = match_res;
K = 4;% the K value for the estimation
LV = 1;
p = res.perm_result.sprob(LV);
short_labs = {'Pre';'Post'};
heading = sprintf('PLS Pre-Post: Music, p = %d',p);

figure
subplot(1,2,1)
z = res.boot_result.orig_usc(:,LV);
bar(z)
hold on
yneg = res.boot_result.llusc(:,LV);
ypos = res.boot_result.ulusc(:,LV);
errorbar(1:length(z),z,yneg-z,ypos-z,'.')
colorbar off
xticklabels(short_labs)
xtickangle(45)
grid on
title(heading,'FontSize',12)

x = res.boot_result.compare_u(:,LV);
x = reshape(x,K,[]);
%x(abs(x)<3) = 0;% the result can be thresholded for readability
subplot(1,2,2)
pbaspect([1 1 1])
shower_tile_plot(x);
xticks([])
yticks([])
ylabel('States')
colorbar
title(sprintf('LV %d',LV),'FontSize',12)

clear x y z p yneg ypos LV res ans option numsubs numcond heading

%% Now the correlation analyses:
% For the next part of the analysis, I'm using correlation matrices I extracted
% between behavioural liking and familiarity data, and the HMM measures

tempdata = fMRI_hmm4_music.LF_PLS;
like_pre = tempdata.Liking.like_mat(match_pre,:);
like_post = tempdata.Liking.like_mat(match_post,:);
fam_pre = tempdata.Fam.fam_mat(match_pre,:);
fam_post = tempdata.Fam.fam_mat(match_post,:);

%% PLS Time!
%addpath(genpath('Pls'))
% Here, I want to look for intervention differences in liking and familiarity
% correlations

clear option

indata_like = {[like_pre;like_post]};
indata_fam = {[fam_pre;fam_post]};
indata_all = {[like_pre;fam_pre;like_post;fam_post]};

option.method = 1;
option.num_perm = 500;
option.num_boot = 100;
nparts = 15;% number of participants
nconds = 2;% number of conditions

match_like = pls_analysis(indata_like,nparts,nconds,option);
match_fam = pls_analysis(indata_fam,nparts,nconds,option);
match_likefam = pls_analysis(indata_all,nparts,nconds*2,option);%nconds here is 4 because we're looking at 2 metrics*2 conditions
clear option indata*

%% Plotting
res = match_likefam;
LV = 1;
K = 4;
p = res.perm_result.sprob(LV);
short_labs = {'Pre';'Post'};
mvr_labs = {'Like Pre';'Fam Pre';'Like Post';'Fam Post'};
heading = sprintf('PLS Pre-Post: Like & Fam, p = %d',p);

figure
subplot(1,2,1)
z = res.boot_result.orig_usc(:,LV);
bar(z)
hold on
yneg = res.boot_result.llusc(:,LV);
ypos = res.boot_result.ulusc(:,LV);
errorbar(1:length(z),z,yneg-z,ypos-z,'.')
colorbar off
xticklabels(mvr_labs)
xtickangle(45)
grid on
title(heading,'FontSize',12)

x = res.boot_result.compare_u(:,LV);
x = reshape(x,4,[]);
%x(abs(x)<3) = 0;
subplot(1,2,2)
pbaspect([1 1 1])
shower_tile_plot(x);
xticks([])
yticks([])
ylabel('States')
colorbar
title(sprintf('LV %d',LV),'FontSize',12)

clear x y z p yneg ypos LV res ans option numsubs numcond heading

%% Plot the indata:
figure
subplot(2,2,1)
imagesc(like_pre)
colorbar
title('Liking Pre','FontSize',12)

subplot(2,2,2)
imagesc(like_post)
colorbar
title('Liking Post','FontSize',12)

subplot(2,2,3)
imagesc(fam_pre)
colorbar
title('Familiarity Pre','FontSize',12)

subplot(2,2,4)
imagesc(fam_post)
colorbar
title('Familiarity Post','FontSize',12)

%% Save it all:
match_PLS.FO.Input.rest_pre = rest_pre;
match_PLS.FO.Input.rest_post = rest_post;
match_PLS.FO.Input.music_pre = music_pre;
match_PLS.FO.Input.music_post = music_post;
match_PLS.FO.Results.match_rest_res = match_rest_res;
match_PLS.FO.Results.match_music_res = match_music_res;
match_PLS.FO.Results.match_mvr_res = match_mvr_res;
match_PLS.FO.Results.match_stim_res = match_stim_res;
match_PLS.FO.Results.match_rest_stim_res =match_rest_stim_res;

match_PLS.corr.Input.like_pre = like_pre;
match_PLS.corr.Input.like_post = like_post;
match_PLS.corr.Input.fam_pre = fam_pre;
match_PLS.corr.Input.fam_post = fam_post;
match_PLS.corr.Results.match_like_res = match_like;
match_PLS.corr.Results.match_fam_res = match_fam;
match_PLS.corr.Results.match_likefam_res = match_likefam;

fMRI_hmm4_music.match_PLS.match_dex.match_pre = match_pre;
fMRI_hmm4_music.match_PLS.match_dex.match_post = match_post;
fMRI_hmm4_music.match_PLS = match_PLS;
