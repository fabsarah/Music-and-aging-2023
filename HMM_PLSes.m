%% Some initial checks with K = 5
tempages = HMMfMRI.sub_info.Age;
tempages([9,26,30,34,69,86]) = []; %remove the noisy folks
fMRI_hmm4_music.Ages = tempages;
%fMRI_hmm5_music.Ages = tempages;
%fMRI_hmm7_new.Ages = tempages;
%fMRI_hmm8_new.Ages = tempages;
%%
tempdata = musbidtidytable(1:(24*86),8:11);
block_dex = cell(86,1);
start = 1;
stop = start+23;
for part = 1:86
    tempdex = tempdata(start:stop,:);
    block_dex{part} = tempdex;
    start = stop+1;
    stop = start+23;
end
clear start stop part tempdex
block_dex([9,26,30,34,69,86]) = [];
fMRI_hmm4_music.block_dex = block_dex;
%musbidtidytable([9,26,30,34,69,86]) = [];
%% Look at the input data
%indata_rest = 
%indata_music = fMRI_hmm5.Metrics.FO(81:end,:);
plotdata = fMRI_hmm4_music;
K = plotdata.options.K;
my_map = parula(K);
state_labs = cell(K,1);
for i = 1:K
    tempstr = sprintf('State %d',i);
    state_labs{i} = tempstr;
end
clear i tempstr

figure
subplot(1,2,1)
imagesc(plotdata.Metrics.FO(1:80,:))
%colormap(my_map)
title('FO Rest','FontSize',12)
colorbar
yticks([])
ylabel('Participants')
xticks(1:5)
xlabel('State')

subplot(1,2,2)
imagesc(plotdata.Metrics.FO(81:end,:))
%colormap(my_map)
title('FO Music','FontSize',12)
colorbar
yticks([])
ylabel('Participants')
xticks(1:5)
xlabel('State')
%%
% FO PLS
addpath(genpath('Pls'))
ages = HMMfMRI.sub_info.Age;
indata_rest = fMRI_hmm4_music.Metrics.FO(1:80,:);
indata_music = fMRI_hmm4_music.Metrics.FO(81:end,:);
indata_all = fMRI_hmm4_music.Metrics.FO;
%%
clear option
option.method = 3;
option.num_perm = 500;
option.num_boot = 100;
%option.stacked_behavdata = fMRI_hmm5_new.Ages;
option.stacked_behavdata = [fMRI_hmm4_music.Ages;fMRI_hmm4_music.Ages];

%rest_res = pls_analysis({indata_rest},80,1,option);
%music_res = pls_analysis({indata_music},80,1,option);
rvm_bres = pls_analysis({indata_all},80,2,option);
%rvm_res = pls_analysis({[indata_rest;indata_music]},80,2,option);
%% Plot it!
res = rest_res;
LV = 1;
p = res.perm_result.sprob(LV);

figure
x = res.boot_result.compare_u(:,LV);
x = reshape(x,5,[]);
subplot(1,2,1)
pbaspect([1 1 1])
imagesc(x)
colorbar
title(sprintf('Age-variable PLS, Rest, p = %d',p),'FontSize',16)
clear z x

p = music_res.perm_result.sprob(LV);
x = music_res.boot_result.compare_u(:,LV);
x = reshape(x,5,[]);
subplot(1,2,2)
pbaspect([1 1 1])
imagesc(x)
colorbar
title(sprintf('Age-variable PLS, Music, p = %d',p),'FontSize',16)

clear z x
%% Plot with Bars
res = rvm_res;
LV = 1;
p = res.perm_result.sprob(LV);

figure
subplot(1,2,1)
z = res.boot_result.orig_usc(:,LV);
bar(z)
hold on
yneg = res.boot_result.llusc(:,LV);
ypos = res.boot_result.ulusc(:,LV);
errorbar(1:length(z(:,LV)),z(:,LV),yneg-z(:,LV),ypos-z(:,LV),'.')
colorbar off
xticklabels({'Rest';'Music'})
xtickangle(45)
grid on
title(sprintf('Mean-centred PLS (FO), Rest vs Music, p = %d',p),'FontSize',12)

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
%%
idx_old = fMRI_hmm4_music.Ages>40;
idx_young = fMRI_hmm4_music.Ages<40;
FO_yr = indata_rest(idx_young==1,:);
FO_or = indata_rest(idx_old==1,:);
FO_ym = indata_music(idx_young==1,:);
FO_om = indata_music(idx_old==1,:);
%%
clear option
indata_ab = [FO_yr(1:39,:);FO_or;FO_ym(1:39,:);FO_om];
option.method = 1;
option.num_perm = 500;
option.num_boot = 100;

ab_res = pls_analysis({indata_ab},39,4,option);
%% Plot it!
res = ab_res;
LV = 1;
p = res.perm_result.sprob(LV);

%% Plot with Bars
res = ab_res;
LV = 3;
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
xticklabels({'Rest Y';'Rest O';'Music Y';'Music O'})
xtickangle(45)
grid on
title(sprintf('Beh. PLS (FO), Rest vs Music, p = %d',p),'FontSize',12)

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
%% Behavioural PLS LV syntax:
%% Age-variable behavioural PLS: rest vs music plot:
res = rvm_bres;
LV = 2;
p = res.perm_result.sprob(LV);

figure
subplot(1,2,1)
z = res.boot_result.orig_corr(:,LV);
bar(z)
hold on
yneg = res.boot_result.llcorr(:,LV);
ypos = res.boot_result.ulcorr(:,LV);
errorbar(1:length(z),z,yneg-z,ypos-z,'.')
colorbar off
xticklabels({'Rest';'Music';'Music Y';'Music O'})
xtickangle(45)
grid on
title(sprintf('Behavioural PLS (FO), Rest vs Music, p = %d',p),'FontSize',12)

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
%% Figure out the plotting!
State_means = getMean(fMRI_hmm4_music.hmm);
thresh_means = nan(size(State_means));

for i = 1:5
    tempdata = State_means(:,i);
    index = prctile(abs(tempdata),[75:99],1);
    minvalue = min(min(index));
    tempdata(abs(tempdata)<minvalue) = nan;
    thresh_means(:,i) = tempdata;
end
clear i tempdata minvalue index
fMRI_hmm4_music.Means.thresh = thresh_means;
fMRI_hmm4_music.Means.raw = State_means;

%% Look at it all
figure
subplot(1,2,2)
imagesc(thresh_means)
colorbar
subplot(1,2,1)
imagesc(State_means)
colorbar
%% Save it!

fMRI_hmm4_music.PLS.FO.results.rest_res = rest_res;
fMRI_hmm4_music.PLS.FO.results.music_res = music_res;
fMRI_hmm4_music.PLS.FO.results.rvm_res = rvm_res;
fMRI_hmm4_music.PLS.FO.results.ab_res = ab_res;
fMRI_hmm4_music.PLS.FO.results.rvm_bres = rvm_bres;


