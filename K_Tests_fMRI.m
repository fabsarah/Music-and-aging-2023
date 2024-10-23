%% Load everything
load('fMRI_hmm3_music.mat')
load('fMRI_hmm4_music.mat')
load('fMRI_hmm5_music.mat')
load('fMRI_hmm6_music.mat')
load('fMRI_hmm7_music.mat')
load('fMRI_hmm8_music.mat')
load('fMRI_hmm9_music.mat')
load('fMRI_hmm10_music.mat')
load('fMRI_hmm12_music.mat')
load('fMRI_hmm15_music.mat')
load('fMRI_hmm20_music.mat')
load('HMMfMRI.mat')
%% Step 1: get our giant timeseries vector
music_data = HMMfMRI.ts.ts_music;
music_data([9,26,30,34,69,86]) = [];
Indata = music_data;
for i = 1:length(Indata)
    Indata{i} = Indata{i}';
end
clear i music_data
%% Free Energy Recalculate
addpath(genpath('HMM-MAR-master'))
dex = 1;
fepc = nan(11,1);
llc = nan(11,1);
%% Do this manually. Like a fool. 
data = fMRI_hmm20_music;
[fe,ll] = hmmfe(Indata,data.T,data.hmm,data.Gamma,data.Xi);
fepc(dex) = fe;
llc(dex) = mean(ll);
dex = dex+1;
clear data fe ll
%%
K_tests.fe = fepc;
K_tests.llc = llc;

%% Plot!
klabs = {'K3','K4','K5','K6','K7','K8','K9','K10','K12','K15','K20'};
figure
%subplot(1,2,1)
bar(normalize(fepc))
grid on
xticklabels(klabs)
title('Free Energy K Comparison')
clear dex llc
%% Assemble everything into a new object
K_tests.FO{1} = fMRI_hmm3_music.Metrics.FO;
K_tests.FO{2} = fMRI_hmm4_music.Metrics.FO;
K_tests.FO{3} = fMRI_hmm5_music.Metrics.FO;
K_tests.FO{4} = fMRI_hmm6_music.Metrics.FO;
K_tests.FO{5} = fMRI_hmm7_music.Metrics.FO;
K_tests.FO{6} = fMRI_hmm8_music.Metrics.FO;
K_tests.FO{7} = fMRI_hmm9_music.Metrics.FO;
K_tests.FO{8} = fMRI_hmm10_music.Metrics.FO;
K_tests.FO{9} = fMRI_hmm12_music.Metrics.FO;
K_tests.FO{10} = fMRI_hmm15_music.Metrics.FO;
K_tests.FO{11} = fMRI_hmm20_music.Metrics.FO;
K_tests.labels = {'K3';'K4';'K5';'K6';'K7';'K8';'K9';'K10';'K12';'K15';'K20'};
%% Average them
% for i = 1:11
%     tempdata = K_tests.FO{i};
%     tempdata = cat(3,tempdata{:});
%     tempdata = nanmean(tempdata,3);
%     K_tests.FO{i} = tempdata;
% end
%% Run FO PLSes
addpath(genpath('Pls'))
%% Get the age indices
idx_pre = fMRI_hmm4_music.pre_dex==1;
idx_post = fMRI_hmm4_music.pre_dex==0;
idx_old = fMRI_hmm4_music.Ages>40;
idx_young = fMRI_hmm4_music.Ages<40;
%% Slice it up
Pre_FO_young = cell(11,1);
Pre_FO_old = cell(11,1);

for i = 1:11
    Pre_FO_young{i} = K_tests.FO{i}(idx_pre==1&idx_young==1);
    Pre_FO_old{i} = K_tests.FO{i}(idx_pre==1&idx_old==1);
end
K_tests.Pre_FO_young = Pre_FO_young;
K_tests.Pre_FO_old = Pre_FO_old;
clear i idx*
%% Appetizer PLS
addpath(genpath('Pls'))
FO_data = K_tests;
FO_results = cell(11,1);
for k = 1:11
    indata = {K_tests.Pre_FO_young{k};K_tests.Pre_FO_old{k}};
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

    K_res = pls_analysis(indata,numsubs,1,option);
    FO_results{k} = K_res;
end
clear i k indata numsubs tempdata option K_res
%% Assemble and plot p values
FO_p = nan(1,11);
for k = 1:11
    temp = FO_results{k}.perm_result.sprob;
    FO_p(k) = temp;
end
clear k temp
%%
figure
bar(FO_p)
grid on
xticklabels(K_tests.labels)
title('PLS p values variable K')
%% Now correlate usc and vsc
LV_corr = nan(1,11);
for k = 1:11
    u = FO_results{k}.usc;
    v = FO_results{k}.vsc;
    tempcorr = corrcoef(u,v);
    LV_corr(k) = tempcorr(2,1);
end
clear k tempcorr u v
%%
figure
bar(LV_corr)
ylabel('r value')
grid on
xticklabels(K_tests.labels)
title('USC VSC Correlations, LV1 Variable K')

%% STATIS tests
FCs = [fMRI_hmm4_music.FC,fMRI_hmm7_music.FC];
FCs = cat(3,FCs{:});
%%
FC_res = distatis2(FCs);
%%
figure
shower_tile_plot(FC_res.C);
colorbar
title('Similarity matrix, K = 4 and K = 7')
% a little difficult to interpret. State 6 in K = 7 not similar to anything
%% Try the dot product
addpath(genpath('HMM-MAR-master'))
res7 = fMRI_hmm7_music;
m7 = getMean(res7.hmm);
res4 = fMRI_hmm4_music;
m4 = getMean(res4.hmm);
clear res7 res4
%%
m7 = normalize(m7,1);
m4 = normalize(m4,1);
%%
figure
subplot(1,2,1)
imagesc(m4)
title('K = 4')
ylabel('Region')
xlabel('State')
colorbar

subplot(1,2,2)
imagesc(m7)
colorbar
title('K = 7')
ylabel('Region')
xlabel('State')
%%
mdot = m7'*m4;
figure
imagesctxt(mdot)
colorbar
title('Dot Product: K = 7, K = 4')
%% TP Matrices: Assemble everything into a new object
K_tests.TP{1} = getTransProbs(fMRI_hmm3_music.hmm);
K_tests.TP{2} = getTransProbs(fMRI_hmm4_music.hmm);
K_tests.TP{3} = getTransProbs(fMRI_hmm5_music.hmm);
K_tests.TP{4} = getTransProbs(fMRI_hmm6_music.hmm);
K_tests.TP{5} = getTransProbs(fMRI_hmm7_music.hmm);
K_tests.TP{6} = getTransProbs(fMRI_hmm8_music.hmm);
K_tests.TP{7} = getTransProbs(fMRI_hmm9_music.hmm);
K_tests.TP{8} = getTransProbs(fMRI_hmm10_music.hmm);
K_tests.TP{9} = getTransProbs(fMRI_hmm12_music.hmm);
K_tests.TP{10} = getTransProbs(fMRI_hmm15_music.hmm);
K_tests.TP{11} = getTransProbs(fMRI_hmm20_music.hmm);
K_tests.labels = {'K3';'K4';'K5';'K6';'K7';'K8';'K9';'K10';'K12';'K15';'K20'};
%%
addpath(genpath('BCT'))
STR_test = cell(10,1);
DEG_test = cell(10,1);
CC_test = cell(10,1);
for part = 1:10
    tempdata = K_tests.TP{part};
    thresh = prctile(abs(tempdata),50,1);
    tempdata(abs(tempdata)<thresh) = 0;
    str = strengths_dir(tempdata);
    [indeg,outdeg] = degrees_dir(tempdata);
    cc = clustering_coef_wd(tempdata);
    STR_test{part} = str;
    DEG_test{part} = [indeg;outdeg];
    CC_test{part} = cc;
end
clear part tempdata thresh str indeg outdeg cc
%%
part = 2;
figure
%subplot(1,2,1)
bar(normalize(STR_test{part}))
grid on
ylabel('Strength')
xlabel('State')
title(sprintf('State Strength, %s',K_tests.labels{part}))

% subplot(1,2,2)
% bar(DEG_test{part})
% grid on
% ylabel('Strength')
% xlabel('State')
% title(sprintf('State Strength, %s',K_tests.labels{part}))

%% Assemble everything into a new object
K_tests.FC{1} = fMRI_hmm3_music.FC;
K_tests.FC{2} = fMRI_hmm4_music.FC;
K_tests.FC{3} = fMRI_hmm5_music.FC;
K_tests.FC{4} = fMRI_hmm6_music.FC;
K_tests.FC{5} = fMRI_hmm7_music.FC;
K_tests.FC{6} = fMRI_hmm8_music.FC;
K_tests.FC{7} = fMRI_hmm9_music.FC;
K_tests.FC{8} = fMRI_hmm10_music.FC;
K_tests.FC{9} = fMRI_hmm12_music.FC;
K_tests.FC{10} = fMRI_hmm15_music.FC;
K_tests.FC{11} = fMRI_hmm20_music.FC;
%%
K_tests.mFC{1} = fMRI_hmm3_music.mean_FC;
K_tests.mFC{2} = fMRI_hmm4_music.mean_FC;
K_tests.mFC{3} = fMRI_hmm5_music.mean_FC;
K_tests.mFC{4} = fMRI_hmm6_music.mean_FC;
K_tests.mFC{5} = fMRI_hmm7_music.mean_FC;
K_tests.mFC{6} = fMRI_hmm8_music.mean_FC;
K_tests.mFC{7} = fMRI_hmm9_music.mean_FC;
K_tests.mFC{8} = fMRI_hmm10_music.mean_FC;
K_tests.mFC{9} = fMRI_hmm12_music.mean_FC;
K_tests.mFC{10} = fMRI_hmm15_music.mean_FC;
K_tests.mFC{11} = fMRI_hmm20_music.mean_FC;

%% STATIS the mean FC
mean_FC = cat(3,K_tests.mFC{:});
mean_FC(isinf(mean_FC)) = 0;
%%
addpath(genpath('Pls'))
res_mFC = distatis2(mean_FC);
K_tests.res_MFC = res_mFC;
%% Plot the similarity matrices
figure
shower_tile_plot(res_mFC.C);
caxis([0.5 1])
colorbar
pbaspect([1 1 1])
xticks(1.5:1:11.5)
xticklabels(K_tests.labels)
yticks(1.5:1:11.5)
yticklabels(flipud(K_tests.labels))
title('Similarity Matrix: mean FC for all Ks')
%% Get the dissimilarity matrix
DSM = 1-res_mFC.C;
%% Big ass STATIS
start = 1;
Kdex = [3,4,5,6,7,8,9,10,12,15,20];
Big_FC = nan(220,220,sum(Kdex));

for i = 1:length(Kdex)
    tempdata = K_tests.FC{i};
    tempdata = cat(3,tempdata{:});
    tempdata(isinf(tempdata)) = 0;
    tempdata(isnan(tempdata)) = 0;
    Big_FC(:,:,start:start+Kdex(i)-1) = tempdata;
    start = start+Kdex(i);
end
clear i tempdata start
%%
res_BigFC = distatis2(Big_FC);

figure
imagesc(res_BigFC.C);
yticks([])
xticks([])
colorbar

%% Get the dissimilarity matrix
DSM_BFC = 1-res_BigFC.C;
DSM_BFC(isinf(DSM_BFC)) = 0;
K_index = nan(sum(Kdex),1);
start = 1;
for i = 1:length(Kdex)
    K_index(start:start+(Kdex(i)-1),1) = Kdex(i);
    start = start+Kdex(i);
end
%clear i start

%%
DI = dunns(10,DSM_BFC,K_index);
%% Correlating the right things this time:
res7 = fMRI_hmm7_music.PrePost_PLS.results.Pre;
res4 = fMRI_hmm4_music.PrePost_PLS.results.Pre;

corr7 = corrcoef(res7.usc,res7.vsc);
corr4 = corrcoef(res4.usc,res4.vsc);