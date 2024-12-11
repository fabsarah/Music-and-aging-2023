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

%% Assemble the fe and ll values into a new object
fepc(1) = fMRI_hmm3_music.fe;
fepc(2) = fMRI_hmm4_music.fe;
fepc(3) = fMRI_hmm5_music.fe;
fepc(4) = fMRI_hmm6_music.fe;
fepc(5) = fMRI_hmm7_music.fe;
fepc(6) = fMRI_hmm8_music.fe;
fepc(7) = fMRI_hmm9_music.fe;
fepc(8) = fMRI_hmm10_music.fe;
fepc(9) = fMRI_hmm12_music.fe;
fepc(10) = fMRI_hmm15_music.fe;
fepc(11) = fMRI_hmm20_music.fe;


llc(1) = fMRI_hmm3_music.ll;
llc(2) = fMRI_hmm4_music.ll;
llc(3) = fMRI_hmm5_music.ll;
llc(4) = fMRI_hmm6_music.ll;
llc(5) = fMRI_hmm7_music.ll;
llc(6) = fMRI_hmm8_music.ll;
llc(7) = fMRI_hmm9_music.ll;
llc(8) = fMRI_hmm10_music.ll;
llc(9) = fMRI_hmm12_music.ll;
llc(10) = fMRI_hmm15_music.ll;
llc(11) = fMRI_hmm20_music.ll;

K_tests.fe = fepc;
K_tests.llc = llc;

%% Plot the FE values!
klabs = {'K3','K4','K5','K6','K7','K8','K9','K10','K12','K15','K20'};
figure

%fepc = K_tests.fe
%subplot(1,2,1)
bar(normalize(fepc))
grid on
xticklabels(klabs)
title('Free Energy K Comparison')
clear dex llc

%% Assemble the FO values into the new object
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

%% Run initial PLSes on the FO data
% in this study, we were looking and pre- and post-intervention data, so 
% we sliced up the data based on intervention status. First, we need to 
% get the data into PLS shape:

%% Get the intervention indices
idx_pre = fMRI_hmm4_music.pre_dex==1;
idx_post = fMRI_hmm4_music.pre_dex==0;

%% Slice it up
Pre_FO = cell(11,1);
Post_FO = cell(11,1);

for i = 1:11
    Pre_FO{i} = K_tests.FO{i}(idx_pre==1);
    Post_FO{i} = K_tests.FO{i}(idx_post==1);
end
K_tests.Pre_FO = Pre_FO;
K_tests.Post_FO = Post_FO;
clear i idx*

%% Run a PLS for every FO value
% This will loop through every K estimation's FO data, run a PLS,
% and save the PLS output into a new array:

addpath(genpath('Pls'))% see the ReadMe for where to download this
FO_data = K_tests;
FO_results = cell(11,1);

for k = 1:11
    indata = {K_tests.Pre_FO{k};K_tests.Post_FO{k}};% stack the data so that pre is on the top, post is on the bottom
    numsubs = nan(1,length(indata));
    for i = 1:length(indata)
        tempdata = indata{i};
        tempdata(isnan(tempdata)) = [];% remove any nans
        indata{i} = tempdata;
        numsubs(i) = size(tempdata,1);% get the number of subjects
    end
    clear i option
    option.method = 1;% option 1 for mean-centred PLS
    option.num_perm = 500;
    option.num_boot = 100;
    numcond = 2;% how many conditions we have (here, 2: pre- and post-intervention)

    K_res = pls_analysis(indata,numsubs,numcond,option);% run the PLS
    FO_results{k} = K_res;% save the results
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

% if you have a clear "winner" here, great! Proceed with that K. 
% Life is rarely that easy, though, so let's continue to explore:

%% Now correlate usc and vsc
% these are the row and column scores from the PLS. We can correlate them
% to get a "concensus" matrix showing us each PLS output's similarity:

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

% from here, I had two estimations (4 and 7) that were similar:

%% Try the dot product between similar estimations
% extract the state means for the similar estimations:

addpath(genpath('HMM-MAR-master'))
res7 = fMRI_hmm7_music;
m7 = getMean(res7.hmm);
res4 = fMRI_hmm4_music;
m4 = getMean(res4.hmm);
clear res7 res4

%% Normalize them:
m7 = normalize(m7,1);
m4 = normalize(m4,1);

%% Look at them:
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

%% Calculate the dot product:
mdot = m7'*m4;

% and look at it:
figure
imagesctxt(mdot)
colorbar
title('Dot Product: K = 7, K = 4')

% from here, I was able to conclude looking at the "loadings" in the matrix
% that K4 and K7 had very similar estimations (all of the K7 states were
% represented in K4)

%% For posterity, I also looked at the TP matrices:
% TP Matrices: Assemble everything into a new object
% we already have TP saved in the array, but juuust in case one were to run 
% the HMM estimation without extracting all the metrics, it's quick to re-estimate:

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

%% Interrogating the graph properties of the TP matrices:
% this section is a little more exploratory. The above tests were sufficient for us 
% to find an optimal K, but it's also possible to use graph metrics to explore 
% properties of the TP matrices and the state FC matrices (not shown here)

addpath(genpath('BCT'))
STR_test = cell(11,1);% strength
DEG_test = cell(11,1);% degree
CC_test = cell(11,1);% cclustering coefficient
for est = 1:11
    tempdata = K_tests.TP{est};
    thresh = prctile(abs(tempdata),50,1);%a very generous threshold of the 50th percentile
    tempdata(abs(tempdata)<thresh) = 0;
    str = strengths_dir(tempdata);
    [indeg,outdeg] = degrees_dir(tempdata);
    cc = clustering_coef_wd(tempdata);
    STR_test{est} = str;
    DEG_test{est} = [indeg;outdeg];
    CC_test{est} = cc;
end
clear est tempdata thresh str indeg outdeg cc

%% Look at specific estimations:
est = 2;
figure
%subplot(1,2,1)
bar(normalize(STR_test{est}))
grid on
ylabel('Strength')
xlabel('State')
title(sprintf('State Strength, %s',K_tests.labels{est}))

% subplot(1,2,2)
% bar(DEG_test{est})
% grid on
% ylabel('Strength')
% xlabel('State')
% title(sprintf('State Strength, %s',K_tests.labels{est}))
