%% Start Here: figure out piece-wise T
%addpath(genpath('HMM-MAR-master'))
StimT = [50;10];
StimT = repmat(StimT,(80*24),1);
%%
data = fMRI_hmm4_music;
vpath = data.vpath;%(1:947*80);
T = data.T(1:80);
T = [cell2mat(T);StimT];%this was the key!
hmm = data.hmm;
Gamma = data.Gamma;
Xi = data.Xi;

RestMask = 1:(947*80);
RestMask = reshape(RestMask,[],80)';
RestMask = num2cell(RestMask,2);
start = (947*80);
MusicMask = 1:(1440*80);
MusicMask = MusicMask+start;
MusicMask = reshape(MusicMask,[],80)';
MusicMask = num2cell(MusicMask,2);
StimMask = cell(length(StimT),1);
%% Get the vpath indices into piece order
tempMask = cell(length(MusicMask),1);
for i = 1:length(MusicMask)
    tempdata = MusicMask{i};
    tempdata = reshape(tempdata,60,[])';%reshape the data row-wise
    tempMask{i} = tempdata;
end
tempMask = cell2mat(tempMask);
clear i tempdata
%% Break them into music+ratings
stimdex = 1:2:length(StimMask);
for i = 1:length(tempMask)
    tempdata = tempMask(i,:);
    StimMask{stimdex(i)} = tempdata(1:50);
    StimMask{stimdex(i)+1} = tempdata(51:60);
end
clear i stimdex
%% Get the masked TP Matrices
outMask = [RestMask;StimMask];
PStim = getMaskedTransProbMats(vpath,T,hmm,outMask,Gamma,Xi);
PStim = PStim';
clear vpath Gamma Xi T hmm RestMask outMask MusicMask
%% Just music and ratings
fMRI_hmm4_music.TP.TP_Music = PStim(81:2:end);
fMRI_hmm4_music.TP.TP_Rate = PStim(82:2:end);
%% Slice 'em up!
newmat = reshape(TP_Music,24,80);
x = newmat{1,3};
y = TP_Music{49};
isequal(x,y)
%%
%here, we're going to use the block_dex table (tells us which block is
%which for stimuli, ratings, etc) to mke some new matrices
part_TPs = cell(80,1);
for part = 1:80
    tempdata = newmat(:,part);
    vecmat = nan(24,25);%preallocating the vectorized matrix
    for tune = 1:length(tempdata)
        vecmat(tune,:) = reshape(tempdata{tune}',1,[]);
    end
    part_TPs{part} = vecmat;
end
clear part tune vecmat tempdata
%% Test! The! Transform!
% ie did I vectorize what I think I vectorized
part = 80;
piece = 23;
x = newmat{piece,part};
y = part_TPs{part}(piece,:);
y = reshape(y,5,5)';

if isequal(x,y)
    disp('Cake and tea!')
else
    disp('Suffering and woe!')
end
clear x y part piece
%%
tempdata = fMRI_hmm4_music;
%transform cells into participant-wise matrices! 
stim_TP = cell(80,4);
for part = 1:80
    index = tempdata.block_dex{part};
    block = table2array(index(:,2));
    tempTP = part_TPs{part};%need to arrange this into blocks
    for stim = 1:4
        newmat = tempTP(block==stim,:);
        stim_TP{part,stim} = newmat;
    end
end
fMRI_hmm4_music.TP.stim_TP = stim_TP;
fMRI_hmm4_music.TP.part_TPs = part_TPs;
clear part index block tempTP stim newmat
% NB cell 65 has no data (NANs in the input)
%% Now, average and assemble into four more matrices
% still needs to be checked for indexing errors!
tempdata = fMRI_hmm4_music.TP.stim_TP;
stim_TP_all = cell(4,1);
for stimcat = 1:4
    stimmat = nan(80,25);
    for part = 1:80
        parttp = tempdata{part,stimcat};
        if isempty(parttp)
            continue
        else
                tempTP = nanmean(parttp);
                stimmat(part,:) = tempTP;
        end
    end
    stim_TP_all{stimcat} = stimmat;
end
fMRI_hmm4_music.TP.stim_TP_all = stim_TP_all;
clear stimcat part stimmat part parttp state tempTP
%% Slice it up again!
idx_pre = fMRI_hmm4_music.pre_dex==1;
idx_post = fMRI_hmm4_music.pre_dex==0;
idx_old = fMRI_hmm4_music.Ages>40;
idx_young = fMRI_hmm4_music.Ages<40;
match_pre = [1;3;5;8;10;9;11;14;15;21;26;27;29;33;37];
match_post = [2;4;7;12;13;16;20;23;24;25;28;30;32;35;38];
%% Keep on dicing!
Time_Stim = cell(4,3);
for i = 1:4
    Time_Stim{i,1} = stim_TP_all{i}(idx_pre==1&idx_young==1,:);%pre young
    Time_Stim{i,2} = stim_TP_all{i}(idx_pre==1&idx_old==1,:);%pre old
    Time_Stim{i,3} = stim_TP_all{i}(idx_post==1,:);%post (old)
end
Match_TP = cell(4,2);
for i = 1:4
    Match_TP{i,1} = stim_TP_all{i}(match_pre,:);
    Match_TP{i,2} = stim_TP_all{i}(match_post,:);
end
clear i idx_pre idx_post idx_young idx-old x match_pre match_post ans
%%
%addpath(genpath('Pls'))
clear option
%indata = {cell2mat(Time_Stim(:,1))};%young pre
%indata = {cell2mat(Time_Stim(:,2))};%old pre
%indata = {cell2mat(Time_Stim(:,3))};%old post
%indata = {cell2mat(Time_Stim(:,1));cell2mat(Time_Stim(:,2))};%young vs old pre
indata = {[Time_Stim{1,1};Time_Stim{1,2};Time_Stim{2,1};Time_Stim{2,2};...
    Time_Stim{3,1};Time_Stim{3,2};Time_Stim{4,1};Time_Stim{4,2}]};
for i = 1:length(indata)
    tempdata = indata{i};
    tempdata(any(isnan(tempdata),2),:)=[];
    indata{i} = tempdata;
end

option.method = 1;
option.num_perm = 500;
option.num_boot = 100;

%TP_young_res = pls_analysis(indata,40,4,option);
%TP_old_res = pls_analysis(indata,24,4,option);
%TP_old_post_res = pls_analysis(indata,15,4,option);
%TP_all_res = pls_analysis(indata,[40,24],4,option);
TP_over_res = pls_analysis(indata,64,4,option);
%% Plot with Bars

res = TP_over_res;
LV = 1;
p = res.perm_result.sprob(LV);
title_var = (sprintf('PLS TP: Stimuli All Ages (Pre), p = %d',p));

figure
subplot(1,2,1)
z = res.boot_result.orig_usc(:,LV);
bar(z)
hold on
yneg = res.boot_result.llusc(:,LV);
ypos = res.boot_result.ulusc(:,LV);
errorbar(1:length(z),z,yneg-z,ypos-z,'.')
colorbar off
xticklabels({'Self-Select';'Popular';'Novel';'BP'})
%xticklabels({'Self-Select Y';'Popular Y';'Novel Y';'BP Y';...
%    'Self-Select O';'Popular O';'Novel O';'BP O'})
xtickangle(45)
grid on
title(title_var,'FontSize',12)

x = res.boot_result.compare_u(:,LV);
x = reshape(x,5,5)';
subplot(1,2,2)
shower_tile_plot(x);
pbaspect([1 1 1])
xticks([])
yticks([])
ylabel('States')
colorbar
title(sprintf('LV %d',LV),'FontSize',12)
% Alright, well there is an effect. Hot damn! One significant LV for each
% age group as well
%%
%addpath(genpath('Pls'));
clear option indata

indata = {[cell2mat(Match_TP(:,1));cell2mat(Match_TP(:,2))]};
option.method = 1;
option.num_perm = 500;
option.num_boot = 100;

match_res = pls_analysis(indata,15,8,option);
%% Go up a block to plot young and old. Stay here to plot everyone together
res = match_res;
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
xticklabels({'Pre SS';'Pre Pop';'Pre Nov';'Pre BP';...
    'Post SS';'Post Pop';'Post Nov';'Post BP'})
xtickangle(45)
grid on
title(sprintf('PLS: Stim Pre and Post, p = %d',p),'FontSize',12)

x = res.boot_result.compare_u(:,LV);
x = reshape(x,5,5)';
%x(abs(x)<6) = 0;
subplot(1,2,2)
shower_tile_plot(x)
pbaspect([1 1 1])
xticks([])
yticks([])
ylabel('States')
colorbar
title(sprintf('LV %d',LV),'FontSize',12)
%% Start Here Plot some TP matrices
labels = {'Self-Select Y';'Popular Y';'Novel Y';'BP Y';...
    'Self-Select O';'Popular O';'Novel O';'BP O';...
    'Self-Select O Post';'Popular O Post';'Novel O Post';'BP O Post'};

figure
start = 1;
for time = 1:3
    for task = 1:4
    subplot(3,4,start)
    plotdata = nanmean(Time_Stim{task,time});
    shower_tile_plot(reshape(plotdata,5,5));
    pbaspect([1 1 1])
    caxis([0 0.1])
    colorbar
    title(labels{start},'FontSize',12)
    start = start+1;
    end
end
clear i start plotdata
%%
subplot(2,2,2)
plotdata = mean(rest_match_post);
shower_tile_plot(reshape(plotdata,5,5));
pbaspect([1 1 1])
caxis([0 0.1])
colorbar
title('Post Rest')

subplot(2,2,3)
plotdata = mean(music_match_pre);
shower_tile_plot(reshape(plotdata,5,5));
pbaspect([1 1 1])
caxis([0 0.1])
colorbar
title('Pre Music')

subplot(2,2,4)
plotdata = mean(music_match_post);
shower_tile_plot(reshape(plotdata,5,5));
pbaspect([1 1 1])
caxis([0 0.1])
colorbar
title('Post Music')
%% Save the data to the array
fMRI_hmm4_music.TP.Rest = P(1:80);
fMRI_hmm4_music.TP.Music = P(81:end);
fMRI_hmm4_music.TP.PLS.results.all_res = all_res;
fMRI_hmm4_music.TP.PLS.results.young_res = young_res;
fMRI_hmm4_music.TP.PLS.results.old_res = old_res;
fMRI_hmm4_music.TP.PLS.results.old_post_res = old_post_res;
fMRI_hmm4_music.TP.PLS.results.match_res = match_res;
%% Plotting young and old TP matrices
figure
subplot(2,3,1)
plotdata = mean(Pre_Rest_young);
shower_tile_plot(reshape(plotdata,5,5));
pbaspect([1 1 1])
caxis([0 0.1])
colorbar
title('Pre Rest Young')

subplot(2,3,2)
plotdata = mean(Pre_Rest_old);
shower_tile_plot(reshape(plotdata,5,5));
pbaspect([1 1 1])
caxis([0 0.1])
colorbar
title('Pre Rest Old')

subplot(2,3,4)
plotdata = mean(Pre_Music_young);
shower_tile_plot(reshape(plotdata,5,5));
pbaspect([1 1 1])
caxis([0 0.1])
colorbar
title('Pre Music Young')

subplot(2,3,5)
plotdata = mean(Pre_Music_old);
shower_tile_plot(reshape(plotdata,5,5));
pbaspect([1 1 1])
caxis([0 0.1])
colorbar
title('Pre Music Old')

subplot(2,3,3)
plotdata = mean(Post_Rest_old);
shower_tile_plot(reshape(plotdata,5,5));
pbaspect([1 1 1])
caxis([0 0.1])
colorbar
title('Post Rest Old')

subplot(2,3,6)
plotdata = mean(Post_Music_old);
shower_tile_plot(reshape(plotdata,5,5));
pbaspect([1 1 1])
caxis([0 0.1])
colorbar
title('Post Music Old')
