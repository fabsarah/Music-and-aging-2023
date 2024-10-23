%% HMM ESTIMATION
% This script describes extracting HMM and preliminary data extraction,
% and visualization from HMM outputs. There are a lot of different ways to
% runn HMM that depend on the type of data and how it's organized, and the
% HMM-MAR Toolbox Wiki is quite helpful: 
% https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide

%% Step 1: get our giant timeseries vector
% To runHMM, you need a timeseries divided into cells. Cells can be 
% participaants, tasks, participants*task. For this analysis, I ran HMM on
% participant-level data (each person's continuous timeseries for the
% entire scan run). 

load('Indata.mat');%get the data
%% 1.1: Parameter Selection
clear T
K = 4; % number of states for this HMM run
ndim = 220; % number of channels
N = 1; % number of trials (songs not repeated, so N = 1 for me)
T = cell(length(Indata),1);
for i = 1:length(Indata)
    T{i} = size(Indata{i},1);
end
clear i
%% 1.2: Setting the initialization parameters
clear options
options = struct();
options.K = K; 
options.Fs = 2.1053; %sampling rate
options.covtype = 'full';
options.zeromean = 0;
options.standardise = 0;
options.standardise_pc = 1;
options.DirichletDiag = 10;
options.pca = 0.95;% 95% of variance captured
options.useParallel=1;% requires the parallel computing toolbox
options.verbose = 1;
%% Step 2: Run the HMM! 
%This takes several hours locally and will kick up an error if T and Indata
%are not the same length
addpath(genpath('HMM-MAR-master'))
[hmm, Gamma, Xi, vpath,~,~,~] = hmmmar(Indata,T,options);
% Sending the output to a structure array
fMRI_hmm4_music.hmm = hmm;
fMRI_hmm4_music.Gamma = Gamma;
fMRI_hmm4_music.vpath = vpath;
fMRI_hmm4_music.Xi = Xi;
filename = '/path_to_directory/fMRI_hmm4_music.mat';
save(filename, 'fMRI_hmm4_music');
%% Populating the structure with our input parameters and region labels
fMRI_hmm4_music.T = T;
fMRI_hmm4_music.options = options;
fMRI_hmm4_music.labs = HMMfMRI.ST_labs;
%% Step 3: Get te data!
% The rest of this script is extracting code and running initial checks and
% visualization on HMM outputs. We'll be sending everything into a giant
% structure array...SAVE THIS OFTEN.
%% 3.1: Get the FC and TP matrices!
FC = cell(1,fMRI_hmm4_music.options.K);%initialize the output cell. We need K cells

for i = 1:length(FC)
    [~,S] = getFuncConn(fMRI_hmm4_music.hmm,i); %get that sweet, sweet FC matrix
    FC{i} = S; %populate the FC cell
end

TProbs = getTransProbs(fMRI_hmm4_music.hmm);

fMRI_hmm4_music.FC = FC; %add FC to the structure array
fMRI_hmm4_music.Tprobs = TProbs;%add TP to the array
clear i S FC TProbs % make look nice
%% PLOTTING: TP and FC matrices:

K = fMRI_hmm4_music.options.K;
dims = [2,3];%dimensions of the figure

figure 
subplot(dims(1),dims(2),1)%first item in the figure space
TProbs = fMRI_hmm4_music.TProbs;
imagesc(TProbs)
pbaspect([1 1 1])%lock aspect ratio...make it nice and square
colorbar
title(sprintf('Transitional Probabilities, %d states',K),'FontSize',14)
xticks(1:K)
yticks(1:K)
xlabel('State')
ylabel('State')

for i = 1:K
subplot(dims(1),dims(2),i+1)%start at 2 to not over-write the TP plot
imagesc(fMRI_hmm4_music.FC{i})
pbaspect([1 1 1])
colorbar
yticks(55:110:220)
yticklabels({'Left','Right'})
xticks(55:110:220)
xticklabels({'Left','Right'})
title(sprintf('State %d',i),'FontSize',14)
end
clear i S res K %make look nice

%% 3.2: Switching Rate and Fractional Occupancy:
data = fMRI_hmm4_music.Gamma;%getting the relevant HMM output
options = fMRI_hmm4_music.options;
T = fMRI_hmm4_music.T;
SR = getSwitchingRate(data,T,options);
FO = getFractionalOccupancy(data,T,options,2);
%This will give us long column vectors of data which we can break up later
fMRI_hmm4_music.Metrics.SR = SR;% put it in the array!
fMRI_hmm4_music.Metrics.FO = FO;
clear data options T temp* SR FO start i dex array% make look nice

%% 3.3: Grand mean of FC states, and mean centering each state
K = fMRI_hmm4_music.options.K;
FC = cat(3,fMRI_hmm4_music.FC{:});%rearranging the cell into a 3-d array
for i = 1:K
FC(:,:,i) = atanh(FC(:,:,i));
end
clear i
fMRI_hmm4_music.mean_FC = mean(FC,3);%mean
fMRI_hmm4_music.std_FC = std(FC,0,3);%standard deviation

mc_FC = nan(220,220,K);%setting the mean output matrices
for i = 1:K
    data = atanh(fMRI_hmm4_music.FC{1,i});
    tempnew = (data-fMRI_hmm4_music.mean_FC)./fMRI_hmm4_music.std_FC;%mean centering!
    mc_FC(:,:,i) = tempnew;
end
fMRI_hmm4_music.mc_FC = mc_FC;% put it in the structure array!

clear data i tempnew mc_FC FC
%% PLOTTING: Grand mean and mean-centred FC matrices!

K = fMRI_hmm4_music.options.K;
figure 

subplot(3,5,1)
imagesc(fMRI_hmm4_music.mean_FC)
pbaspect([1 1 1])
colorbar
title(sprintf('Grand Mean FC, K = %d',K),'FontSize',14)
yticks(55:110:220)
yticklabels({'Left','Right'})
xticks(55:110:220)
xticklabels({'Left','Right'})

for i = 1:K
subplot(3,5,i+1)
imagesc(fMRI_hmm4_music.mc_FC(:,:,i))
pbaspect([1 1 1])
colorbar
yticks(55:110:220)
yticklabels({'Left','Right'})
xticks(55:110:220)
xticklabels({'Left','Right'})
title(sprintf('State %d, Mean-Centred',i),'FontSize',14)
end
clear i S res K

%% 3.4: Vpath and Gamma by participant and task:
% Reorganizing the vpath and Gamma
T = fMRI_hmm4_music.T;
nparts = length(T);% the number of participants/cells in your original T file

vpath_Sub = cell(nparts,1);
vdata = fMRI_hmm4_music.vpath;
gdata = fMRI_hmm4_music.Gamma;

start = 1;%this is a long vector, so we need a counter to keep track of where we are
for i = 1:length(vpath_Sub)
    stop = start+(T{i}-1);%get the stopping point from the T vector
    tempdata = vdata(start:stop);%get the relevant slice of data
    vpath_Sub{i} = tempdata';%add it to the new cell
    start = stop+1;%update the starting point
end
clear i start stop tempdata

%% Break Gamma into a cell array
% now we're going to do the same thing for Gamma
Gamma_Sub = cell(nparts,1);
start = 1;
for j = 1:nparts
    stop = start+(T{j}-1);
    tempdata = gdata(start:stop,:);
    start = stop+1;
    Gamma_Sub{j} = tempdata';
    disp(j)% this can take a while
end
clear i G start stop j tempdata ntask

fMRI_hmm4_music.Gamma_Sub = Gamma_Sub;
fMRI_hmm4_music.vpath_Sub = vpath_Sub;

%% Getting the vpath divided up by task
% Here, I convert the vpath_Sub cell array into a matrix 
music_vpath = cell2mat(vpath_Sub);
fMRI_hmm4_music.music_vpath = music_vpath;
%% DATA CHECK:
% check the data!
x = reshape(music_vpath',1,[]);
x = x';
y = vdata;
if isequal(x,y)
    disp('The data are good!')
else
    disp('The data is full of spiders - audit this!')
end
clear s piece* x y
clear vpath_Sub vdata start i dex stop newdata index 

%% 3.5: Life Times!
data = fMRI_hmm4_music;
Gamma = data.vpath;
T = data.T;
options = data.options;
lifetimes = getStateLifeTimes(Gamma, T, options, [],[],0);

fMRI_hmm4_music.LT = lifetimes;

%% 3.6 Gamma and FO by piece
Gamma = fMRI_hmm4_music.Gamma_Sub;
K = fMRI_hmm4_music.options.K;
T = fMRI_hmm4_music.T;
options = fMRI_hmm4_music.options;
nparts = length(T);
rest_FO = nan(nparts,K);

for part = 1:length(Gamma)
    FO = getFractionalOccupancy(Gamma{part}',T{part},options);
    rest_FO(part,:) = FO;
end
clear part FO

%% Using the above output to get task-wise data
Gamma_music = Gamma;
block_FO = cell(nparts,1);
ntask = 24;
time = 60;%the length of the excerpt as an integer

for part = 1:length(Gamma_music)
    tempmat = nan(ntask,K);
    tempdata = Gamma_music{part};
    tempdata = reshape(tempdata,K,time,ntask);
    for piece = 1:ntask
        FO = getFractionalOccupancy(tempdata(:,:,piece)',time,options);
        tempmat(piece,:) = FO;
    end
    Block_FO{part} = tempmat;
end

%% DATA CHECK
part = 69;
x = Gamma_music{part};
y = reshape(x,K,time,24);
piece = 24;
if isequal(x(:,1440-(time-1):1440),y(:,:,piece))%manually update this...super crunchy, sorry :D
    disp('Acceptable!')
else
    disp('Unacceptable!')
end
% we're good to go!
clear x y piece part ntask time

%% If the data are good, add them to the array:
fMRI_hmm4_music.Gamma_music = Gamma_music;
fMRI_hmm4_music.block_FO = Block_FO;

clear data K gdata hmm Gamma* dims filename FO Indata lifetimes music_vpath... 
    N ndim npieces options rest_FO T temp* TProbs vpath Xi Block_FO

%% PLOTTING: Vpath
plotdata = fMRI_hmm4_music.music_vpath;
K = fMRI_hmm4_music.K;

tick_labs = {'State 1';'State 2';'State 3';'State 4'};
figure
imagesc(plotdata)
colormap(parula(K))
ch = colorbar;
ch.Ticks = (1:K);
ch.TickLabels = tick_labs;
ylabel('Participants')
xlabel('Time')
title('Music VPath','FontSize',12)

clear plotdata K ch
%% DON'T FORGET TO SAVE THE ARRAY.
% We now have most of what we need to do the analysis. See the subsequent
% scripts for the analysis code. 

