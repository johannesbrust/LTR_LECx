%-------------------- plots_EX -------------------------------------------%
%
% Script to generate performance profiles for the experiment outcomes
%
%-------------------------------------------------------------------------%
% 12/11/20, J.B., Initial setup

% Load perf_ext_fnc
%addpath(fullfile(pwd,'..','..','/auxiliary'));
addpath(genpath('../../auxiliary'));
addpath(genpath('../../ssget'));

% Paths
currentpath = pwd;
datapath = fullfile(currentpath,'..','..','/data/');
figpath = fullfile(currentpath,'..','..','/figs/');

% Labels
leg={   'TR2',...
        'TR2H',...
        'TR2L',...
        'TR1',...       
        'TR1H',...
        'TR1L',...        
        'IPOPT'};
                
types.colors    = ['b' 'r' 'm' 'g' 'k' 'c' 'y']; %'k' 'y'
types.lines     = {'-.', '-', '-','-.','-','-','-'}; %'-',   '-'
types.markers   = ['o' 'o' 'o' 'o' 'o' 'o' 'o']; %'s' 's'

%indAlg          = [1 2 3 4 5 6 7];
indAlg          = [4 5 6 1 2 3 7];

% Load data from experiment
%datapath = fullfile(pwd,'..','..','/data/EXPERIMENT_III_EXT');
dataEX = load([datapath,'EXPERIMENT_III_EXT']);

% Data stored as (only rows 1-50 of the problems were used in experiment):
% dataEX = 
% 
%     convs: [50x7 double]
%       exs: [50x7 double]
%       its: [50x7 double]
%       nbs: [50x7 double]
%      numf: [50x7 double]
%      numg: [50x7 double]
%     times: [50x7 double]

exs = dataEX.exs;
its = dataEX.its;
times = dataEX.times;

% Patch certain experiments with data from rerun
dataEX_CHK = load([datapath,'EXPERIMENT_III_EXT_IPOPTCHK']);

% Checklist of problems
clist={'beacxc';...
    'lp_cre_d';...
    'lp_fit2d';...
    'lp_scsd1'};
data = load('EVEN_N_PROBS');
probIdx = data.probIdx;
nprob = 50;
selIdx = zeros(nprob,1,'int8');
% 
for i = 1:nprob
    % Loading problem and some data    
    Prob= ssget(probIdx(i));
    
    pname_L = Prob.name;
    pname_SP= strsplit(pname_L,'/');
    pname = pname_SP{2};
    
    if sum(strcmp(pname,clist))==0
        continue; % Skip problem if not in clist
    else
        selIdx(i) = 1;
    end
end

% Patch selected problems
selIdxF = find(selIdx);
exs(selIdxF,:)=dataEX_CHK.exs(selIdxF,:);
its(selIdxF,:)=dataEX_CHK.its(selIdxF,:);
times(selIdxF,:)=dataEX_CHK.times(selIdxF,:);

% Preparing and call to performance plot

% if isempty(legLocation) == true
%     legLocation='SouthEast';
% end
legLocation='SouthEast';
ticks   = -1;
ticke   = 5;
XTick   = 2.^(ticks:ticke);
XLim    = [XTick(1),XTick(end)];
%perf_fnc(exs(:,indAlg),its(:,indAlg),leg(indAlg),1,types,'SouthEast',taumax);

% Iterations
perf_ext_fnc(exs(:,indAlg),its(:,indAlg),leg(indAlg),1,types,legLocation,XTick,XLim);
title('ITER');
% Saving/formatting plots

fig                     = gcf;
fig.PaperPositionMode   = 'auto';
fig_pos                 = fig.PaperPosition;
fig.PaperSize           = [fig_pos(3) fig_pos(4)];

figname ='EXPERIMENT_III_ITER';
figname = [figname,'.pdf'];

print(fig,'-dpdf',fullfile(figpath,figname));

% Times
perf_ext_fnc(exs(:,indAlg),times(:,indAlg),leg(indAlg),1,types,legLocation,XTick,XLim);
title('TIME');
% Saving/formatting plots

fig                     = gcf;
fig.PaperPositionMode   = 'auto';
fig_pos                 = fig.PaperPosition;
fig.PaperSize           = [fig_pos(3) fig_pos(4)];

figname ='EXPERIMENT_III_TIME';
figname = [figname,'.pdf'];

print(fig,'-dpdf',fullfile(figpath,figname));
