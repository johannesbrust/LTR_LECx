%-------------------- plots_EX_V -----------------------------------------%
%
% Script to generate performance profiles for the experiment outcomes
%
%-------------------------------------------------------------------------%
% 12/11/20, J.B., Initial setup
% 06/04/21, J.B., Update for most recent experiment ("Va")
% 06/07/21, J.B., Update to generate plots for release
% 06/09/21, J.B., Prepration of file for release

% Load perf_ext_fnc
%addpath(fullfile(pwd,'..','..','/auxiliary'));
addpath(genpath('../../auxiliary'));
%addpath(genpath('../../ssget'));

% Paths
currentpath = pwd;
datapath = fullfile(currentpath,'..','..','/data/');
figpath = fullfile(currentpath,'..','..','/figs/');
probpath = fullfile(currentpath,'..','/..','/auxiliary/');

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
dataEX = load([datapath,'EXPERIMENT_V_EXT']);

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

exs = dataEX.ex;
its = dataEX.numit;
times = dataEX.t_aver;


% Patch iteration outcomes "0"
its(its(:,7)==0,7) = 0.5;

% Preparing and call to performance plot
legLocation='SouthEast';
ticks   = -2;
ticke   = 10;
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

% Axis annotation
ax = gca;
ax.XTick = [2^(-2), 1, 2^2, 2^4, 2^(6), 2^8, 2^10];
%XTickLabel = {'$2^{-8}$', '$2^{-4}$', '$1$', '$2^{4}$', '$2^{24}$'};
XTickLabel = {'2^{-2}', '1', '2^{2}', '2^{4}', '2^{6}', '2^{8}' , '2^{10}'};
set(ax,'XTickLabel',XTickLabel);

figname ='EXPERIMENT_Va_ITER';
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

% Axis annotation
ax = gca;
ax.XTick = [2^(-2), 1, 2^2, 2^4, 2^(6), 2^8, 2^10];
%XTickLabel = {'$2^{-8}$', '$2^{-4}$', '$1$', '$2^{4}$', '$2^{24}$'};
XTickLabel = {'2^{-2}', '1', '2^{2}', '2^{4}', '2^{6}', '2^{8}' , '2^{10}'};
set(ax,'XTickLabel',XTickLabel);


figname ='EXPERIMENT_Va_TIME';
figname = [figname,'.pdf'];

print(fig,'-dpdf',fullfile(figpath,figname));
