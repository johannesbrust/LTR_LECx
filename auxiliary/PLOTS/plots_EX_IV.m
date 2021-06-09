%-------------------- plots_EX_IV ----------------------------------------%
%
% Script to generate performance profiles for the experiment outcomes
% (Experiment IV)
%-------------------------------------------------------------------------%
% 12/11/20, J.B., Initial setup
% 06/08/21, J.B., Patch experiment data with reruns
% 06/09/21, J.B., Preparation for release

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
dataEX = load([datapath,'EXPERIMENT_IV_EXT']);
dataS  = load([datapath,'EXPERIMENT_IV_EXT_IPOPT_RERUN_SELECTED']);
dataSa = load([datapath,'EXPERIMENT_IV_EXT_IPOPT_RERUN_SELECTEDa']);
dataSb = load([datapath,'EXPERIMENT_IV_EXT_IPOPT_RERUN_SELECTEDb']);

% Patching column of IPOPT
colIP = 7;
for p=1:(size(dataEX.ex,1))
    
    if dataS.ex(p,colIP) == 1
        
       dataEX.ex(p,colIP) = dataS.ex(p,colIP);
       dataEX.numit(p,colIP) = dataS.numit(p,colIP);
       dataEX.t_aver(p,colIP) = dataS.t_aver(p,colIP);       
       
    elseif dataSa.ex(p,colIP) == 1
        
       dataEX.ex(p,colIP) = dataSa.ex(p,colIP);
       dataEX.numit(p,colIP) = dataSa.numit(p,colIP);
       dataEX.t_aver(p,colIP) = dataSa.t_aver(p,colIP);
       
    elseif dataSb.ex(p,colIP) == 1
        
       dataEX.ex(p,colIP) = dataSb.ex(p,colIP);
       dataEX.numit(p,colIP) = dataSb.numit(p,colIP);
       dataEX.t_aver(p,colIP) = dataSb.t_aver(p,colIP); 
       
    end
    
end

% Do not include one "smaller" problem
idxL = dataEX.nms(:,2)>500;

exs = dataEX.ex(idxL,:);
its = dataEX.numit(idxL,:);
times = dataEX.t_aver(idxL,:);

% Preparing and call to performance plot
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

% Axis annotation
ax = gca;
ax.XTick = [2^(-1), 1, 2, 2^2, 2^3, 2^(4), 2^5];
%XTickLabel = {'$2^{-8}$', '$2^{-4}$', '$1$', '$2^{4}$', '$2^{24}$'};
XTickLabel = {'2^{-1}', '1', '2', '2^{2}', '2^{3}', '2^{4}' , '2^{5}'};
set(ax,'XTickLabel',XTickLabel);

figname ='EXPERIMENT_IV_ITER';
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
ax.XTick = [2^(-1), 1, 2, 2^2, 2^3, 2^(4), 2^5];
%XTickLabel = {'$2^{-8}$', '$2^{-4}$', '$1$', '$2^{4}$', '$2^{24}$'};
XTickLabel = {'2^{-1}', '1', '2', '2^{2}', '2^{3}', '2^{4}' , '2^{5}'};
set(ax,'XTickLabel',XTickLabel);

figname ='EXPERIMENT_IV_TIME';
figname = [figname,'.pdf'];

print(fig,'-dpdf',fullfile(figpath,figname));
