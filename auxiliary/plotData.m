%%------------------------------- plotData ------------------------------%%
%
% Script to generate figures from precomputed data files.
% Indicies swapped for TR1<->TR2.
%
%--------------------------------------------------------------------------    
% Initial version: 08/13/19, J.B.

clc;
clear;

addpath ../data

wtest = warning('off','all');
 
% Directories for data storage
currentpath     = pwd;
figpath         = fullfile(currentpath,'..','/figs/');

% EXPERIMENT I 
% Selection of algorithms
% 1 : PDCO, 2: RSQP, 3: LMTR_SC, 4: LMTR_L2, 5: fmincon-ldl, 6: fmincon-cg,
% 7 : ipopt
indAlg          = [4 3 5 6 7];

leg={   'PDCO',...
        'RSQP',...
        'TR2',...
        'TR1',...       
        'FMIN-LDL', ...
        'FMIN-CG',...
        'IPOPT'};
                
types.colors    = ['k' 'm' 'b' 'r' 'm' 'c' 'g']; types.colors = types.colors(indAlg);
types.lines     = {'-', '-.', '-', '-.', '-', '-.', '-'}; types.lines = types.lines(indAlg);
types.markers   = ['o' 'o' 'o' 'o' 'o' 'o' 'o']; types.markers = types.markers(indAlg);

leglocation = 'NorthEast';

load('EXPERIMENT_I_L');
perf(exs(:,indAlg),times(:,indAlg),leg(indAlg),1,types,leglocation);    
print(gcf, '-dpsc2', fullfile(figpath,'time_EX_I_L.eps'));

perf(exs(:,indAlg),its(:,indAlg),leg(indAlg),1,types,leglocation);
print(gcf, '-dpsc2', fullfile(figpath,'iter_EX_I_L.eps'));

load('EXPERIMENT_I_S');
perf(exs(:,indAlg),times(:,indAlg),leg(indAlg),1,types,leglocation);    
print(gcf, '-dpsc2', fullfile(figpath,'time_EX_I_S.eps'));

perf(exs(:,indAlg),its(:,indAlg),leg(indAlg),1,types,leglocation);
print(gcf, '-dpsc2', fullfile(figpath,'iter_EX_I_S.eps'));

close all;

% EXPERIMENT II
leg={   'TR2',...
        'TR1',...
        'FMIN-LDL', ...
        'FMIN-CG', ...
        'IPOPT'};

types.colors    = ['b' 'r' 'm' 'c' 'g'];
types.lines     = {'-.', '-', '-.', '-' '-.'};
types.markers   = ['o' 'o' 'o' 'o' 'o'];

indAlg          = [2 1 3 4 5];
types.colors    = types.colors(indAlg);

leglocation     = 'NorthEast';

load('EXPERIMENT_II');
% Modify location of legend to have the same output as in article, i.e.
% change loaction to 'NorthEast' from 'SouthEast'.
perf(exs(:,indAlg),times(:,indAlg),leg(indAlg),1,types,leglocation);
print(gcf, '-dpsc2', fullfile(figpath,'time_EX_II.eps'));

perf(exs(:,indAlg),its(:,indAlg),leg(indAlg),1,types,leglocation);
print(gcf, '-dpsc2', fullfile(figpath,'iter_EX_II.eps'));

close all;

% EXPERIMENT III
leg={   'TR2',...
        'TR1',...       
        'FMIN-LDL',...
        'IPOPT'};
                
types.colors    = ['r' 'b' 'm' 'g']; %'k' 'y'
types.lines     = {'-.', '-', '-.','-'}; %'-',   '-'
types.markers   = ['o' 'o' 'o' 'o']; %'s' 's'

%indAlg          = [1 2 3 4 5]; %

indAlg          = [2 1 3 4];

load('EXPERIMENT_III');
perf(exs(:,indAlg),times(:,indAlg),leg(indAlg),1,types);
print(gcf, '-dpsc2', fullfile(figpath,'time_EX_III.eps'));

perf(exs(:,indAlg),its(:,indAlg),leg(indAlg),1,types);
print(gcf, '-dpsc2', fullfile(figpath,'iter_EX_III.eps'));

close all;

% EXPERIMENT IV
leg={   'TR2',...
        'TR1',...       
        'FMIN-LDL',...
        'IPOPT'};

% leg={'TR-$(\mathbf{P},\infty)$',...
%        'TR-$\ell_2$',...
%        'fmincon-ldl'};
                
types.colors    = ['r' 'b' 'm' 'g']; %'k' 'y'
types.lines     = {'-.', '-', '-.','-'}; %'-',   '-'
types.markers   = ['o' 'o' 'o' 'o']; %'s' 's'

%indAlg          = [1 2 3 4 5]; %

indAlg          = [2 1 3 4];

load('EXPERIMENT_IV');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1,types);
print(gcf, '-dpsc2', fullfile(figpath,'time_EX_IV.eps'));

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1,types);
print(gcf, '-dpsc2', fullfile(figpath,'iter_EX_IV.eps'));

% Comparison between proposed solvers only
indAlg          = [2 1];

perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1,types);
print(gcf, '-dpsc2', fullfile(figpath,'time_EX_IV_SEL.eps'));

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1,types);
print(gcf, '-dpsc2', fullfile(figpath,'iter_EX_IV_SEL.eps'));

close all;
    