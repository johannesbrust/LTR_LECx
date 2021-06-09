%-------------------- READEME_EX_III.txt ---------------------------------%
%
% Readme file for Experiment III of this Software (called Experiment I) 
% in the article
%
% "Large-Scale Optimization With Linear Equality Constraints Using
% Reduced Compact Representation", J.J.Brust, R.F.Marica, C.G.Petra and
% M.A.Saunders (2021)
%
%-------------------------------------------------------------------------%
% 06/09/21, J.B., Initial version

EXPERIMENT
To run this experiment you can type 

    >> EXPERIMENT_III

OUTPUTS
To generate the figures from the manuscript you can use

    auxiliary/PLOTS/plots_EX_III.m

To generate the tables you can type

    auxiliary/TABLES/table_EX_III.m

By default the figures are stored in the /figs folder with the corresonding
experiment name. A text file with the name 'table_*' is stored in
the /data folder.

ADDITIONAL
Information on dense columns for the test problems from this 
experiment are stored in "auxiliary/TABLES/DENSE_COLS.mat". 
(The script tests/PROBS_DENSECOLS.m was used to generate this dense
column data)
"EVEN_N_PROBS.mat" stores the ID's of the problems used in this experiment.
