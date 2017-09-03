close all
clear
clc

format shortEng
format compact

% Numer of cores to use for FEA
num_cores_fea = 8;

% Number of cores for other tasks (Pareto search)
% Since there may not be a break here, using all 8 cores will make machine
% unusable for anything else
num_cores_other = 7;
globallog.tstart = timeStamp(0);
%% Composite Analysis
%%%%%%%%%%%%%%%%%%%%%
% Objective:
% MAXIMIZE Stiffness (Ex)
% MINIMIZE Weight (Density)
%
% Design Vars:
% h
% FVF_i, i = 1:6
% MaterialFiber_i, i = 1:6
% MaterialMatrix_i, i = 1:6

disp(['Starting sublaminate analysis'])
timeStamp(1); % Display timestamp

% Beware the 'real' LayerOrientation is hard coded inside optimizationInputFunction
LayerOrientation=[45 -45 0 0 -45 45];
NumberOfLayers=length(LayerOrientation);

% Row 1 = 'a' (lower bound), Row 2 = 'b' (upper bound
% n columns for n dimensions. Weird but goes with PSO code better that way.
variable_limits = [1,8; 1,8; 1,8; 1,8; 1,8; 1,8; 1,3; 1,3; 1,3; 1,3; 1,3; 1,3; ...
    0.3,0.7; 0.3,0.7; 0.3,0.7; 0.3,0.7; 0.3,0.7; 0.3,0.7;0.1e-3,0.3e-3]';
% Add 3rd row for step? No is taken care of by filterFeasible

param_types = [1,1,1,1,1,1, 1,1,1,1,1,1, 0,0,0,0,0,0, 0];

max_iterations = 1000;%100;
pop_size = 500;%200;
pareto_p = .9;

%C = [2.1, 2]; % Keep C(1) + C(2) >= 4 for stability
C = [3.9,0.2];

% Weights of objectives
W = [1, 1];

% Maximize the first objective (Ex), minimize second (Density)
ObjT = [1, -1];
globallog.CompAnalysisStart = timeStamp(0);
% Run optimizer for comp analysis
[g, p, settings] = PSOpy_extendedV2(@compositeInputFunction ...
    ,variable_limits, param_types, C, W, ObjT, pop_size, max_iterations, pareto_p);
timeStamp(1);
globallog.CompAnalysisFinish = timeStamp(0);
globallog.Paretosearch1Start = timeStamp(0);

% Find pareto points from first analysis
p = paretoFinder_para(g, p, num_cores_other);

% Find old 'pareto points' for shits and giggs
p = paretoFinder(g, p, settings, .9);
globallog.Paretosearch1Finish = timeStamp(0);
% Plot results
paretoPlotter( g, p, 'Ex', 'Density', 1,1,1)
globallog.Paretoplot1Finish = timeStamp(0);

timeStamp(1);
globallog.BeamAnalysisStart = timeStamp(0);
disp(['Starting beam analysis'])
%% Beam Analysis
%%%%%%%%%%%%%%%%
% Objective:
% MINIMIZE Mass
% MINIMIZE MaxDeflection
%
% Design Vars:
% LaminateComposition (from above)
% NumberOfLaminates: N_j, j = 1:10



% Extract valid sublaminates from first results, pass to optimizer which
% passes to second wrapper function (feaInputFunction)
globallog.secondStudyFunctionStart = timeStamp(0);
[Sublaminate ] = secondStudyFunction( g, p, settings); 
globallog.secondStudyFunctionFinish = timeStamp(0);
variable_limits = [1,length(Sublaminate); 1,9; 1,9; 1,9; 1,9; 1,9; 1,9]';

param_types = [1, 3,3,3,3,3,3];

max_iterations=100;%10;
pop_size=50;%6;
% Rest use same as before, hence using settings from output

% Both objectives (Mass, MaxDeflection) should be minimized
ObjT = [-1, -1]; 


% Run optimizer for beam w/ FEA
timeStamp(1);
globallog.FEAStart = timeStamp(0);
% Start matlabpool if not already
matlabpoolStart(num_cores_fea)
[g2, p2, settings2] = PSOpy_extendedV2(@feaInputFunction ...
    ,variable_limits, param_types, settings.C, settings.W, ObjT, pop_size, max_iterations, settings.ParetoPercentage, Sublaminate);
% Stop matlabpool
matlabpoolStop();
timeStamp(1);
globallog.FEAStop = timeStamp(0);

% Find Pareto points
globallog.Paretosearch2Start = timeStamp(0);
p2 = paretoFinder_para(g2, p2, num_cores_other); 

% Do old 'pareto points' for shits and giggs
p2 = paretoFinder(g2, p2, settings2, .9);%poverride = .8 instead of .99
globallog.Paretosearch2Finish = timeStamp(0);
%paretoPlotter(g2, p2, 'Mass', 'MaxDeflection', 1, 1, 1)
paretoPlotter(g2, p2, 'MaxDeflection', 'Mass', 1, 1, 1)
globallog.Paretoplot2Finish = timeStamp(0);

timeStamp(1);
globallog.tstop = timeStamp(0);