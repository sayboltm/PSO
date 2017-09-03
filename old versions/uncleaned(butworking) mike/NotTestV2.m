close all
clear
clc

format shortEng
format compact

pth=pwd;
compfolder='\Composite Analysis\';
FEfolder='\Finite Element Analysis\';




%% Composite Analysis
%%%%%%%%%%%%%%%%%%%%%
% Objective:
% MAXIMIZE Stiffness (Ex)
% MINIMIZE Weight (Density)
%
% Design Vars:
% h
% FVF_i
% MaterialFiber_i
% MaterialMatrix_i, i = 1:6

% Beware the 'real' LayerOrientation is hard coded inside optimizationInputFunction
LayerOrientation=[45 -45 0 0 -45 45];
NumberOfLayers=length(LayerOrientation);

%OptimizationInputScript; % This script is where most of the input variables are defined

% Row 1 = 'a' (lower bound), Row 2 = 'b' (upper bound
% n columns for n dimensions. Weird but goes with PSO code better that way.
variable_limits = [1,8; 1,8; 1,8; 1,8; 1,8; 1,8; 1,3; 1,3; 1,3; 1,3; 1,3; 1,3; ...
    0.3,0.7; 0.3,0.7; 0.3,0.7; 0.3,0.7; 0.3,0.7; 0.3,0.7;0.1e-3,0.3e-3]';
% Add 3rd row for step?

%init_pts = [1,1,2,2,1,1, 1,1,2,2,1,1, 0.5,0.5,0.5,0.5,0.5,0.5, 0.15e-3];
init_pts = [1,1,1,1,1,1, 1,1,1,1,1,1, 0.3,0.3,0.3,0.3,0.3,0.3, 0.3e-3];
param_types = [1,1,1,1,1,1, 1,1,1,1,1,1, 0,0,0,0,0,0, 0];
%max_opt_iteration = 10;
max_PSO_iteration = 1000;%100;
pop_size = 500;%200;
pareto_p = .9;
dimensions = length(init_pts);
%C = [2.1, 2]; % Keep C(1) + C(2) >= 4 for stability
C = [3.9,0.2];
W = [1, 1];

ObjT = [1, -1]; % PSO maximizes, want Max first obj, but min second

[g, p, settings] = PSOpy_extendedV2(@compositeInputFunction ...
    ,variable_limits, param_types, C, W, ObjT, pop_size, max_PSO_iteration, pareto_p);
p = paretoFinder(g, p, settings, .9);
paretoPlotter( g, p, 'Ex', 'Density', 1,1,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATLABPOOL START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verLessThan('matlab','8.2.0')
    if matlabpool('size')==0
        set(findResource(),'ClusterSize',8);
        matlabpool open 8
    end
else
    if isempty(gcp('nocreate'))
        parpool;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Starting Beam Analysis'])
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
% passes to second wrapper function (feaInputFunctionV2)
[Sublaminate ] = secondStudyFunction( g, p, settings); 
variable_limits = [1,length(Sublaminate); 1,9; 1,9; 1,9; 1,9; 1,9; 1,9]';
param_types = [1, 3,3,3,3,3,3];
max_iterations=100;%10;
pop_size=50;%6;

ObjT = [-1, -1]; % Rest use same as before, hence using settings from output
[g2, p2, settings2] = PSOpy_extendedV2(@feaInputFunction ...
    ,variable_limits, param_types, settings.C, settings.W, ObjT, pop_size, max_iterations, settings.ParetoPercentage, Sublaminate);

p2 = paretoFinder_para(g2, p2, 8, settings2, .9); %poverride = .8 instead of .99

% Do old 'pareto points' for shits and giggs
p2 = paretoFinder(g2, p2, settings2, .9);

%paretoPlotter(g2, p2, 'Mass', 'MaxDeflection', 1, 1, 1)
paretoPlotter(g2, p2, 'MaxDeflection', 'Mass', 1, 1, 1)
%
% [g, p] = PSOpy_extendedV2(@feaInputFunctionV2 ...
%     ,variable_limits, param_types, C, W, ObjT, pop_size, max_PSO_iteration, pareto_p);

%[Un,MaxStress,MaxDeflection,Mass]=FiniteAnalysisExecution(Sublaminate,BeamRegions); % This is the FEA function
% The Outputs of this function are fairly self-explanatory

% fprintf('Maximum Deflection = '),disp([num2str(MaxDeflection*10^3),' (mm)'])
% fprintf('Maximum Stress = '),disp([num2str(MaxStress*10^(-6)),' (MPa)'])


%The following two lines of code are simply to accomodate the use of HEEDS.
% mass=Mass;
% deflection=MaxDeflection;

if verLessThan('matlab','8.2.0')
    if matlabpool('size') > 0
        %set(findResource(),'ClusterSize',8);
        matlabpool close 
    end
else
    if isempty(gcp('nocreate'))
        delete(gcp('nocreate'));
    end
end

% for i = 1:length(p2.properties)
%    A(i,:) = [p2.properties(i).base_response.MaxDeflection, p2.properties(i).base_response.Mass];
% end
% 
% for i = 1:length(p.properties)
%     B(i,1) = p.properties(i).Ex;
%     B(i,2) = p.properties(i).Density;
%     D(i,:) = p.properties(i).position;
% end
