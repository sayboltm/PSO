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

cd([pth,compfolder]);

% Row 1 = 'a' (lower bound), Row 2 = 'b' (upper bound
% n columns for n dimensions. Weird but goes with PSO code better that way.
variable_limits = [1,8; 1,8; 1,8; 1,8; 1,8; 1,8; 1,3; 1,3; 1,3; 1,3; 1,3; 1,3; ...
    0.3,0.7; 0.3,0.7; 0.3,0.7; 0.3,0.7; 0.3,0.7; 0.3,0.7;0.1e-3,0.3e-3]';
% Add 3rd row for step?

%init_pts = [1,1,2,2,1,1, 1,1,2,2,1,1, 0.5,0.5,0.5,0.5,0.5,0.5, 0.15e-3];
init_pts = [1,1,1,1,1,1, 1,1,1,1,1,1, 0.3,0.3,0.3,0.3,0.3,0.3, 0.3e-3];

%max_opt_iteration = 10;
max_PSO_iteration = 1000;
pop_size = 200;
pareto_p = .99;
dimensions = length(init_pts);
%[globbest_fitness, globbest_position, p_fitness] = PSOpy_extended(@optimizationInputFunction, ...
%        dimensions, variable_limits, max_PSO, pop_size, 1, init_pts);
[best_Ex, best_Dens, base_Ex, base_Density, globbest_position, pbest ...
    , randinitpos, globbest_fitness, pfitness] = PSOpy_extended(@optimizationInputFunction, ...
    dimensions, variable_limits, max_PSO_iteration, pop_size, pareto_p);
% i = 1;
% while (i <= max_opt_iteration)
%     [globbest_fitness(i,:), globbest_position(i,:), p_fitness(i)] = PSOpy_extended(@optimizationInputFunction, ...
%         dimensions, variable_limits, max_PSO, pop_size, 1);
%     disp(['Iterations complete: ', num2str(i),'/',num2str(max_opt_teration)])
%     i = i + 1;
% end
cd(['..']);
return     

% call PSOpy_ext @optInputFunc

%////STOPPED 4/24/16/////%
% Arent all FvF always the same?  Layer orientation always the same? Y
% Why FvF have sq brackets
%[Ex, Ey, Gxy, vxy,vyx, ABD, Thickness, Density] = CompAnalysis(
% {FibMat1:FibMat6}, {MtrxMat1:MtrxMat6}, [FvF1:FvF6], LayerThickness,
% [layer_orientation]
%Sublaminate=CompositeAnalysis(FiberMaterial,MatrixMaterial,FvF,LayerThickness,LayerOrientation); % This is the laminate analysis function
% The output "Sublaminate" is a structure containing the following fields:
%   Ex, Ey, Gxy, vxy,vyx, ABD, Thickness, and Density;
% The most important outputs are "Ex", "Thickness", and "Density". This
% should be used in the composite optimization Study

cd(pth)

%% Beam Analysis
%%%%%%%%%%%%%%%%
% Objective:
% MINIMIZE Stress
% MINIMIZE Weight
%
% Design Vars:
% LaminateComposition (from above)
% NumberOfLaminates: N_j, j = 1:10

cd([pth,FEfolder])

[Un,MaxStress,MaxDeflection,Mass]=FiniteAnalysisExecution(Sublaminate,BeamRegions); % This is the FEA function
% The Outputs of this function are fairly self-explanatory

% fprintf('Maximum Deflection = '),disp([num2str(MaxDeflection*10^3),' (mm)'])
% fprintf('Maximum Stress = '),disp([num2str(MaxStress*10^(-6)),' (MPa)'])
cd(pth)

%The following two lines of code are simply to accomodate the use of HEEDS.
mass=Mass;
deflection=MaxDeflection;