close all
clear
clc

format shortEng
format compact

pth=pwd;
compfolder='\Composite Analysis\';
FEfolder='\Finite Element Analysis\';

%% Composite Analysis
LayerOrientation=[45 -45 0 0 -45 45];
NumberOfLayers=length(LayerOrientation);

OptimizationInputScript; % This script is where most of the input variables are defined

cd([pth,compfolder]);

Sublaminate=CompositeAnalysis(FiberMaterial,MatrixMaterial,FvF,LayerThickness,LayerOrientation); % This is the laminate analysis function
% The output "Sublaminate" is a structure containing the following fields:
%   Ex, Ey, Gxy, vxy,vyx, ABD, Thickness, and Density;
% The most important outputs are "Ex", "Thickness", and "Density". This
% should be used in the composite optimization Study

cd(pth)

%% Beam Analysis
cd([pth,FEfolder])

[Un,MaxStress,MaxDeflection,Mass]=FiniteAnalysisExecution(Sublaminate,BeamRegions); % This is the FEA function
% The Outputs of this function are fairly self-explanatory

% fprintf('Maximum Deflection = '),disp([num2str(MaxDeflection*10^3),' (mm)'])
% fprintf('Maximum Stress = '),disp([num2str(MaxStress*10^(-6)),' (MPa)'])
cd(pth)

%The following two lines of code are simply to accomodate the use of HEEDS.
mass=Mass;
deflection=MaxDeflection;