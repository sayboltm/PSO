function [ Objectives, Constraints, NonObjectives ] = feaInputFunction( args_in, Sublaminate_in)
%UNTITLED Simplifies entry into PSO algorithm
% 19 dimensional search space (LayerOrientation not included)
% In a nutshell, this provides a modified version of the objective
% function, including only what is needed or what will be varied.

% Wrapper function is passed the resulting Sublaminate, containing valid
% sublaminates that made up the Pareto front in the first optimization

% Args in should be = length(Sublaminate_in.variable_limits)
% for i = 1:length(Sublaminate_in)
%     switch i
%         case 1
%             Sublaminate(i)
%         
%     end
% Sublaminate input selected from first args_in which *should* be a valid
% sublaminate because optimizer was fed Sublaminate_in.variable_limits
Sublaminate = Sublaminate_in(args_in(1));
% Stuff from end of OptimizationInputScript.m
Reg1=args_in(2);%1;
Reg2=args_in(3);%1;
Reg3=args_in(4);%1;
Reg4=args_in(5);%1;
Reg5=args_in(6);%1;
Reg6=args_in(7);%1;

% Combines the sublaminate numbers into a single vector
BeamRegions=[Reg1 Reg2 Reg3 Reg4 Reg5 Reg6];


%Sublaminate = Ex, Density, Ey, Gxy, vxy, vyx, ABD, Thickness, [fitness]

[Un,MaxStress,MaxDeflection,Mass,Domain]=FiniteAnalysisExecution(Sublaminate,BeamRegions);
NonObjectives.Un = Un;
NonObjectives.MaxStress = MaxStress;
Objectives.MaxDeflection = MaxDeflection;
Objectives.Mass = Mass;
NonObjectives.Domain = Domain;

% Sublaminate=CompositeAnalysis(FiberMaterial,MatrixMaterial,FvF,LayerThickness, LayerOrientation);
% Objectives.Ex = Sublaminate.Ex;
% NonObjectives.Ey = Sublaminate.Ey;
% NonObjectives.Gxy = Sublaminate.Gxy;
% NonObjectives.vxy = Sublaminate.vxy;
% NonObjectives.vyx = Sublaminate.vyx;
% NonObjectives.ABD = Sublaminate.ABD;
% NonObjectives.Thickness = Sublaminate.Thickness;
% Objectives.Density = Sublaminate.Density;
 Constraints.null = 0;


end

