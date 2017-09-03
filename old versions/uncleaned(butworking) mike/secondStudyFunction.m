function [Sublaminate ] = secondStudyFunction( g, p, settings) 
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Take the results from 1st study, and any other inputs, and feed into
% optimizer.  Output results

% Would be cool to plot the globalbest from the first one at the end
% First, build a set of laminates from the Pareto front points:

% sublam_fields = fields(p.properties);
% pop_size = length(p.properties);
% j = 0;
% for i = 1:pop_size;
%     if p.properties(i).Pareto == 1;
%         j = j + 1;
%         for k = 1:length(sublam_fields)
%             Sublaminate(j).(sublam_fields{k}) = p.properties(i).(sublam_fields{k});
%         end
%     end
% end

% Version for actual Pareto points
ObjFields = fields(p.storage.trial(1).Objectives);
NonObjFields = fields(p.storage.trial(1).NonObjectives);
num_nonobj = length(NonObjFields);
num_obj = length(ObjFields);
j = 0;
for i = 1:length(p.storage.trial)
    if p.storage.trial(i).paretoactual == 1;
        j = j + 1;
        for k = 1:length(ObjFields)
        Sublaminate(j).(ObjFields{k}) = p.storage.trial(i).Objectives.(ObjFields{k});
        end
        for k = 1:(num_nonobj)
            Sublaminate(j).(NonObjFields{k}) = p.storage.trial(i).NonObjectives.(NonObjFields{k});
        end
        i
        Sublaminate(j).position = p.storage.trial(i).position;
    end  
end



% Sublaminate.variable_limits = [1,length(Sublaminate); 1,9; 1,9; 1,9; 1,9; 1,9; 1,9]';
% 
% Sublaminate.param_types = [1, 3,3,3,3,3,3];

%Sublaminate.ObjT = [-1, -1];

% PSOpy_extendedV2(@feaInputFunctionV2 ...
%     ,variable_limits, param_types, settings.C, settings.W, ObjT, pop_size, settings.MaxIterations, settings.ParetoPercentage);


end
