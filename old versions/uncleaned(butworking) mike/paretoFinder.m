function [ p ] = paretoFinder( g, p, settings, poverride )
%UNTITLED3 Finds the pareto points and outputs modified p with a
%p.properties.Pareto object where 1 signifies a Pareto point
% Has override if you want a different ParetoPercentage than in settings
if ((nargin <= 3) || (isempty(poverride)))
    pareto_coef = settings.ParetoPercentage;
else
    pareto_coef = poverride;
end

pareto_p = 1-pareto_coef; % This is NOT settings.ParetoPercentage

pop_size = length(p.properties);
paretocounter = 0;

%%%%%%%%%%%%%%%%%%%%%%%
% Get objectives, since Pareto is front of best of them
% Currently only supporst 2 obj at a time
ObjFields = fields(p.storage.trial(1).Objectives); % 1 chosen arbitrarily


% Actual defn of Pareto (since we have competing objectives)
tot_points = length(p.storage.trial);
ten_percent = tot_points/10; % Might speed up by only doing division once
for i = 1:tot_points
    p.storage.trial(i).paretoactual = 1;
    for j = 1:tot_points
        if ((p.storage.trial(i).fitness_separated.(ObjFields{1}) < p.storage.trial(j).fitness_separated.(ObjFields{1})) && (p.storage.trial(i).fitness_separated.(ObjFields{2}) < p.storage.trial(j).fitness_separated.(ObjFields{2})))
           p.storage.trial(i).paretoactual = 0; % IS NOT a pareto 
        end
    end
    if (mod(i,(ten_percent)) == 0)
        disp(['paretoFinder: Points processed: ', num2str(i), '/', num2str(tot_points)])
    end
end
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Pareto points
for j = 1:pop_size
    p.properties(j).Pareto = 0; % Create new parameter for each particle
    
    %if p.properties(j).fitness >= pareto_coef*g.bestfitness
    if abs(p.properties(j).fitness) >= abs(g.bestfitness) - abs(pareto_p*g.bestfitness)
        p.properties(j).Pareto = 1;
        paretocounter = paretocounter + 1;
    end
    
end

if paretocounter == 0;
    disp(['WARNING: NO Pareto points found.  Increase tolerance'])
end

end

