function [ p ] = paretoFinder( g, p, settings, poverride )
%UNTITLED3 Finds the pareto points and outputs modified p with a
%p.properties.Pareto object where 1 signifies a Pareto point
% Has override if you want a different ParetoPercentage than in settings

% WARNING: This is depricated and WRONG.  Not actual definition of a Pareto
% point.

if ((nargin <= 3) || (isempty(poverride)))
    pareto_coef = settings.ParetoPercentage;
else
    pareto_coef = poverride;
end

pareto_p = 1-pareto_coef; % This is NOT settings.ParetoPercentage

pop_size = length(p.properties);
paretocounter = 0;

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

