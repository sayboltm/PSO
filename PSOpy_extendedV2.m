function [ g, p, settings ] = PSOpy_extendedV2( f , search_bounds, param_types, C, W, ObjT, pop_size, N, pareto_p, input_struct, init_pts)
%%
% Matlab port of Particle Swarm Optimization made in Python
% Universal Modular Version
%
% Usage:
%
% [g, p, settings]=PSOpy_extendedV2(f, search_bounds, param_types, C, W, ObjT, pop_size, N, pareto_p, [input_struct], [init_pts])
%
% Inputs:
%       o f  - The objective function.  Pass with @function into here
%       o search_bounds  - Range of search space for each var
%       o param_types  - Type of values variables are valid for
%       o C  - Local and Global weight constants
%       o W  - Weight factor for each objective
%       o ObjT  - Objective type for each objective (-1 min, 1 max)
%       o pop_size  - Number of particles in swarm
%       o N  - Numer of iterations (optional, defaults to 100)
%       o pareto_p [Depricated]  -  pareto percentage var
%       o input_struct [Optional]  - Any additional info objective function
%       needs should be passed through this input struct
%       o init_points [Depricated] - Init point array 
%
% Outputs:
%       o g  - Struct of global results
%       o p  - Struct of particle results and results of every computation
%       o settings - Struct of settings used
%
%%

% Some constants for velocity updating
% Local and Global 'correction factors'
c1 = C(1); c2 = C(2);
%c1 = 2.1;
%c2 = 2;
settings.C = C;
settings.W = W;
settings.ObjT = ObjT;
settings.MaxIterations = N;
settings.ParetoPercentage = pareto_p;
settings.PopSize = pop_size;
% Get dimensions (NEED TO CHANGE IF NOT ALL VARS HAVE BOUNDS)
dimensions = length(search_bounds);

% Weight factors for outputs
%w1 = W(1); w2 = W(2); % Fix

% Percentage to show up on Pareto Front
%pareto_p = .8;

% Modifications to objective type: (since program minimizes)
%ObjT1 = 1;%-1; % Since we want to max, not min
%ObjT2 = -1;%1; % Keep as 1 to minimize Density (second output param)

%num_responses = length(W); %same as objectives....
% Error criterion for known functions (ignore/unused)
err_crit = 0.00001;

% Initialize the particles
p.pos = zeros(pop_size,dimensions);% 'x values' (pos) for each particle
p.fitness([1:1:pop_size]) = -Inf;% =  zeros(pop_size,1);% fitness values for each particle
p.v = ones(pop_size,dimensions).*rand(pop_size,dimensions);% velocit of each particle
%p.bestpos = zeros(pop_size,dimensions);% best position for each particle

% Do random or given initial start?
% If initial start isn't given, init to random
if nargin<=10 || isempty(init_pts)
    for i = 1:pop_size
        % FOR each particle in swarm, init randpos in search space
        
        for j = 1:dimensions
            a = search_bounds(1,j);
            b = search_bounds(2,j);
            ppostemp(i,j) = (b-a).*rand(1,1) + a;
            %a = search_bounds(1,1), b = search_bounds(2,1)
            %ppos(i,:) = (b-a).*rand(1,dimensions) + a;
        end
        p.pos(i,:) = filterFeasibleV2(ppostemp(i,:), search_bounds, param_types, 'Sticky');
    end
    p.randominit = p.pos;
    
else
    for i = 1:pop_size
        p.pos(i,:) = init_pts;%(b-a).*rand(1,dimensions) + a;
    end
end

%%% NOTE: the fitness elements are a bit weird.  Because it is desireable
%%% to group all globalbests together, they are the same size as position
%%% and velocity, although the rest will be zeros. Unless num output vars
%%% is equal to num input vars (dimensions)

% let the first particle be the global best
%gbest = [ppos(1,:); pfitness(1,:); pv(1,:)];          
%             for k = 1:num_objectives
%                 g.bestobj(k) = Objectives.(ObjFields{k})
%             end
p.bestpos = p.pos;
g.bestpos = p.pos(1,:);
g.bestfitness = -Inf; %struct.empty; %p.fitness(1);
g.bestvelocity = struct.empty; %p.v(1,:);

err = 999999999; % Error, ignore it

progress = 0;
disp('Progress: 0%')
disp('[..........]')
m = 0; % Total iteration counter
i = 0;
while( i < N )
    % While less than itneration max
    
    % For each particle
    for j = 1:pop_size
        m = m + 1;
        if i == 0
            % If first run, get a baseline, compute fitness differently
            % Also store input_struct if exists;
            %fitness(:) = f(ppos(j,:)); %some fitness = evalfuncat some xparams
            if ((nargin<=9) || (isempty(input_struct)))
                [Objectives, ~, NonObjectives] = f(p.pos(j,:));
            else
                % The objective function needs some struct
                [Objectives, ~, NonObjectives] = f(p.pos(j,:),input_struct);
                p.storage.input_struct = input_struct;
            end
            
            num_objectives = length(fields(Objectives));
           
            num_nonobjectives = length(fields(NonObjectives));
            
            ObjFields = fields(Objectives);
            NonObjFields = fields(NonObjectives);
            for k = 1:num_objectives
                base_response(k) = Objectives.(ObjFields{k});
                % Save base_response for each particle
                p.properties(j).base_response.(ObjFields{k}) = base_response(k);
                
                % Even a base response might be on Pareto front %
                % Commented, itll get this after the if statements
                %p.storage.trial(m).Objectives.(ObjFields{k}) = base_response(k);
                fitmat(k) = W(k)*ObjT(k)*base_response(k)/base_response(k);
                p.storage.trial(m).fitness_separated.(ObjFields{k}) =  fitmat(k);
            end
            %Ex = base_Ex; Density = base_Density;
            %WARNING: This will need to be changed around when constraints
            %added
%             for k = 1:num_nonobjectives
%                 p.storage.trial(m).NonObjectives.(NonObjFields{k}) = NonObjectives.(NonObjFields{k});
%             end
            %fitness = w1*ObjT1*base_Ex/base_Ex + w2*ObjT2*base_Density/base_Density;
            fitness = sum(fitmat);
            
        else
            % Else runs, use baseline to normalize
            %[Ex, Density] = f(ppos(j,:));
            %fitness = w1*ObjT1*Ex/base_Ex + w2*ObjT2*Density/base_Density;           
            if ((nargin<=9) || (isempty(input_struct)))
                [Objectives, ~, NonObjectives] = f(p.pos(j,:));
            else
                % The objective function needs some struct
                [Objectives, ~, NonObjectives] = f(p.pos(j,:),input_struct);
            end
            for k = 1:num_objectives
                response(k) = Objectives.(ObjFields{k});
                %p.storage.trial(m).Objectives.(ObjFields{k}) = response(k);
            end
            %Ex = base_Ex; Density = base_Density;
            %WARNING: This will need to be changed around when constraints
            %added
            for k = 1:num_objectives
                fitmat(k) = W(k)*ObjT(k)*response(k)/base_response(k);
                p.storage.trial(m).fitness_separated.(ObjFields{k}) =  fitmat(k);
            end

            
            %fitness = w1*ObjT1*base_Ex/base_Ex + w2*ObjT2*base_Density/base_Density;
            fitness = sum(fitmat);
        end
        
        
        % Store Everything for current iteration: (fitness_separated
        % already stored)
        p.storage.trial(m).fitness = fitness;
        p.storage.trial(m).pfitness = p.fitness(j);
        p.storage.trial(m).position = p.pos(j,:);
        
        for k = 1:num_objectives
            p.storage.trial(m).Objectives.(ObjFields{k}) = Objectives.(ObjFields{k});
        end
        for k = 1:num_nonobjectives%num_objectives+1:num_nonobjectives+num_objectives
            p.storage.trial(m).NonObjectives.(NonObjFields{k}) = NonObjectives.(NonObjFields{k});
        end
       
        
        
        if fitness > p.fitness(j) % if this fitness better than particle fitness
            p.fitness(j) = fitness; % update the particle fitness
            %p.storage.trial(m).pfitness = p.fitness(j);
            
            % Record all shit in particle's 'properties'
            p.properties(j).fitness = fitness;
            
            for k = 1:num_objectives
                %p.properties(j,k) = Objectives.(ObjFields{k});
                p.properties(j).(ObjFields{k}) = Objectives.(ObjFields{k});
                %p.storage.trial(i+j).Objectives.(ObjFields{k}) = Objectives.(ObjFields{k});
            end
            for k = 1:num_nonobjectives%num_objectives+1:num_nonobjectives+num_objectives
                %p.properties(j,k) = NonObjectives.(NonObjFields{k-num_objectives});
                p.properties(j).(NonObjFields{k}) = NonObjectives.(NonObjFields{k});
                %p.storage.trial(i+j).NonObjectives.(NonObjFields{k}) = NonObjectives.(NonObjFields{k});
            end
            % And the position
            p.properties(j).position = p.pos(j,:);
            %p.storage.trial(i+j).position = p.pos(j,:);
            
            %pproperties(j,2:3) = [Ex, Density];
            p.bestpos(j,:) = p.pos(j,:); % update particle's best POS with this value
            %p.storage.trial(i+j).bestpos = p.pos(j,:);
            %plot(pparams(1,1),pparams(1,2),'.')
            %text(pparams(1,1),pparams(1,2),num2str(j), 'Position', [pparams(1,1)+.1, pparams(1,2)+.1])
        end
        
        if fitness > g.bestfitness %gbest(2,1:length(fitness)) % if this fitness better than groups'
            %gbest = p % update group best
            %gbest = [ppos(j,:); pfitness(j,:); pv(j,:)];
            g.bestpos = p.pos(j,:);
            g.bestfitness = p.fitness(j);
            g.bestvelocity = p.v(j,:);
            
            
            
            for k = 1:num_objectives
                %g.bestproperties(k) = Objectives.(ObjFields{k});
                g.bestproperties.(ObjFields{k}) = Objectives.(ObjFields{k});
            end
            
            for k = 1:num_nonobjectives
                %g.bestproperties(k) = NonObjectives.(NonObjFields{k-num_objectives});
                g.bestproperties.(NonObjFields{k}) = NonObjectives.(NonObjFields{k});
            end
            g.bestproperties.fitness = fitness; % Or p.fitness.. should be same thing
            g.bestproperties.position = p.pos(j,:);
            %best_Ex = Ex; 
            %best_Dens = Density;
            
            %plot(ppos(1,1),ppos(1,2),'rx')
            %text(ppos(1,1),ppos(1,2),num2str(i), 'Position', [ppos(1,1)+.1, ppos(1,2)+.1])
            %text(pparams(1,1),pparams(1,2),num2str(j), 'Position', [pparams(1,1)+.1, pparams(1,2)+.1])
        end
        
        % Try K from:
        %http://web.itu.edu.tr/etaner/ea11.pdf
        phi = c1 + c2;
        K = 2/abs(2-phi-sqrt(phi^2-4*phi));
        % Update velocity
        p.v(j,:) = K*(p.v(j,:) + c1.*rand(1,dimensions).*(p.bestpos(j,:)-p.pos(j,:)) + c2.*rand(1,dimensions).*(g.bestpos - p.pos(j,:)));
        
        % Update position
        p.pos(j,:) = filterFeasibleV2(p.pos(j,:) + p.v(j,:), search_bounds, param_types, 'Sticky');
    end
    
    % Unused
    if err < err_crit
        break
    end
    %progress bar. '.' = 10%
    if (mod(i,(N/10)) == 0)
        progress = progress+1;
        progressBar(progress)
    end
    i  = i+1;
    
end

%% Get points that made the cut to belong on Pareto Front and plot %%
% This is old, useless crap but is being kept because $/GB is cheap and
% might be useful (and faster) 'Pareto' definition for noncompeting
% objectives.

% figure()
% hold on
% % Since we can only visualize 2-3 variables, this is hardcoded
% for j = 1:pop_size
%     p.properties(j).Pareto = 0; % Create new parameter for each particle
%     if p.properties(j).fitness >= pareto_p*g.bestfitness
%         plot(p.properties(j).Ex,p.properties(j).Density,'x')
%         p.properties(j).Pareto = 1;
%     end
% end
% plot(base_response(1), base_response(2), 'r.')
% plot(g.bestproperties.Ex, g.bestproperties.Density, 'g.')
% 
% hold off
% xlabel('Ex')
% ylabel('Density')

% Since above commented out, at least assign the 'Pareto' property. scratch
% that


%% Local functions
% Check out the switch statement in the paretoPlotter for some real crafty
% comparisons!

    function progressBar(progress)
        switch progress
            case 1
                disp('[|.........]')
            case 2
                disp('[||........]')
            case 3
                disp('[|||.......]')
            case 4
                disp('[||||......]')
            case 5
                disp('[|||||.....]')
            case 6
                disp('[||||||....]')
            case 7
                disp('[|||||||...]')
            case 8
                disp('[||||||||..]')
            case 9
                disp('[|||||||||.]')
            case 10
                disp('[||||||||||]')
                disp('Progress: 100%!')
            otherwise
                disp('progressBar Error: Invalid percentage')
        end
    end


end



