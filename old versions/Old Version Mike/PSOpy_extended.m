function [ best_Ex, best_Dens, base_Ex, base_Density, globbest_position ...
    , pbest, randominit, globbest_fitness, pproperties ] = PSOpy_extended( f ...
    , dimensions, search_bounds, N, pop_size, pareto_p, init_pts)
%%
% Matlab port of Particle Swarm Optimization made in Python
%
% Usage:
%
% [globbest_fitness, globbest_position, N]=PSOpy(f, dimensions, a, b, N, pop_size)
%
% Inputs:
%       o f  - The objective function.  Pass with @function into here
%       o dimensions  - How many dimensions is the objective function?
%       o a  - Lower bound of search space
%       o b  - Upper bound of search space
%       o N  - Numer of iterations (optional, defaults to 100)
%       o pop_size  - Size of swarm of particles
%       o show_progress  - Show a graph? (2D functions ONLY)
%
% Outputs:
%       o globbest_fitness  - best function value obtained
%       o globbest_position  - best particle position
%       o N  - Iterations.  Shows default if none specified
%
%%

% Some constants for velocity updating
% Local and Global 'correction factors'
c1 = 2.1;
c2 = 2;

% Weight factors for outputs
w1 = 1;
w2 = 1;

% Percentage to show up on Pareto Front
%pareto_p = .8;

% Modifications to objective type: (since program minimizes)
ObjT1 = 1;%-1; % Since we want to max, not min
ObjT2 = -1;%1; % Keep as 1 to minimize Density (second output param)

% Error criterion for known functions (ignore/unused)
err_crit = 0.00001;

% Initialize the particles
ppos = zeros(pop_size,dimensions);% 'x values' (pos) for each particle
% Note for below: only need 'dimensions' rows because matlab vertcat issues
pfitness = zeros(pop_size,dimensions);% evaluated values for each particle
pv = zeros(pop_size,dimensions);% velocit of each particle
pbest = zeros(pop_size,dimensions);% best position for each particle

% Do random or given initial start?
% If initial start isn't given, init to random
if nargin<=6 || isempty(init_pts)
    for i = 1:pop_size
        % FOR each particle in swarm, init randpos in search space
        
        for j = 1:dimensions
            a = search_bounds(1,j);
            b = search_bounds(2,j);
            ppostemp(i,j) = (b-a).*rand(1,1) + a;
            %a = search_bounds(1,1), b = search_bounds(2,1)
            %ppos(i,:) = (b-a).*rand(1,dimensions) + a;
        end
        ppos(i,:) = filterFeasible(ppostemp(i,:), search_bounds, 'Hard');
    end
    randominit = ppos;
else
    for i = 1:pop_size
        ppos(i,:) = init_pts;%(b-a).*rand(1,dimensions) + a;
    end
end

%%% NOTE: the fitness elements are a bit weird.  Because it is desireable
%%% to group all globalbests together, they are the same size as position
%%% and velocity, although the rest will be zeros. Unless num output vars
%%% is equal to num input vars (dimensions)

% let the first particle be the global best
gbest = [ppos(1,:); pfitness(1,:); pv(1,:)];
err = 999999999; % Error, ignore it

progress = 0;
disp('Progress: 0%')
disp('[..........]')
i = 0;
while( i < N )
    % While less than itneration max
    
    % For each particle
    for j = 1:pop_size
        if i == 0
            % If first run, get a baseline
            %fitness(:) = f(ppos(j,:)); %some fitness = evalfuncat some xparams
            [base_Ex, base_Density] = f(ppos(j,:));
            Ex = base_Ex; Density = base_Density;
            fitness = w1*ObjT1*base_Ex/base_Ex + w2*ObjT2*base_Density/base_Density;
        else
            % Else runs, use baseline to normalize
            [Ex, Density] = f(ppos(j,:));
            fitness = w1*ObjT1*Ex/base_Ex + w2*ObjT2*Density/base_Density;
        end
        
        
        %fitness(:) = optimizationInputFunction(ppos(j,:))
        
        
        if fitness > pfitness(j,1:length(fitness)) % if this fitness LOWER (better) than particle fitness
            pfitness(j,1:length(fitness)) = fitness; % update the particle fitness
            pproperties(j,1) = fitness;
            pproperties(j,2:3) = [Ex, Density];
            pbest(j,:) = ppos(j,:); % update particle's best POS with this value
            
            %plot(pparams(1,1),pparams(1,2),'.')
            %text(pparams(1,1),pparams(1,2),num2str(j), 'Position', [pparams(1,1)+.1, pparams(1,2)+.1])
        end
        
        if fitness > gbest(2,1:length(fitness)) % if this fitness better than groups'
            %gbest = p % update group best
            gbest = [ppos(j,:); pfitness(j,:); pv(j,:)];
            best_Ex = Ex; 
            best_Dens = Density;
            
            %plot(ppos(1,1),ppos(1,2),'rx')
            %text(ppos(1,1),ppos(1,2),num2str(i), 'Position', [ppos(1,1)+.1, ppos(1,2)+.1])
            %text(pparams(1,1),pparams(1,2),num2str(j), 'Position', [pparams(1,1)+.1, pparams(1,2)+.1])
        end
        
        % Try K from:
        %http://web.itu.edu.tr/etaner/ea11.pdf
        phi = c1 + c2;
        K = 2/abs(2-phi-sqrt(phi^2-4*phi));
        % Update velocity
        %v = pv(j,:) + c1.*rand(1,dimensions).*(pbest(j,:)-ppos(j,:)) + c2.*rand(1,dimensions).*(gbest(1,:) - ppos(j,:));
        %pv(j,:) = pv(j,:) + c1.*rand(1,dimensions).*(pbest(j,:)-ppos(j,:)) + c2.*rand(1,dimensions).*(gbest(1,:) - ppos(j,:));
        pv(j,:) = K*(pv(j,:) + c1.*rand(1,dimensions).*(pbest(j,:)-ppos(j,:)) + c2.*rand(1,dimensions).*(gbest(1,:) - ppos(j,:)));
        
        % Update position
        %ppos(j,:) = ppos(j,:) + v;
        % Update position with feasible options
        %pposUnfiltered(j,:) = ppos(j,:) + v;
        %ppostest(j,:) = pposUnfiltered(j,:);
        %ppostest(j,:) = filterFeasible(pposUnfiltered(j,:), search_bounds, 'Hard');
        %ppos(j,:) = filterFeasible(ppos(j,:) + v, search_bounds, 'Hard');
        ppos(j,:) = filterFeasible(ppos(j,:) + pv(j,:), search_bounds, 'Hard');
    end
    
    
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

globbest_fitness = gbest(2,1:length(fitness));
globbest_position = gbest(1,:);
pfitnessout = pfitness(:,1);
%particle_fitness = pfitness(:,1:length(fitness));

%% Get points that made the cut to belong on Pareto Front and plot %%
%pareto_cut_index = 0;
figure()
hold on
for i = 1:pop_size
    if pproperties(i,1) >= pareto_p*globbest_fitness
        %pareto_cut_index = pareto_cut_index + 1;
        %pareto_cut(pareto_cut_index,:) = pproperties(i,:);
        plot(pproperties(i,2),pproperties(i,3),'x')
    end
end
plot(base_Ex, base_Density, 'r.')
plot(best_Ex, best_Dens, 'g.')
hold off
xlabel('Ex')
ylabel('Density')
%plot(pareto_cut(:,2), pareto_cut(:,3))

%% Local functions


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



