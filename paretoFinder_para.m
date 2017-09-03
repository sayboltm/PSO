function [ p ] = paretoFinder_para( g, p, num_cores)
%UNTITLED3 Finds the pareto points and outputs modified p with a
% Parallelized, updated version

% TODO: Update definitions

% Threshold of [total points to process]^2 to use parallel processing
para_threshold = 5000; 

% Total number of points to process
tot_points = length(p.storage.trial);

% Get objectives, since Pareto is front of best of them
% Currently only supporst 2 obj at a time
ObjFields = fields(p.storage.trial(1).Objectives); % 1 chosen arbitrarily

para = 0;
if tot_points^2 >= para_threshold
    % Use parallel processing
    para = 1;
    timeStamp(1);
    disp(['paretoFinder: Using ', num2str(num_cores), ' cores for parallel searching'])
    
end

if para == 0
    % Use non-parallelized code
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
    
    
    
elseif para == 1
    % Use parallel processing
    
    % Start matlabpool
    matlabpoolStart(num_cores)
    
    % Init Pareto status to one
    % Matlab won't prealloc this for some reason, but can slowly allocate
    % it (see below loop to set values)
%     pparetoactual = ones(tot_points);
%     
%     pfitObjFld1 = zeros(tot_points);
%     pfitObjFld2 = zeros(tot_points);
    
   
    % parfor loops don't like structures, so load data into regular double
    for i = 1:tot_points
        pparetoactual(i) = 1;
        pfitObjFld1(i) = p.storage.trial(i).fitness_separated.(ObjFields{1});
        pfitObjFld2(i) = p.storage.trial(i).fitness_separated.(ObjFields{2});
    end
    
    timeStamp(1);
    disp(['paretoFinder: Data loaded, processing beginning....'])
    
    parfor i = 1:tot_points
        for j = 1:tot_points
            if ((pfitObjFld1(i) < pfitObjFld1(j)) && (pfitObjFld2(i) < pfitObjFld2(j)))
                pparetoactual(i) = 0; % IS NOT a pareto
            end
        end
    end
    
    timeStamp(1);
    disp(['paretoFinder: Processing complete! Assigning Pareto points...'])
    
    % Assign the acquired values back
    for i = 1:tot_points
        p.storage.trial(i).paretoactual = pparetoactual(i);
    end
    
    % Stop matlabpool
    matlabpoolStop();
    timeStamp(1);
    disp(['paretoFinder: Pareto points assigned.'])
    
end

end

