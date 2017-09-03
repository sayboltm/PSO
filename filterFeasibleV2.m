function [out] = filterFeasible(in, search_bounds, param_type, boundary_type)
% V2.0: takes param_type, where is length of search_bounds and specifies if
% int or not. 1 for int, 3 for oddint, 0/else for float.
%
% Usage:
%
% [out] = filterFeasible(in, search_bounds, boundary_type)
%
% Inputs:
%       o in  - 1xN array of positions to 'fix'
%       o search_bounds  - 2xN array of limits of the positions
%       o boundary_type  - String matching one of the boundary types below
%
% Outputs:
%       o out  - Returns modified 'in' array
%       o globbest_position  - best particle position
%       o N  - Iterations.  Shows default if none specified
%
% Boundary_types:
%       o 'Sticky' Boundary:    Move to boundary
%       o 'Periodic' Boundary:  Make search space dimensions 'round' and
%               come through the other side if a particle passes over.
%       o 'Hard' Boundary: Bounce off by amount of violation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out = zeros(size(in));
switch boundary_type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% STICKY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Sticky'
        for i = 1:length(in)
            % Process and make feasible
            if param_type(i) == 1
                % param_type 1 must be type int AND within range (a,b)
                a = search_bounds(1,i); b = search_bounds(2,i);
                intver = round(in(i));
                if ((intver < a) || (intver > b)) % If out of bounds,
                    if intver < a
                        increaseby = a-intver;
                        out(i) = intver + increaseby;
                    elseif intver > b
                        decreaseby = b-intver;
                        out(i) = intver + decreaseby;
                    end
                else
                    out(i) = intver;
                end
                
            elseif param_type(i) == 3
                % param_type 3: ODD INT
                a = search_bounds(1,i); b = search_bounds(2,i);
                intver = round_oddint(in(i));
                if ((intver < a) || (intver > b)) % If out of bounds,
                    if intver < a
                        increaseby = a-intver;
                        out(i) = intver + increaseby;
                    elseif intver > b
                        decreaseby = b-intver;
                        out(i) = intver + decreaseby;
                    end
                else
                    out(i) = intver;
                end
            else
                a = search_bounds(1,i); b = search_bounds(2,i);
                if ((in(i) < a) || (in(i) > b)) % If out of bounds,
                    if in(i) < a
                        increaseby = a-in(i);
                        out(i) = in(i) + increaseby;
                    elseif in(i) > b
                        decreaseby = b-in(i);
                        out(i) = in(i) + decreaseby;
                    end
                else
                    out(i) = in(i);
                end
                
            end
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% PERIODIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Periodic'
        
        for i = 1:length(in)
            % Process and make feasible
            if param_type(i) == 1
                
                a = search_bounds(1,i); b = search_bounds(2,i);
                intver = round(in(i));
                if ((intver < a) || (intver > b)) % If out of bounds,
                    modify = 1;
                    jj = 0;
                    %coef_restitution = 1; % Can't be used for int
                    int_c_rest = 1;
                    %situations
                    while modify == 1
                        if intver < a
                            amount_out = a - intver;
                            intver = b - amount_out*int_c_rest;
                            if ((intver >= a) && (intver <= b))
                                modify = 0;
                                out(i) = intver;
                            end
                            
                        elseif intver > b
                            amount_out = intver - b;
                            intver = a + amount_out*int_c_rest;
                            if ((intver >= a) && (intver <= b))
                                modify = 0;
                                out(i) = intver;
                            end
                        end
                        jj = jj + 1;
                        if (jj >= 10)
                            disp(['Search bounds for var: ',...
                                num2str(i),...
                                ' might be too tight. Cannot fix an ',...
                                'int type. Breaking....'])
                            int_c_rest = 1/2;
                            %coef_restitution = 9/jj;
                        end
                    end
                else % If not out of bounds, just use rounded version
                    out(i) = intver;
                end
            elseif param_type(i) == 3
                
                a = search_bounds(1,i); b = search_bounds(2,i);
                intver = round_oddint(in(i));
                if ((intver < a) || (intver > b)) % If out of bounds,
                    modify = 1;
                    jj = 0;
                    %coef_restitution = 1; % Can't be used for int
                    int_c_rest = 1;
                    %situations
                    while modify == 1
                        if intver < a
                            amount_out = a - intver;
                            intver = b - amount_out*int_c_rest;
                            if ((intver >= a) && (intver <= b))
                                modify = 0;
                                out(i) = intver;
                            end
                            
                        elseif intver > b
                            amount_out = intver - b;
                            intver = a + amount_out*int_c_rest;
                            if ((intver >= a) && (intver <= b))
                                modify = 0;
                                out(i) = intver;
                            end
                        end
                        jj = jj + 1;
                        if (jj >= 10)
                            disp(['Search bounds for var: ',...
                                num2str(i),...
                                ' might be too tight. Cannot fix an ',...
                                'int type. Breaking....'])
                            int_c_rest = 1/2;
                            %coef_restitution = 9/jj;
                        end
                    end
                else % If not out of bounds, just use rounded version
                    out(i) = intver;
                end
                
                % 13-19 just need to be within bounds, but can stay float
            else
                a = search_bounds(1,i); b = search_bounds(2,i);
                if ((in(i) < a) || (in(i) > b)) % If out of bounds,
                    modify = 1;
                    coef_restitution = 1;
                    jj = 0;
                    while modify == 1
                        if in(i) < a
                            amount_out = a - in(i);
                            in(i) = b - amount_out*coef_restitution;
                            if ((in(i) >= a) && (in(i) <= b))
                                modify = 0;
                                out(i) = in(i);
                            end
                            
                        elseif in(i) > b
                            amount_out = in(i) - b;
                            in(i) = a + amount_out*coef_restitution;
                            if ((in(i)>= a) && (in(i) <= b))
                                modify = 0;
                                out(i) = in(i);
                            end
                        end
                        jj = jj + 1;
                        coef_restitution = instabilityMitigation(jj);
                    end
                else % If not out of bounds, just use rounded version
                    out(i) = in(i);
                end
            end
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HARD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Hard'
        
        for i = 1:length(in)
            % Process and make feasible
            if param_type(i) == 1
                % 1-12 must be type int AND within range (a,b)
                a = search_bounds(1,i); b = search_bounds(2,i);
                intver = round(in(i));
                if ((intver < a) || (intver > b)) % If out of bounds,
                    modify = 1;
                    jj = 0;
                    int_c_rest = 1;
                    while modify == 1
                        if intver < a
                            increaseby = a-intver;
                            intver = intver + 2*increaseby*int_c_rest;
                            if ((intver >= a) && (intver <= b))
                                modify = 0;
                                out(i) = intver;
                            end
                        elseif intver > b
                            decreaseby = b-intver;
                            intver = intver + 2*decreaseby*int_c_rest;
                            if ((intver >= a) && (intver <= b))
                                modify = 0;
                                out(i) = intver;
                            end
                        end
                        jj = jj + 1;
                        if (jj >= 10)
                            disp(['Search bounds for var: ',...
                                num2str(i),...
                                ' might be too tight. Cannot fix an ',...
                                'int type. Breaking....'])
                            int_c_rest = 1/2;
                            %coef_restitution = 9/jj;
                        end
                    end
                else
                    out(i) = intver;
                end
                % 13-19 just need to be within bounds, but can stay float
            elseif param_type(i) == 3
                % 1-12 must be type int AND within range (a,b)
                a = search_bounds(1,i); b = search_bounds(2,i);
                intver = round_oddint(in(i));
                if ((intver < a) || (intver > b)) % If out of bounds,
                    modify = 1;
                    jj = 0;
                    int_c_rest = 1;
                    while modify == 1
                        if intver < a
                            increaseby = a-intver;
                            intver = intver + 2*increaseby*int_c_rest;
                            if ((intver >= a) && (intver <= b))
                                modify = 0;
                                out(i) = intver;
                            end
                        elseif intver > b
                            decreaseby = b-intver;
                            intver = intver + 2*decreaseby*int_c_rest;
                            if ((intver >= a) && (intver <= b))
                                modify = 0;
                                out(i) = intver;
                            end
                        end
                        jj = jj + 1;
                        if (jj >= 10)
                            disp(['Search bounds for var: ',...
                                num2str(i),...
                                ' might be too tight. Cannot fix an ',...
                                'int type. Breaking....'])
                            int_c_rest = 1/2;
                            %coef_restitution = 9/jj;
                        end
                    end
                else
                    out(i) = intver;
                end
                % 13-19 just need to be within bounds, but can stay float
            else
                a = search_bounds(1,i); b = search_bounds(2,i);
                
                if ((in(i) < a) || (in(i) > b)) % If out of bounds,
                    modify = 1;
                    jj = 0;
                    coef_restitution = 1;
                    while modify == 1
                        if in(i) < a
                            increaseby = a-in(i);
                            in(i) = in(i) + 2*increaseby*coef_restitution;
                            if ((in(i) >= a) && (in(i) <= b))
                                modify = 0;
                                out(i) = in(i);
                            end
                        elseif in(i) > b
                            decreaseby = b-in(i);
                            in(i) = in(i) + 2*decreaseby*coef_restitution;
                            if ((in(i) >= a) && (in(i) <= b))
                                modify = 0;
                                out(i) = in(i);
                            end
                        end
                        jj = jj + 1;
                        coef_restitution = instabilityMitigation(jj);
                    end
                else
                    out(i) = in(i);
                end
                
            end
        end
        
        
end

    function [coef_rest] = instabilityMitigation(iterator)
        if ((iterator >= 5) && (iterator < 10))
            coef_rest = 4/iterator;
        elseif iterator >= 10
            coef_rest = 4/iterator^2;
            disp(['Search bounds for var: ', num2str(i),...
                ' might be too tight. While loop has tried to modify >10 times.'])
        else
            coef_rest = 1;
            
        end
    end

    function [out] = round_oddint(x)
        out = 2*floor(x/2)+1;
    end

end