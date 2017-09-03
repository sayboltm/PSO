function [ null ] = paretoPlotter( g, p, field1, field2, pltbaseresponse, pltbestresponse, pltnearbest )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Extrapoints is a numextrapoints x fieldn object
% so to plot a set of points for all particles, extrapoints would be a
% pop_size x fields object
pop_size = length(p.properties);
figure()
hold on

% If desired, plot base response points first
if pltbaseresponse == 1
    for i = 1:pop_size;
        bapoints = plot(p.properties(i).base_response.(field1), p.properties(i).base_response.(field2),'rx');
%         plot(p.properties(i).(field1), p.properties(i).(field2), 'rx')
        %plot(extrapoints(i).(field1), extrapoints(i).(field2), 'rx')
    end
end

%Plot initial and best points for shits and giggs
for i = 1:length(p.storage.trial);
    
    if ((isempty(p.storage.trial(i).Objectives.(field1))) || (isempty(p.storage.trial(i).Objectives.(field2))))
        continue
    else
        plot(p.storage.trial(i).Objectives.(field1), p.storage.trial(i).Objectives.(field2),'rx');
        %plot(p.storage.trial(i).Objectives.(field1), p.storage.trial(i).Objectives.(field2),'mx');
    end
    %disp(['Shit: ', num2str(i)])
end

% Plot points near best global fitness (fake pareto points)
if pltnearbest == 1
    for i = 1:pop_size;
        if p.properties(i).Pareto == 1;
            %j = j + 1;
            
            %fppoints = plot(p.properties(i).(field1), p.properties(i).(field2),'mx');
            fppoints = plot(p.properties(i).(field1), p.properties(i).(field2),'rx');
            
            %         for k = 1:length(sublam_fields)
            %             Sublaminate(j).(sublam_fields{k}) = p.properties(i).(sublam_fields{k});
            %         end
        end
    end
end

% Plot Pareto points
for i = 1:length(p.storage.trial)
    if p.storage.trial(i).paretoactual == 1;
        ppoints = plot(p.storage.trial(i).Objectives.(field1), p.storage.trial(i).Objectives.(field2), 'o');
        %disp(['trial: ', num2str(i), ', Objective.', field1, ' = ', num2str(p.storage.trial(i).Objectives.(field1))])
    end
end

% Plot globalbest if desired
if pltbestresponse == 1
    bpoint = plot(g.bestproperties.(field1), g.bestproperties.(field2), 'g.');
end



hold off
xlabel(field1)
ylabel(field2)

plt_response_options = strcat(num2str(pltbaseresponse), num2str(pltbestresponse), num2str(pltnearbest));
switch plt_response_options
    case '000'
        title(['Pareto Front of: ',field1, ' vs ', field2])
        legend('Pareto Points')
    case '100'
        title(['Pareto Front of: ',field1, ' vs ', field2, ', with base points'])
        legend([ppoints bapoints], {'Pareto Points', 'Base Response'})
    case '010'
        title(['Pareto Front of: ',field1, ' vs ', field2, ', with "best" point']) 
        legend([ppoints bpoint], {'Pareto Points', '"Best" Point'})
    case '110'
        title(['Pareto Front of: ',field1, ' vs ', field2, ', with base and "best" points'])
        %title(['Pareto Front of: ',field1, ' vs ', field2, ', with ALL POINTS'])
        legend([ppoints bapoints bpoint], {'Pareto Points', 'Base Responses', '"Best" Point'})
    case '111'
        title(['Pareto Front of: ',field1, ' vs ', field2, ', with point progression and "best" point']) 
        %legend([ppoints bapoints bpoint fppoints], {'Pareto Points', 'Point Progression' , '"Best" Point', 'Near Best'})
        legend([ppoints bapoints bpoint], {'Pareto Points', 'Point Progression' , '"Best" Point'})
	case '011'
		title(['Pareto Front of: ',field1, ' vs ', field2, ', with "best" point, and near best']) 
        legend([ppoints bpoint fppoints], {'Pareto Points', '"Best" Point', 'Near Best'})
    otherwise
        disp(['Option not yet supported. Plotting Pereto points only.'])
        title(['Pareto Front of: ',field1, ' vs ', field2])
        legend('Pareto Points')
end


end

