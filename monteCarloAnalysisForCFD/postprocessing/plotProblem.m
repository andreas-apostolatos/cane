function [] = plotProblem (problemStruct,varargin)

    nDoF = 3;
    
    % Init
    options = containers.Map(                                           ...
            {'fldMsh', 'homDBC', 'inhomDBC'},                           ...
            { true,    true,      true  });
    if ~isempty(varargin)
        % Set values to false
        key = keys(options);
        for k=1:length(options)
            options(key{k}) = false;
        end
        % Set specific values to true
        for k=1:length(varargin)
            if isKey(options,varargin{k})
                options(varargin{k}) = true;
            end
        end
    end
    drawn = [];
    pltLegend = [];
    
    % Plot all elements
    figure(1)
    hold on
    if options('fldMsh')
        for k=1:length(problemStruct.fldMsh.elements)
            nodeList = problemStruct.fldMsh.elements(k,:);
            nodeList = [nodeList,nodeList(1)];
            for j=2:length(nodeList)
                % Sorted pair
                if nodeList(j-1)>nodeList(j)
                    nodePair = nodeList(j-1:j);
                else
                    nodePair = [nodeList(j),nodeList(j-1)];
                end
                % Check if pair was already drawn
                if isempty(drawn)
                    drawn(1,1:2) = nodePair;
                    plot(                                                   ...
                        problemStruct.fldMsh.nodes( nodePair,1 ),           ...
                        problemStruct.fldMsh.nodes( nodePair,2 ),           ...
                        'b');
                else
                    if sum( drawn(:,1)==nodePair(1) & drawn(:,2)==nodePair(2) ) == 0
                        drawn(end+1,1:2) = nodePair;
                        plot(                                               ...
                            problemStruct.fldMsh.nodes( nodePair,1 ),       ...
                            problemStruct.fldMsh.nodes( nodePair,2 ),       ...
                            'b');
                    end
                end
            end

        end
    end
    
    % Plot nodes with homDBC (boundary nodes -> black, bodies -> red)
    if options('homDBC')
        
        bodyNodes   = getBodyNodes(problemStruct.fldMsh,problemStruct.homDBC);

        for k=1:length(problemStruct.homDBC)
            index = ceil(problemStruct.homDBC(k)/nDoF);
            % Check if point is on the boundary
            if ismember(index,bodyNodes)
                marker = 'r+';
            else
                marker = 'ko';
            end

            plot(                                                       ...
                problemStruct.fldMsh.nodes(index, 1),                   ...
                problemStruct.fldMsh.nodes(index, 2),                   ...
                marker                                                  ...
            )
        end
    end
    
    % Plot nodes with inhomDBC
    if options('inhomDBC')
        nrm = (     0.1/max(problemStruct.valuesInhomDBC)*              ...
                    max( problemStruct.fldMsh.nodes(:,1))               ...
                    -                                                   ...
                    min(problemStruct.fldMsh.nodes(:,2))                ...
               );
        for k=1:length(problemStruct.inhomDBC)
            index = ceil(problemStruct.inhomDBC(k)/nDoF);
            x = problemStruct.fldMsh.nodes(index, 1);
            y = problemStruct.fldMsh.nodes(index, 2);
            plot(                                                       ...
                [x,x+nrm*problemStruct.valuesInhomDBC(k)],              ...
                [y,y],                                                  ...
                'g'                                                     ...
            )
        end
    end
    
    
    % Legend workaround
    pltLegend(end+1) = plot(NaN,NaN,'b');
    pltLegend(end+1) = plot(NaN,NaN,'ko');
    pltLegend(end+1) = plot(NaN,NaN,'r+');
    pltLegend(end+1) = plot(NaN,NaN,'g');
    legend(pltLegend, 'Elements','Nodes with homDBC on boundaries', 'Nodes with homDBC on bodies','inhomDBC');
    
    axis equal
    hold off

end