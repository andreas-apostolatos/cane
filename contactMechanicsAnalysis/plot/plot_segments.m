function [] = plot_segments(segments)
%PLOTSEGMENTS Plot segments on current graph
%   Plot segments on current graph
 
%% 1. Add the segments

% create marks coordinates
if(ndims(segments)>2) % ie more than one segments possible
    for i=1:size(segments,1)
        contWall = createWall(squeeze(segments(i,:,:)));
        hold on;
        if ~isempty(segments)
            plot(squeeze(segments(i,:,1)),squeeze(segments(i,:,2)),'Linewidth',2,'Color','black');
        end
        for k =1:length(contWall.xw(:,1))
            plot(contWall.xw(k,:),contWall.yw(k,:),'LineWidth',1,'Color','black');
        end
        hold off;
    end
else
    contWall = createWall(segments);
    hold on;
    if ~isempty(segments)
        plot(segments(:,1),segments(:,2),'Linewidth',2,'Color','black');
    end
    for k =1:length(contWall.xw(:,1))
        plot(contWall.xw(k,:),contWall.yw(k,:),'LineWidth',1,'Color','black');
    end
    hold off;
end

end

