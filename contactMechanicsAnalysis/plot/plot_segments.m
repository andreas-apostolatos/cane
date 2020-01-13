function [] = plot_segments(segments)
%PLOTSEGMENTS Plot segments on current graph
%   Plot segments on current graph
 
%% 1. Add the segments

% create marks coordinates
if(size(segments,3)>1) % if more than one segments exist
    for i=1:size(segments,3)
        contWall = createWall(squeeze(segments(:,:,i)));
        hold on;
        if ~isempty(segments)
            plot(squeeze(segments(:,1,i)),squeeze(segments(:,2,i)),'Linewidth',2,'Color','black');
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

