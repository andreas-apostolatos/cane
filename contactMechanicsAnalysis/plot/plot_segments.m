function [] = plot_segments(segments)
%   Plot segments on current graph
 
%% 
% check if segmets exist
if ~isempty(segments.points)
    % loop through the number of segments
    for i=1:size(segments.points,3)
        % create contour
        contourWall = createWall(segments.points(:,:,i));

        hold on;
        % plot a line form start to end coordiante x0 to x1
        startPoint = segments.points(:,1,i);
        endPoint = segments.points(:,2,i);
        plot(startPoint,endPoint,'-','Linewidth',2,'Color','black');
        
        %plot(startPoint,endPoint,':','Linewidth',1,'Color','black');
       
        
        % plot contours below the line
%         for k =1:length(contourWall.xw(:,1))
%             plot(contourWall.xw(k,:),contourWall.yw(k,:),'LineWidth',1,'Color','black');
%         end
        hold off;
    end
end

end

