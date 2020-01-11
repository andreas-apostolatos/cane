function create_control_polygon(CP)
%% Function documentation
%
% Plots control points and control polygon for the given set of Control
% Points
%
%   Input :
%      CP : Set of Control Point coordinates
%
%  Output : Graphics
%
%% Function main body

for k=1:length(CP(1,:,1))-1
    for l=1:length(CP(:,1,1))-1
        plot3(CP(l,k:k+1,1),CP(l,k:k+1,2),CP(l,k:k+1,3),'--or');
        plot3(CP(l:l+1,k,1),CP(l:l+1,k,2),CP(l:l+1,k,3),'--or');
    end
    % Update counter
    l=l+1;
    plot3(CP(l,k:k+1,1),CP(l,k:k+1,2),CP(l,k:k+1,3),'--or');
end
% Update counter
k=k+1;
for l =1:length(CP(:,1,1))-1
    plot3(CP(l:l+1,k,1),CP(l:l+1,k,2),CP(l:l+1,k,3),'--or');
end

end