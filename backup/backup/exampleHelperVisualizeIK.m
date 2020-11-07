function exampleHelperVisualizeIK(I, jointData, poseData, sawyer, imageOrigin, scale)
%exampleHelperVisualizeIK Function to visualize outputs of Inverse Kinematics example

% Copyright 2018 The MathWorks, Inc.

% Remove unnecessary meshes for faster visualization
clearMeshes(sawyer);

% Data for mapping image
[m,n] = size(I);

[X,Y] = meshgrid(0:m,0:n);
X = imageOrigin(1) + X*scale;
Y = imageOrigin(2) + Y*scale;

Z = zeros(size(X));
Z = Z + imageOrigin(3);

% Initialize end effector visualization
f = figure;
set(f, 'Visible','on');
show(sawyer, jointData(1,:)');
hold on

% Initialize end effector plot position
p = plot3(0,0,0,'.');
warp(X,Y,Z,I');

% Change view angle and axis
view(65,45)
axis([-.25 1 -.25 .75 0 0.75])

% Iterate through the outputs at 10-sample intervals to visualize the results
for j = 1:10:length(jointData)
    % Display manipulator model
    show(sawyer,jointData(j,:)', 'Frames', 'off', 'PreservePlot', false);
    
    % Get end effector position from homoegenous transform output
    pos = poseData(1:3,4,j);
    
    % Update end effector position for plot
    p.XData = [p.XData pos(1)];
    p.YData = [p.YData pos(2)];
    p.ZData = [p.ZData pos(3)];
    
    % Update figure
    drawnow
end

end
