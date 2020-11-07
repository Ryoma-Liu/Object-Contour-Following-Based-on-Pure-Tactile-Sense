% Import the manipulator as a rigidBodyTree Object
sawyer = importrobot('sawyer.urdf');
sawyer.DataFormat = 'column';

% Define end-effector body name
eeName = 'right_hand';

% Define the number of joints in the manipulator
numJoints = 8;

% Visualize the manipulator
show(sawyer);
xlim([-1.00 1.00])
ylim([-1.00 1.00]);
zlim([-1.02 0.98]);
view([128.88 10.45]);


I = imread('coins.png');
bwBoundaries = imread('coinBoundaries.png');

figure
subplot(1,2,1)
imshow(I,'Border','tight')
title('Original Image')

subplot(1,2,2)
imshow(bwBoundaries,'Border','tight')
title('Processed Image with Boundary Detection')


load boundaryData.mat boundaries
whos boundaries