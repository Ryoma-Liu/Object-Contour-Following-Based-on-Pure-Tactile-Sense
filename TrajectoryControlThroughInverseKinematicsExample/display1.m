data1 = poseData1;
data2 = poseData2;
% plot3(targetPoint(1), targetPoint(2), targetPoint(3), 'go', 'MarkerSize', 5);
% hold on;
figure(1);
for i = 1:100
   
   u1 = data1(1:3,4,i);
plot3(u1(1), u1(2), u1(3), 'ro', 'MarkerSize', 5);
hold on
   u2 = data2(1:3,4,i);
plot3(u2(1), u2(2), u2(3), 'go', 'MarkerSize', 5);

end