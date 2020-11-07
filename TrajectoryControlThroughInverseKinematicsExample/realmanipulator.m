robot = importrobot('demo.urdf');
showdetails(robot);
robot.NumBodies ;
robot.Bodies;

body3 = getBody(robot,'link_1')

body3.Mass = 5;
body3.CenterOfMass = [0.5 0 0];
body3.Inertia = [0.167 0.001 0.167 0 0 0];

body3 = getBody(robot,'link_2')

body3.Mass = 4;
body3.CenterOfMass = [0.5 0 0];
body3.Inertia = [0.167 0.001 0.167 0 0 0];

body3 = getBody(robot,'link_5')

body3.Mass = 3;
body3.CenterOfMass = [0.5 0 0];
body3.Inertia = [0.167 0.001 0.167 0 0 0];

body3 = getBody(robot,'link_6')

body3.Mass = 2;
body3.CenterOfMass = [0.5 0 0]; 
body3.Inertia = [0.167 0.001 0.167 0 0 0];


body3 = getBody(robot,'link_7')

body3.Mass = 1;
body3.CenterOfMass = [0.5 0 0];
body3.Inertia = [0.167 0.001 0.167 0 0 0];
showdetails(robot);
body3 = getBody(robot,'link_6')
body3 = getBody(robot,'link_7')