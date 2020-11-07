function LMDemo7
M = csvread('1-6-3.csv');
M=M(1:3160,:);
for j=1:3160
dz=2;r=18;
 [loc, normal(j,1:3),tforce(j,1:3),forcemag(j,1:3)]=process_values(M(j,:));
 x(j)=loc(1);
 y(j)=loc(2);
 z(j)=loc(3);
end
close all
figure(1)
plot3(x,y,z,'.');
figure(2)
plot(normal);
figure(3)
plot(tforce); 
figure(4)
plot(forcemag)
end
