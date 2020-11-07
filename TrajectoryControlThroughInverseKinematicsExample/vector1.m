point = zeros(100,3);
fn = zeros(100,3);
ft = zeros(100,3);

for i = 1 : 100
    point(i,1) = 0;
%     point(i,2) = -1+0.02*i;
    point(i,2) = -1.57+0.0314*i;
%     point(i,3) = sqrt(1-point(i,2)*point(i,2));
    point(i,3) = -cos(point(i,2));
    fn(i,1) = 0;
    fn(i,2) = 0-point(i,2);
    fn(i,3) = 0-point(i,3);
    ft(i,2:3) = fn(i,2:3)*[0  1
                          -1 0];

end
fn
ft