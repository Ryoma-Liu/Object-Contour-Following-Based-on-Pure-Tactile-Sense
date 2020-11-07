%find the contact location (inputs raw values, outputs location, finger surface normal, and force magnitude)
function [loc, normal,tforce,forcemag]=process_values(values)
 RADIUS =18;% //radius of fingertip contact cylinder in mm
 LENGTH =12; %  //(estimated!) length from zero point to center of sphere
 ACTUAL_COUNTS_PER_FORCE= 1;
 ACTUAL_COUNTS_PER_TORQUE= 1;
    normal=[0 0 0];  tforce=[0 0 0];
  Fx=values(1)/ACTUAL_COUNTS_PER_FORCE;
  Fy=values(2)/ACTUAL_COUNTS_PER_FORCE;
  Fz=values(3)/ACTUAL_COUNTS_PER_FORCE;
  Tx=values(4)/ACTUAL_COUNTS_PER_TORQUE;
  Ty=values(5)/ACTUAL_COUNTS_PER_TORQUE;
  Tz=values(6)/ACTUAL_COUNTS_PER_TORQUE;
  %fprintf("%d %d %d %d %d %d\n",values(1),values(2),values(3),values(4),values(5),values(6));
  fprintf("%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n",Fx,Fy,Fx,Tx,Ty,Tz,ACTUAL_COUNTS_PER_FORCE,ACTUAL_COUNTS_PER_TORQUE);
 x=0;y=0;z=0;
  %find contact point and normal to finger surface
  theta=atan2(Fy,Fx);
  cosT=cos(theta);
  sinT=sin(theta);
  if (theta<0) 
    angle=360+theta*180/pi;
  else 
    angle=theta*180/pi;
  end
	forcemag = sqrt(Fx*Fx + Fy*Fy + Fz*Fz);
    zangle=0;
	if(forcemag > 0.0001)
        %angle from z-axis
		zangle = atan2(sqrt(Fx*Fx + Fy*Fy), -Fz);
    %contact on the cyl
	if(zangle >= 1.57)
			x = -RADIUS*cosT;
			y = -RADIUS*sinT;
			normal(1) = -cosT;
			normal(2) = -sinT;
			normal(3) = 0;
			
			 TxminusFz = Tx - Fz*RADIUS*cosT;
			 TyminusFz = Ty + Fz*RADIUS*sinT;
			z = sqrt(TxminusFz*TxminusFz + TyminusFz*TyminusFz) / sqrt(Fx*Fx + Fy*Fy);
			%z should be < 0
			if((Fx>0 && Fy>0 && Tx>0 && Ty<0) ||...
				 (Fx<0 && Fy<0 && Tx<0 && Ty>0) ||...
				 (Fx<0 && Fy>0 && Tx>0 && Ty>0) ||...
				 (Fx>0 && Fy<0 && Tx<0 && Ty<0))
                  z = -z;			
		%contact on the sphere
            else
                x = -RADIUS*cosT*sin(zangle);
                y = -RADIUS*sinT*sin(zangle);
                z = LENGTH + RADIUS*cos(zangle);
                normal(1) = -cosT*sin(zangle);
                normal(2) = -sinT*sin(zangle);
                normal(3) = cos(zangle);
            end
    end
	else 
		x=0;
		y=0;
        z=0;
        end
	loc(1) = x;
	loc(2) = y;
	loc(3) = z; 
    tforce(1)=Fx-normal(1)*forcemag;
    tforce(2)=Fy-normal(2)*forcemag;
    tforce(3)=Fz-normal(3)*forcemag;
end
