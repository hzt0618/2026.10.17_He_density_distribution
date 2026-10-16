function [r,theta,phi]=xyz_rtp(x,y,z)

r=sqrt(x^2+y^2+z^2);
theta=acos(z/r);
phi=mod(atan2(y,x),2*pi);

end
