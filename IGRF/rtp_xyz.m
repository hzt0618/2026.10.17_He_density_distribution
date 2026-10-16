function [x,y,z]=rtp_xyz(r,theta,phi)
if theta<0 || theta>pi
    disp('Error: theta beyond the boundary')
    return
end

if phi<0 || phi>2*pi
    disp('Error: phi beyond the boundary')
end

x=r*sin(theta)*cos(phi);
y=r*sin(theta)*sin(phi);
z=r*cos(theta);

end