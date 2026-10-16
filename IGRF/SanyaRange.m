%% 配套MagLine,MagLinePlot,给出三亚的磁力线和观测范围
figure(1);

%% 给出观测范围的角度
ZA=1*pi/180;

r3=r0+1000;
theta3=theta0;
lambda3=lambda0;
[x3,y3,z3]=rtp_xyz(r3,theta3,lambda3);
R=1000*tan(ZA);
n=linspace(0,2*pi,100);
m=atan(-z3./(x3*cos(n)+y3*sin(n)))-(sign(-z3./(x3*cos(n)+y3*sin(n)))-1)/2*pi;
x4=x3+R*sin(m).*cos(n);
y4=y3+R*sin(m).*sin(n);
z4=z3+R*cos(m);

plot3(x4,y4,z4);
hold on;

len=length(x4);
[x0,y0,z0]=rtp_xyz(r0,theta0,lambda0);
num=100;
for i=1:len
    x5=(x0*(0:num)+x4(i)*(num:-1:0))/num;
    y5=(y0*(0:num)+y4(i)*(num:-1:0))/num;
    z5=(z0*(0:num)+z4(i)*(num:-1:0))/num;
    plot3(x5,y5,z5)
    hold on
end


[X,Y,Z] = sphere(100);
r = 6371.2;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;

s=surf(X2,Y2,Z2);
s.EdgeColor = 'none';
axis equal
hold on