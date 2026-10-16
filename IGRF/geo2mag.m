% 将大地经纬度转化为地磁经纬度
% 磁极点取：经度72.2°W，纬度80.0°N（IGRF，2011）
% 输入及输出值为弧度制单位（rad）
% 磁极点经纬度

b=80.0*pi/180;
l=-72.2*pi/180;

Bg = theta1;
Lg = lambda1;
% 以下公式基于刘长建，2011

sb=sin(Bg)*sin(b)+cos(Bg)*cos(b)*cos(Lg-l);

Bm=asin(sb);

sl=cos(Bg)*sin(Lg-l)/cos(Bm);

cl=(sin(b)*sb-sin(Bg))/(cos(b)*cos(Bm));

Lm=atan2(sl,cl);

qqq = 90 - Bm*180/pi;
aaa = Lm*180/pi;
disp(qqq);
disp(aaa);