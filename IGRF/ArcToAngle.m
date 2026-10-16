%% 弧度转角度的函数
function [deg1,min,sec]=ArcToAngle(arc)
arc0=arc*180/pi;
deg1=fix(arc0);     %向0的方向取整
min=fix((arc0-deg1)*60);
sec=((arc0-deg1)*60-min)*60;
end