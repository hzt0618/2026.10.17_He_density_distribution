%% 是一个函数，输入地心坐标，可以得到此处在地心坐标系下的磁场方向矢量。
% 考虑地球为球形

function [B_r,B_theta,B_lambda]=IGRF_vs2(r,theta,lambda)
persistent gmn hmn nmax Re

if isempty(gmn)
    [gmn, hmn, nmax, Re] = load_igrf_coeffs();
end

B_r=0;
B_theta=0;
B_lambda=0;

for n=1:nmax
    for m=0:n
        B_r=B_r+((-n-1)*(Re./r).^(n+2)).*(gmn(n,m+1).*cos(m*lambda)+...
            hmn(n,m+1).*sin(m*lambda)).*Pnm(n,m,theta);
    end
end

for n=1:nmax
    for m=0:n
        B_theta=B_theta+(Re./r).^(n+2).*(gmn(n,m+1).*cos(m*lambda)+...
            hmn(n,m+1).*sin(m*lambda)).*dPnm_dsita(n,m,theta);
    end
end

for n=1:nmax
    for m=0:n
        B_lambda=B_lambda-(Re./r).^(n+2).*(m*gmn(n,m+1).*sin(m*lambda)-...
            m*hmn(n,m+1).*cos(m*lambda))./sin(theta).*Pnm(n,m,theta);
    end
end

end

function [gmn, hmn, nmax, Re] = load_igrf_coeffs()
filename='.\igrf13coeffs.txt';
fileID=fopen(filename,'r');
if fileID == -1
    error('IGRF coefficient file not found: %s', filename)
end

year=2020;
ilen=round(-(year-2020)/5)+1;
data=textscan(fileID,['%s %n %n',repmat('%*f',[1,25-ilen]) ...
    ,'%f',repmat('%*f',[1,ilen])],'headerlines',4);
fclose(fileID);

Re=6371.2;
nmax=max(data{2});
mmax=max(data{3})+1;
hmn=zeros(nmax,mmax);
gmn=zeros(nmax,mmax);

for i=1:size(data{1},1)
    if data{1}{i}=='h'
        hmn(data{2}(i),data{3}(i)+1)=data{4}(i);
    end
    if data{1}{i}=='g'
        gmn(data{2}(i),data{3}(i)+1)=data{4}(i);
    end
end

end


