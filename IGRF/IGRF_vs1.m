%% 按照球坐标均分全球，得到的网格点处的磁场数据。

clear;
close all;
filename='.\igrf13coeffs.txt';
fileID=fopen(filename,'r');
year=2020;
ilen=round(-(year-2020)/5)+1;
data=textscan(fileID,['%s %n %n',repmat('%*f',[1,25-ilen]) ...
    ,'%f',repmat('%*f',[1,ilen])],'headerlines',4);
fclose(fileID);
Re=6371.2;
rnum=100;
L=5;
pi=3.1415926535897932384626433832795028841;
snum=rnum;
fnum=2*rnum;
gui='norm';

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

pi_snum=3*pi/(4*snum);
pi_fnum=3*pi/(4*fnum);

%磁场分量
B_r     =  zeros(snum,fnum);
B_sita  =  zeros(snum,fnum);
B_faif  =  zeros(snum,fnum);

%大地坐标
h       =  0;
lamd    =  linspace(0+pi_fnum,2*pi-pi_fnum,fnum);
sita0   =  linspace(0+pi_snum,pi-pi_snum,snum);

a=6378.16;
b=6356.775;

%球坐标网格、地心坐标
r       =  sqrt(h*(h+2*sqrt(a^2*sita0.^2))+(a^4*sin(sita0).^2+ ...
               b^4*cos(sita0).^2)./(a^2*sin(sita0).^2+b^2*cos(sita0).^2));
sita    =  sita0-acos((h+sqrt(a^2*sin(sita0).^2+b^2*cos(sita0).^2))./r);

r       =  reshape(r,[snum,1]);
r1      =  repmat(r,[1,fnum]);

sita       =  reshape(sita,[snum,1]);
sita1   =  repmat(sita,[1,fnum]);

faif    =  lamd;
faif1   =  repmat(faif,[snum,1]);

for n=1:nmax
    for m=0:n
        B_r=B_r+((-n-1)*(Re./r1).^(n+2)).*(gmn(n,m+1).*cos(m*faif1)+...
            hmn(n,m+1).*sin(m*faif1)).*Pnm(n,m,sita1);
    end
end

for n=1:nmax
    for m=0:n
        B_sita=B_sita+(Re./r1).^(n+2).*(gmn(n,m+1).*cos(m*faif1)+...
            hmn(n,m+1).*sin(m*faif1)).*dPnm_dsita(n,m,sita1);
    end
end

for n=1:nmax
    for m=0:n
        B_faif=B_faif+(Re./r1).^(n+2).*(-m*gmn(n,m+1).*sin(m*faif1)+...
            m*hmn(n,m+1).*cos(m*faif1)).*Pnm(n,m,sita1);
    end
end

B_r=-B_r;
B_sita=-B_sita;
B_faif=-B_faif;
