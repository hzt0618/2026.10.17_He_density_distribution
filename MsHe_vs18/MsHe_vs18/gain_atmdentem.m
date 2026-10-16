function [Tem,He,O1,N2,O2,den]=gain_atmdentem(z,glat,glong,Yr,DoY,ut,f107a,f107,ap)
%Implement mathematical representation of 2001 United States 
%Naval Research Laboratory Mass Spectrometer 
%and Incoherent Scatter Radar Exosphere
%输入z，单位为米，0-1,000,000
[T,DEN]=atmosnrlmsise00(z*1000,glat,glong,Yr,DoY,ut,f107a,f107,ap);
%{
T，Array of N-by-2 values of temperature, in kelvin. 
The first column is exospheric temperature, in kelvin. 
The second column is temperature at altitude, in kelvin.
DEN,rho,An array of N-by-9 values of densities (kg/m3 or 1/m3) 
in selected density units. 
The column order is:
Density of He, in 1/m3
Density of O, in 1/m3
Density of N2, in 1/m3
Density of O2, in 1/m3
Density of Ar, in 1/m3
Total mass density, in kg/m3
Density of H, in 1/m3
Density of N, in 1/m3
Anomalous oxygen number density, in 1/m3
density(6), total mass density, 
is the sum of the mass densities of He, O, N2, O2, Ar, H, and N. 
Optionally, density(6) can include the mass density of 
anomalous oxygen making density(6), 
the effective total mass density for drag. 
%}
len=size(T,1);
Tem=reshape(T(:,2),[1,len]);    %大气温度,单位开尔文
He=reshape(DEN(:,1),[1,len])* 10.^(-6);  %cm-3;      %He的密度,单位1/m^3
O1=reshape(DEN(:,2),[1,len])* 10.^(-6);  %cm-3;          %O_1的密度,单位1/m^3
N2=reshape(DEN(:,3),[1,len])* 10.^(-6);  %cm-3;         %N_2的密度,单位1/m^3
O2=reshape(DEN(:,4),[1,len])* 10.^(-6);  %cm-3;         %O_2的密度,单位1/m^3
den=reshape(DEN(:,6),[1,len])* 10.^(-6);  %cm-3;            %H大气的密度,单位1/m^3
end