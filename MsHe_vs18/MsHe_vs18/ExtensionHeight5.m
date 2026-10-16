function [den_e,ze2]=ExtensionHeight5(ze,interval,sza,Tem,result)

%% 把亚稳态氦的结果向更高高度外推
    usehu_e=3000;
    usehd_e=100;
    interval_e=interval;
    nze=length(ze);
    nze2=(usehu_e-usehd_e)/interval_e+1;
    ze2=linspace(usehd_e,usehu_e,nze2);  %km
    den_e=zeros(1,nze2);
    den_e(1:nze)=result;
    [Hn_e,g_e,a_e,gstar_e]=he_Hn2(ze2,Tem(end),sza);
    for i=nze+1:nze2
        den_e(i)=result(nze).*exp(-(ze2(i)-ze(end))*10^5/Hn_e(i));
    end
end