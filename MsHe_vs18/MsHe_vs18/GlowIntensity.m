function It=GlowIntensity(ze,den_e,sza,whre)

%% 把不同高度的散射强度相加得到地面的被动观测结果
% whre 标记参数，0：考虑共振荧光；1：不考虑共振荧光

nze=length(ze);

% calculate interval of ze
dze=zeros(1,nze);
for i=1:nze-1
    dze(i)=ze(i+1)-ze(i);   % km
end
dze(nze)=dze(nze-1);   % km

if whre==0
    %% consider resonance
    % calculate shadow height
    Re=6.371393*10^(6);    %m地球半径
    h_shadow=(Re/1000+30)/sin(pi-sza)-Re/1000;
    
    % calculate density integral along height, consider resonance
    It=0;
    density=zeros(1,nze);
    for i=1:length(ze)
        if (ze(i)+dze(i)/2)<h_shadow
            density(i)=0;
        elseif (ze(i)-dze(i)/2)>h_shadow
            density(i)=den_e(i);
        else
            density(i)=den_e(i)*(ze(i)+dze(i)/2-h_shadow)/dze(i);
        end
        It=It+dze(i)*10^5*density(i);   % photons s-1 cm-2
    end
    
    % Rayleighs: 10^6 phontons s-1 cm-2
    g_photos=16.8;      % s^-1
    sigma_t=g_photos*scafac1(sza); %/4/pi;
    It=It*sigma_t*10^(-6); % Rayleighs
elseif whre==1
    %% no consider resonance
    It=0;
    % 每平方厘米柱上的亚稳态氦原子数
    for i=1:nze
        It=It+den_e(i)*dze(i)*10^5;    % cm-2
    end
    % Rayleighs: 10^6 phontons s-1 cm-2
    g_photos=16.8;      % s^-1
    sigma_t=g_photos*scafac1(sza); %/4/pi;
    It=It*sigma_t*10^(-6); % Rayleighs
end

end