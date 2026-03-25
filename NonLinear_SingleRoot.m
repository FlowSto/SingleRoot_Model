%% Sand no root hairs
clear all
clc
%%
hb= linspace(-100,-10000,10); % bulk soil water potential


xx= [29.06  3      6.0e-4 % LOAM
    3.402   3.6    0.047535 %SAND
    ];
%%
Rroot=1.2261e-11;%hPa cm^-3 s root resistance;
r= 0.03;%+0.0245; % 0.03 cm root radius in cm
r2= 1; %bulk+root+rhizo soil radius in cm
L = 1;% root segment length

hroot = 15000; %soil:root interface potential
%%
for k=2 % 1 = loam; 2 = sand
    for i = 1:10 % iteratons of bulk soil water potential
        h0=-xx(k,1);
        k0=xx(k,3);
        tau=xx(k,2);
        
        Ks_nonlinear(i) = (2.*pi.*r.*L.* k0/(1-tau)/(h0^(-tau))./(r./2- r.*r2^2.*(log(r2)-log(r))./(r2^2-r.^2) ).*(-hb(i)^(1 - tau) - (-hroot)^(1 - tau)))./abs(hb(i)-hroot);
    end
end

%% Loam no root hairs
clear all
clc

%%
hb= linspace(-100,-10000,10); % bulk soil water potential


xx= [29.06  3      6.0e-4 % LOAM
    3.402   3.6    0.047535 %SAND
    ];
%%
Rroot=1.2261e-11;%hPa cm^-3 s;
r= 0.03;%+0.0245; % 0.03 cm root radius in cm
r2= 1; %bulk+root+rhizo soil radius in cm
L = 1;% root segment length

hroot = 15000; %soil:root interface potential
%%
for k=1 % 1 = loam; 2 = sand
    for i = 1:10 % iteratons of bulk soil water potential
        h0=-xx(k,1);
        k0=xx(k,3);
        tau=xx(k,2);
        
        
        Ks_nonlinearL(i) = (2.*pi.*r.*L.* k0/(1-tau)/(h0^(-tau))./(r./2- r.*r2^2.*(log(r2)-log(r))./(r2^2-r.^2) ).*(-hb(i)^(1 - tau) - (-hroot)^(1 - tau)))./abs(hb(i)-hroot);
    end
end

%% Sand with root hairs
clear all
clc
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

%%
hb= linspace(-100,-10000,10); % bulk soil water potential


xx= [29.06  3      6.0e-4 % LOAM
    3.402   3.6    0.047535 %SAND
    ];
%%
Rroot=1.2261e-11;%hPa cm^-3 s;
r= 0.03+0.04; % 0.03 cm root radius in cm
r2= 1; %bulk+root+rhizo soil radius in cm
L = 1;% root segment length

hroot = 15000; %soil:root interface potential
%%
for k=2 % 1 = loam; 2 = sand
    for i = 1:10 % iteratons of bulk soil water potential
        h0=-xx(k,1);
        k0=xx(k,3);
        tau=xx(k,2);
        
        Ks_nonlinear_no(i) = (2.*pi.*r.*L.* k0/(1-tau)/(h0^(-tau))./(r./2- r.*r2^2.*(log(r2)-log(r))./(r2^2-r.^2) ).*(-hb(i)^(1 - tau) - (-hroot)^(1 - tau)))./abs(hb(i)-hroot);
    end
end
%% Loam with root hairs
clear all
clc

%%
hb= linspace(-100,-10000,10); % bulk soil water potential


xx= [29.06  3      6.0e-4 % LOAM
    3.402   3.6    0.047535 %SAND
    ];
%%
Rroot=1.2261e-11;%hPa cm^-3 s;
r= 0.03+0.04; % 0.03 cm root radius in cm
r2= 1; %bulk+root+rhizo soil radius in cm
L = 1;% root segment length

% radial rhizosphere domain
% r=[r0:0.001:r2];r=r(:);
% 
% 
% csoil=0.*r(:);
hroot = 15000; %soil:root interface potential
%%
for k=1 % 1 = loam; 2 = sand
    for i = 1:10 % iteratons of bulk soil water potential
        h0=-xx(k,1);
        k0=xx(k,3);
        tau=xx(k,2);
        
        Ks_nonlinearL_no(i) =(2.*pi.*r.*L.* k0/(1-tau)/(h0^(-tau))./(r./2- r.*r2^2.*(log(r2)-log(r))./(r2^2-r.^2) ).*(-hb(i)^(1 - tau) - (-hroot)^(1 - tau)))./abs(hb(i)-hroot);
    end
end
