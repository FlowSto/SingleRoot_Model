%% Extraction and estimation of radial conductance from Meunier et al. 2018, Maize neutron radiography root hydraulics
%%%%%%%%%%%%%%%%%%%%%
clc; clear all;
L = [40 40 40 30;                                                           % axially measured root length
    30 30 30 15];                                                           % radially measured root length
l = 0.1;                                                                    % segment length
N = round(L./l);                                                            % number of segment
% piecewise function (Meunier) for kx Brace
x1 = linspace(0.78*10^-4/(86400), 2.56*10^-4/(86400), 169);   % 0 - 16.9 cm
x2 = linspace(2.56*10^-4/(86400), 121.8*10^-4/(86400), (27-16.9)/l);  % 16.9 - 27
x3 = linspace(121.8*10^-4/(86400), 195.6*10^-4/(86400), (32.6-27)/l);  % 27 - 32.6
x4 = linspace(195.6*10^-4/(86400), 225.3*10^-4/(86400), (40+l-32.6)/l);  % 32.6 - 40
kx.brace = fliplr([x1 x2 x3 x4]);                                   %cm^4 hpa^-1 s^-1
% piecewise function (Meunier) for kr Brace
x11 = linspace(0.48*10^-4/(86400), 0.48*10^-4/(86400), 19/l);  % 0 - 19 cm
x12 = linspace(0.48*10^-4/(86400), 0.028*10^-4/(86400), 1/l);  % 19 - 20
x13 = linspace(0.028*10^-4/(86400), 0.02*10^-4/(86400), (10+l)/l); % 20 - 30 x3 = linspace(0.028*10^-4/(86400), 0.02*10^-4/(86400), 10/l)
kr.brace = fliplr([x11 x12 x13]);                                   %cm^4 hpa^-1 s^-1
r.brace = 0.067;
% Crown Root
x21 = linspace(0.82*10^-4/(86400), 3.71*10^-4/(86400), 11.2/l);  % 0 - 11.2 cm
x22 = linspace(3.71*10^-4/(86400), 24.44*10^-4/(86400), (28.7-11.2)/l); % 11.2 - 28.7
x23 = linspace(24.44*10^-4/(86400), 73.63*10^-4/(86400), (34.3-28.7)/l); % 28.7 - 34.3
x24 = linspace(73.63*10^-4/(86400), 73.63*10^-4/(86400), (40+l-34.3)/l); % 34.3 - 40
kx.crown = fliplr([x21 x22 x23 x24]); %cm^4 hpa^-1 s^-1
% piecewise function (Meunier) for kr Crown
x31 = linspace(0.39*10^-4/(86400), 0.25*10^-4/(86400), 19/l);  % 0 - 19 cm
x32 = linspace(0.25*10^-4/(86400), 0.17*10^-4/(86400), 1/l);   % 19 - 20
x33 = linspace(0.17*10^-4/(86400), 0.17*10^-4/(86400), (10+l)/l);  % 20 - 30
kr.crown = fliplr([x31 x32 x33]);                                   %cm^4 hpa^-1 s^-1
r.crown = 0.052;
% piecewise function (Meunier) for kx Seminal
x61 = linspace(1.76*10^-4/(86400), 1.91*10^-4/(86400), 7/l);   % 0 - 7 cm
x62 = linspace(1.91*10^-4/(86400), 5.45*10^-4/(86400), 3/l);   % 7 - 10
x63 = linspace(5.45*10^-4/(86400), 14.85*10^-4/(86400), 11/l); % 10 - 21
x64 = linspace(14.85*10^-4/(86400), 15.85*10^-4/(86400), (19+l)/l);% 21 - 40
kx.seminal = fliplr([x61 x62 x63 x64]);                               %cm^4 hpa^-1 s^-1
% piecewise function (Meunier) for kr Seminal
x71 = linspace(0.47*10^-4/(86400), 0.47*10^-4/(86400), 19/l);  % 0 - 19 cm
x72 = linspace(0.47*10^-4/(86400), 0.43*10^-4/(86400), 1/l);   % 19 - 20
x73 = linspace(0.43*10^-4/(86400), 0.28*10^-4/(86400), (10+l)/l);  % 20 - 30
kr.seminal = fliplr([x71 x72 x73]);                                   %cm hpa^-1 s^-1
r.seminal = 0.037;
% piecewise function (Meunier) for kx Lateral
x41 = linspace(0.74*10^-4/(86400), 0.98*10^-4/(86400), 4.8/l);   % 0 - 4.8 cm
x42 = linspace(0.98*10^-4/(86400), 5.04*10^-4/(86400), (15.2-4.8)/l);  % 4.8 - 15.2
x43 = linspace(5.04*10^-4/(86400), 5.99*10^-4/(86400), (30-15.2)/l);  % 15.2 - 30
x44 = linspace(5.99*10^-4/(86400), 5.99*10^-4/(86400), (10)/l);
kx.lateral = fliplr([x41 x42 x43 x44]);                               %cm^4 hpa^-1 s^-1
% piecewise function (Meunier) for kr Lateral
kr.lateral = fliplr(linspace(5.62*10^-5/(86400), 5.62*10^-5/(86400), (15)/l));                                   %cm^4 hpa^-1 s^-1
r.lateral = 0.03;
save('Meunier_Conductances.mat','kx','kr', 'r')