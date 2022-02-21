function [mass] = Mc1EW

numOfStories=6;
numOfWalls=2;

Bf=125;
Lf=125;
wD = 125+25; % psf
Wf = (wD*Bf*Lf); % lb
% wL = 55; % psf
% Wf=0.001*wD*Bf*Lf;
% Wf = 1.0*(wD) + 0.25*(wL); % weight per floor
g=32.2;

mass=eye(numOfStories)*(Wf/g/numOfWalls); %lb