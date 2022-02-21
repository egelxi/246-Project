function [mass] = Mc1NS

numOfStories=6;
numOfBays=4;

Bf=125;
Lf=125;
wD = 125+25; % psf
Wf = (wD*Bf*Lf); % lbs
% wL = 55; % psf
% Wf=0.001*wD*Bf*Lf;
% Wf = 1.0*(wD) + 0.25*(wL); % weight per floor
g=32.2; %ft/s^2

mass=eye(numOfStories)*(Wf/g/numOfBays); %lbs/(ft/s^2)