%% HW5 Problem 2 ELF Procedure
% Turns the 20 by 1 vector of seismic lateral forces calculated per ASCE
% 7-16 Equivalent Lateral Force Procedure.
% Proper factors have been applied to account for Redundancy
% and accidental torsion

% N-S Direction 
clc;
close all;
clear all;

numOfWalls=4; % 4 moment frames
numOfFloors=6;
Hs1=15; % Height of Story 1
Hs=10.5; % Height of Stories 2-6

wD=125+25; % 125psf + 25psf superimposed
Bf=125;
Lf=125;

R=8;
Cd=5.5;
omega=3;
Ie=1.0; % for RC II structures

rho=1.0;
d_mf=125; % distance between moment frames [ft]
M=((d_mf/2)+(0.05*Bf))/d_mf;
accTfactor=M/0.5;

% Calculate displacement for checking drift.
Ht=((numOfFloors-1)*Hs)+Hs1;
Wf=wD*Bf*Lf*0.001;
Wt=numOfFloors*Wf;

% Calculate response spectrum.
Ss=1.97;
S1=0.701;

Fa=1.2;
Fv=1.4;

SMS=Fa*Ss;
SM1=Fv*S1;

SDS=(2/3)*SMS;
SD1=(2/3)*SM1;

TS=1.0*SD1/SDS;
TL=8.0;

T=logspace(-3,2,1001);
pSa=zeros(1,length(T));

for c=1:length(T)
    
    Tn=T(1,c);
    if Tn<TS
        pSan=SDS;
    elseif ((TS<=Tn)&&(Tn<TL))
        pSan=SD1/Tn;
    else
        pSan=SD1*TL/(Tn^2);
    end
    pSa(1,c)=pSan; 
    
end
% semilogx(T,pSa);
% title('Response spectrum');
% xlabel('Period [s]');
% ylabel('Pseudo spectral acceleration [g]');
% xlim([0.01 20]);
% pbaspect([2 1 1]);
% grid on;
% hold on;

% Calculate first mode period form eigen analysis.
[~,eigValues]=eig(Kc1NS,Mc1NS);
numOfFloors=6;

[eigValues6x1,~]=sort(diag(eigValues),'ascend');
w1=sqrt(eigValues6x1(1,1));
T1=(2*pi)/w1;

Ta=0.02*(Ht^0.75);
Cu=1.4;
Tu=Cu*Ta;

if T1<Tu
    TforELFP=T1;
else
    TforELFP=Tu;
end

Cs=(interp1(T,pSa,TforELFP))/(R/Ie);

if Cs<0.044*SDS*Ie
    Cs=0.044*SDS*Ie;
end

if Cs<0.5*S1/(R/Ie)
    Cs=0.5*S1/(R/Ie);
end

% Calculate Base Shear
V_b=Cs*Wt/numOfWalls

% Calculate k factor
if TforELFP<0.5
    k=1;
elseif TforELFP>2.5
    k=2;
else
    k=interp1([0.5 2.5],[1 2],TforELFP);
end

% Distribute base shear over height

sumWiHik=0;
for c=1:numOfFloors
    if c==1
    sumWiHik=sumWiHik+(Wf*(c*Hs1)^k);
    else
    sumWiHik=sumWiHik+(Wf*(c*Hs)^k);
    end 
end

FL=zeros(3*numOfFloors,1);
Cv=zeros(numOfFloors,1);

for c=1:numOfFloors
    if c==1
    Cv(c,1)=(Wf*(c*Hs1)^k)/sumWiHik;
    FL(3*c-2,1)=Cv(c,1)*V_b;
    else
    Cv(c,1)=(Wf*(c*Hs)^k)/sumWiHik;
    FL(3*c-2,1)=Cv(c,1)*V_b;
    end
end

% Display Floor Loads
FL

% Magnify for redundancy.
FL_AT=accTfactor*rho*FL
% Determine Story Displacements
FL1 = [FL(1) FL(4) FL(7) FL(10) FL(13) FL(16)];
K_NS = Kc1NS();
delta_el = FL1/(K_NS);
for c=1:numOfFloors
     delta(c) = Cd*delta_el(c)/Ie;
end

% Determine Story Drifts
drifts = zeros(6,1);

drifts(1) = delta(1);
drifts(2) = delta(2) - delta(1);
drifts(3) = delta(3) - delta(2);
drifts(4) = delta(4) - delta(3);
drifts(5) = delta(5) - delta(4);
drifts(6) = delta(6) - delta(5);
drifts
delta = delta'

% Determine Story Shear
V_s=zeros(numOfFloors,1);
V_s(6)=FL(16,1);
V_s(5)=FL(16,1)+FL(13,1);
V_s(4)=FL(16,1)+FL(13,1)+FL(10,1);
V_s(3)=FL(16,1)+FL(13,1)+FL(10,1)+FL(7,1);
V_s(2)=FL(16,1)+FL(13,1)+FL(10,1)+FL(7,1)+FL(4,1);
V_s(1)=FL_AT(1,1)+FL(4,1)+FL(7,1)+FL(10,1)+FL(13,1)+FL(16,1);
V_sn=V_s./Wt;

V_s
V_sn

% Determine Overturning Moments
H=[Hs1;Hs1+Hs;Hs1+Hs*2;Hs1+Hs*3;Hs1+Hs*4;Hs1+Hs*5];

OTM = zeros(6,1);
OTM(6)=FL(16,1)*H(6);
OTM(5)=FL(16,1)*H(6)+FL(13,1)*H(5);
OTM(4)=FL(16,1)*H(6)+FL(13,1)*H(5)+FL(10,1)*H(4);
OTM(3)=FL(16,1)*H(6)+FL(13,1)*H(5)+FL(10,1)*H(4)+FL(7,1)*H(3);
OTM(2)=FL(16,1)*H(6)+FL(13,1)*H(5)+FL(10,1)*H(4)+FL(7,1)*H(3)+FL(4,1)*H(2);
OTM(1)=FL_AT(1,1)*H(1)+FL(4,1)*H(2)+FL(7,1)*H(3)+FL(10,1)*H(4)+FL(13,1)*H(5)+FL(16,1)*H(6);

OTM

figure(1); % Initialize window of plots
subplot(1,5,1); % Subplot: one row, 5 columns of plots, 1st plot (1,5,1)
hold on; %IDK
grid on; %grid lines in plot
title ('Story Displacements') %title of subplot
ylabel('Floors'), xlabel('Displacement (in)')
yticks([0 1 2 3 4 5 6])
plot([0 delta'],0:6) %x axis = delta as a row with a 0 in front
hold off;

subplot(1,5,2);
hold on;
title ('Story Drift')
ylabel('Floors'), xlabel('Drift (in)')
grid on;
yticks([0 1 2 3 4 5 6])
plot([0 drifts'],0:6)
hold off;

subplot(1,5,3);
hold on;
title ('Story Shears')
ylabel('Floors'), xlabel('Shear (lbs)')
grid on;
yticks([0 1 2 3 4 5 6])
plot([V_s(1) V_s'],0:6)
hold off;

subplot(1,5,4);
hold on;
title ('Story Shears (Normalized)'), xlabel('Shear')
ylabel('Floors')
yticks([0 1 2 3 4 5 6])
grid on;
plot([V_sn(1) V_sn'],0:6)
hold off;

subplot(1,5,5);
hold on;
title ('Story OTM')
ylabel('Floors'), xlabel('Moment (lb-in)')
grid on;
yticks([0 1 2 3 4 5 6])
plot([OTM(1) OTM'],0:6)
hold off;

sgtitle('North-South Structural Responses')  
