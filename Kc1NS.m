function [stiffness] = Kc1NS
% Turns the 6 by 6 condensed stiffness matrix of the moment frame.

numOfStories=6;

% Material Properties
fpc=5000; 
Ec = 57*sqrt(fpc);
% COLUMNS
t_c = 24;
d_c = 32;
EI_col = 0.7*Ec*(1/12)*(t_c)*(d_c^3); % 0.7*Ig for col
% BEAMS
d_b = 30;
L_b = 25*12;
% Calculate Overhang (One Side)
l_n = L_b - d_c; % clear span from col face to col face
s_w = 25*12; % clear distance from web to web
t_overh = min([6*d_b, s_w/2, l_n/12]);
t_b = 24 + t_overh;
EI_beam = 0.35*Ec*(1/12)*(t_b)*(d_b^3); % 0.35*Ig for beams

stiffness = zeros(3*numOfStories+3,3*numOfStories+3);

for story = 1:numOfStories
    if story == 1
        H_c = 15*12;
    else
        H_c = 10.5*12;
    end

    Kelem = [24*EI_col/(H_c^3)    -6*EI_col/(H_c^2)   -6*EI_col/(H_c^2)   -24*EI_col/(H_c^3)   -6*EI_col/(H_c^2)                -6*EI_col/(H_c^2)
             -6*EI_col/(H_c^2)    4*EI_col/H_c        0                   6*EI_col/(H_c^2)     2*EI_col/H_c                     0
             -6*EI_col/(H_c^2)    0                   4*EI_col/H_c        6*EI_col/(H_c^2)     0                                2*EI_col/H_c
             -24*EI_col/(H_c^3)   6*EI_col/(H_c^2)    6*EI_col/(H_c^2)    24*EI_col/(H_c^3)    6*EI_col/(H_c^2)                 6*EI_col/(H_c^2)
             -6*EI_col/(H_c^2)    2*EI_col/H_c        0                   6*EI_col/(H_c^2)     (4*EI_col/H_c)+(4*EI_beam/L_b)   2*EI_beam/L_b
             -6*EI_col/(H_c^2)    0                   2*EI_col/H_c        6*EI_col/(H_c^2)     2*EI_beam/L_b                    (4*EI_col/H_c)+(4*EI_beam/L_b)];

    stiffness(3*story-2:3*story+3,3*story-2:3*story+3) = stiffness(3*story-2:3*story+3,3*story-2:3*story+3) + Kelem;

end

stiffness=stiffness(4:3*numOfStories+3,4:3*numOfStories+3);

%% Swap degrees of freedom
for counter=1:numOfStories

    stiffness=swapDOF(stiffness,counter,3*counter-2);

end
%% Condense stiffness matrix
Ktt=stiffness(              1:  numOfStories,              1:  numOfStories);
Ktr=stiffness(              1:  numOfStories, numOfStories+1:3*numOfStories);
Krt=stiffness( numOfStories+1:3*numOfStories,              1:  numOfStories);
Krr=stiffness( numOfStories+1:3*numOfStories, numOfStories+1:3*numOfStories);

Kt=Ktt-(Ktr/Krr)*Krt;

stiffness=Kt;


