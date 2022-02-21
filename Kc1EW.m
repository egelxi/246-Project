 function [stiffness] = Kc1EW
% Turns the 6 by 6 condensed stiffness matrix of the shear wall.

numOfStories=6;
tw=18;
Lw=27*12;
fpc=5000;
Ec=57*sqrt(fpc);

stiffness=zeros(2*numOfStories+2,2*numOfStories+2);

for story=1:numOfStories

    elemNodesNumber=[story,story+1];
    
    if story==1
        L=15*12;
    else
        L=10.5*12;
    end
    
    EIz=0.8*Ec*(1/12)*(tw)*(Lw^3); % 0.8*Ig

    
    Kelem=[  12*EIz/(L^3),  -6*EIz/(L^2), -12*EIz/(L^3),  -6*EIz/(L^2);
   
             -6*EIz/(L^2),     4*EIz/(L),   6*EIz/(L^2),     2*EIz/(L);
              
            -12*EIz/(L^3),   6*EIz/(L^2),  12*EIz/(L^3),   6*EIz/(L^2);
                
             -6*EIz/(L^2),     2*EIz/(L),   6*EIz/(L^2),     4*EIz/(L)];
    
    
    for c1=1:2
        for c2=1:2
            
            nodeI=elemNodesNumber(c1);
            nodeJ=elemNodesNumber(c2);
            stiffness(2*nodeI-1:2*nodeI,2*nodeJ-1:2*nodeJ)=stiffness(2*nodeI-1:2*nodeI,2*nodeJ-1:2*nodeJ)+...
                                                               Kelem(   2*c1-1:   2*c1,   2*c2-1:   2*c2);
        end  
    end
    
end

stiffness=stiffness(3:2*numOfStories+2,3:2*numOfStories+2);


%% Swap degrees of freedom
for counter=1:numOfStories
    
    stiffness=swapDOF(stiffness,counter,2*counter-1);
    
end


%% Condense stiffness matrix
Ktt=stiffness(              1:  numOfStories,              1:  numOfStories);
Ktr=stiffness(              1:  numOfStories, numOfStories+1:2*numOfStories);
Krt=stiffness( numOfStories+1:2*numOfStories,          1:numOfStories);
Krr=stiffness( numOfStories+1:2*numOfStories, numOfStories+1:2*numOfStories);

Kt=Ktt-(Ktr/Krr)*Krt;

stiffness=Kt;    

