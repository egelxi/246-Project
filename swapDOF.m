function[C]=swapDOF(A,i,j)

[numOfDOF,~]=size(A);

%% Swap rows
swapMatrix=zeros(numOfDOF,numOfDOF);
for c=1:numOfDOF
    
    if c==i
        swapMatrix(c,j)=1;
    elseif c==j
        swapMatrix(c,i)=1;
    else
        swapMatrix(c,c)=1;
    end
    
end
B=swapMatrix*A;

%% Swap columns
C=B*swapMatrix;