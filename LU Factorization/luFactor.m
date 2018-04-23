function [L, U, P] = luFactor(A)
%This function finds the upper triangular matrix (U), the lower triangular
%matrix (L), and the pivot matrix (P) of the coefficient matrix (A) entered.
%   A is a matrix of coefficients
%   U is a upper triagular matrix made of the coefficients obtained through
%       forward elimination
%   L is a lower triangular matrix made of the coefficients obtained
%       through backward substitution
%   P is a matrix that tracks the pivoting of the rows of A throughout the
%       Gauss elimination method
%% Presetting P, L and U matrices
[r,c]=size(A); %computes number of rows and columns
L=zeros(r); %creates L as a zero matrix the same size as A
P=zeros(r); %creates P as a zero matrix the same size as A
n=1;
while n<=r %adds ones across the diagonal of P 
    P(n,n)=1;
    n=n+1;
end
U=A; %creates U as a matrix the same as A
%% Eliminating invalid inputs

if nargin>1 || nargin<1 %only one input needed
    error('Cannot input more than one matrix');
elseif r~=c %matrix must be square
    error('Matrix must be square');
end
%% Checking for pivoting in first column
z=abs(A(1,1)); %absolute value of first variable in first column
y=abs(A(:,1));%absolute value of all other first column values
y=max(y);%max value of y
if y>z %deciding if pivoting is neccesary
    i=2;%row
    j=1;%column
    xold=z;%xold is the absolute value of A(1,1)
    while i<=r %checking which row has the largest value in the first column
        x=abs(A(i,j));
        if x>xold%compares next row value to old row value
            pivotrow=i; %keeps row number if x is greater
            xold=x;%keeps old x so that it can be compared to the value in the next row 
        end
        i=i+1;%changes row in next iteration
    end
P([1 pivotrow],:)=P([pivotrow 1],:);%updating P matrix
U([1 pivotrow],:)=U([pivotrow 1],:);%updating U matrix
L([1 pivotrow],:)=L([pivotrow 1],:);%updating L matrix
end
%% Eliminating first variable
h=2;%row
g=1;%column
t=2;%kept to reset h after inner while loop runs
while g<c %while column number is less than the number of columns in A
    while h<=r
        eliminator=U(h,g)/U(g,g);%determining multiplyer to cancel out the variable
        U(h,:)=U(h,:)-eliminator*U(g,:);%updating U matrix
        L(h,g)=eliminator;%puts eliminator value in the L matrix
        h=h+1;
    end
    t=t+1;%used to reset h
    h=t;%resets h so that its only one larger than the last outside while loop
    g=g+1;%adds one to the number of columns
            z=abs(U(g,g));%absolute value of variable to be eliminated
            y=abs(U(:,g));%absolute value of all other values in the rest of that row
            y=max(y);%max of y
            if y>z %deciding if pivoting is neccesary
                i=h;%rows
                j=g;%columns
                xold=z;%sets xold as the first value in the column of the row that is being eliminated
                while i<=r %checking which row has the largest value in the column
                    x=abs(U(i,j));%absolute value of all the values in that column
                    if x>xold %checks if the value in that row is bigger than the value in the previous row
                        pivotrow=i; %keeps row number with largest value 
                        xold=x;
                    end
                    i=i+1;
                   
                end
                P([g pivotrow],:)=P([pivotrow g],:);%updating P matrix
                U([g pivotrow],:)=U([pivotrow g],:);%updating U matrix
                L([g pivotrow],:)=L([pivotrow g],:);%updating L matrix
            end
    
end
n=1;
while n<=r %adds ones across the diagonal of L
    L(n,n)=1;
    n=n+1;
end
end

