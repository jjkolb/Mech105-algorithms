function [xr,fx,ea,iter] = falseposition(func,xl,xu,es,maxiter)
%This function will find the root of any equation using the false position
%method. This function will only work if the inputted upper and lower
%guesses encapsulate a root. The function returns the estimated root (root),
%the function evaluated at that root (fx), the approximate relative error
%of the root (ea), and the number of iterations performed (iter).

%func is whatever function the user wants evaluated
%xl is the lower guess of the root
%xu is the upper guess of the root
%es is an optional argument where the user can input the desired relative
%   error. The default relative error is 0.0001%
%maxiter is another optional argument where the user can input the desired
%   number of iterations. The default number of iterations is 200


%checking for valid inputs
if nargin==3 %defaults values for error and iterations
    es=0.0001; 
    maxiter=200;
elseif nargin<3 %there has to be at least 3 inputs
    error('Invalid number of inputs');
elseif nargin==4
    maxiter=200;
elseif nargin>5 %there cannot be more than 5
    error('Invalid number of inputs');
elseif es<0
    error('Approximate error cannot be negative');
elseif maxiter<0
    error('Desired number of iterations cannot be negative');
end

%checking the first upper and lower guesses have different signs
prodofguesses=func(xl)*func(xu);
if prodofguesses>0
    error('xl and xu must be estimates that bracket the root');
end

%perform false position method
iter=0; %tracks iterations
ea=100; %tracks approximate error
while ea>es && iter<maxiter
    xr=double(xu-(func(xu)*(xl-xu))/(func(xl)-func(xu))); %false position formula
    
    %calculating new approximate error
    if iter==0
        xrold=0;
    end
    ea=abs((xr-xrold)/xr)*100;

    
    %determining new interval
    prodrootandxl=func(xr)*func(xl);
    if prodrootandxl<0
        xu=xr;
    elseif prodrootandxl>0 && xr>0
        xu=xr;
    else
        xl=xr;
    end
    iter=iter+1; %adds to the number of iterations
    xrold=xr; %keeps old root estimate to determine approximate error in next iteration
end
root=xr;
fx=double(func(root));
end

