function I = Simpson(x,y)
%This function approximates the integral of the equation that cooresponds to
%the values of x and y using the composite Simpson's one-third rule.
%   The approximation is computed using the composite Simpson's one-third rule and a
%   a trapeziod if the inputted x and y have an even number of elements.
%   x=vector containing evenly spaced values that represent x values of a
%       function over a certain interval
%   y=vector containing y values of a function over the same interval
%   I=approximate integral value
%% Checking for invalid inputs
if nargin~=2 %making sure theres only 2 inputs
    error('There can only be 2 inputs')
end
if range(x(2:end)-x(1:end-1))~=0 %checking to make sure x is evenly spaced
    error('x vector must be evenly spaced')
end
sizex=size(x);
sizey=size(y);
if sizex(1,1)~=1 && sizey(1,1)~=1 %making sure x and y only have one row
    error('x and y must be row vectors')
elseif sizex(1,2)~=sizey(1,2)%making sure x and y are the same length
    error('x and y must have the same number of elements')
elseif rem(sizex(1,2),2)==0
    warning('Trapeziodal rule required for last interval')
    trap=1;%keeping track of whether or not to use the trapezoidal rule
else
    trap=0;
end

%% Computing composite simpsons rule
if sizex(1,2)==3 %usessimpsons one-third rule if there are only 3 elements
    I=(x(1,3)-x(1,1))*((y(1,1)+4*y(1,2)+y(1,3))/6);
elseif sizex(1,2)>3 && trap==0 %uses composite simpsons one-third rule
    i=1;
    sumodds=0;
    while i<=sizey(1,2)%this loop adds up all the odd number elements of y
        sumodds=y(1,i)+sumodds;
        i=i+2;
    end
    n=2;
    sumevens=0;
    while n<=sizey(1,2) %this loop adds up all the even number elements
        sumevens=y(1,n)+sumevens;
        n=n+2;
    end
    steps=sizex(1,2)-1;%calulatesnumber of steps in x
    I=(x(1,end)-x(1,1))*((y(1,1)+4*(sumodds)+2*(sumevens)+y(1,end))/(3*steps));
else
    i=1;
    sumodds=0;
    while i<=sizey(1,2)%this loop adds up all the odd number elements of y
        sumodds=y(1,i)+sumodds;
        i=i+2;
    end
    n=2;
    sumevens=0;
    while n<sizey(1,2) %this loop adds up all the even number elements
        sumevens=y(1,n)+sumevens;
        n=n+2;
    end
    steps=sizex(1,2)-1;%calulatesnumber of steps in x
    I=(x(1,end)-x(1,1))*((y(1,1)+4*(sumodds)+2*(sumevens)+y(1,end))/(3*steps))+(x(1,end)-x(1,end-1))*((y(1,end)+y(1,end-1))/2);%calculates integral using composite rule and a trapezoid
end
end

