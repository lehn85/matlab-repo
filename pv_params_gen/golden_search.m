function [ min,x ] = golden_search( func, range, errTolerance)
%GOLDEN_SEARCH Summary of this function goes here
%   Detailed explanation goes here
a=range(1);
b=range(2);
tau=double((sqrt(5)-1)/2);
x1=a+(1-tau)*(b-a);
x2=a+tau*(b-a);
f_x1=0;
f_x2=0;
maxIter=100;
iter=0;
while(iter<maxIter && abs(b-a)>errTolerance)
    iter=iter+1;
    f_x1=func(x1);
    f_x2=func(x2);
    if (f_x1<f_x2)
        b=x2;
        x2=x1;
        x1=a+(1-tau)*(b-a);
    else
        a=x1;
        x1=x2;
        x2=a+tau*(b-a);
    end
end
x=(x1+x2)/2;
min=(f_x1+f_x2)/2;
end

