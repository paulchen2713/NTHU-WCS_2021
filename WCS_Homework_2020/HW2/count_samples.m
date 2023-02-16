function [counts,points] = count_samples(x,N,a,b)% x : input array ,divide [a,b] into N interval, count how many terms lie in each subinterval
if nargin==1
    N=10^4;
    a=min(x);
    b=max(x);
    
elseif nargin==2
    a=min(x);
    b=max(x);
elseif nargin==3
    b=max(x);
end
points=[a:(b-a)/N:b];
counts=zeros(1,N);
for i=1:length(x)
    for j=1:N
        if(x(i)>=points(j)&&x(i)<=points(j+1))
            counts(j)=counts(j)+1;
            break
        end
    end
end
end

