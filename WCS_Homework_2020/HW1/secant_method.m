function [ans,steps,value] = secant_method(f,x0,x1,M,e) %x0,x1 are initial points, M : max steps, e : error
if nargin<4
    M=20;
    e=10^-3;
elseif nargin<5
    e=10^-3;
end
for i=1:M
    x2=x1-f(x1)*(x1-x0)/(f(x1)-f(x0));
    x0=x1;
    x1=x2;
    if(abs(f(x2))<e)
        ans=x2;
        steps=i;
        value=f(x2);      
        return
    end
end
error('Does not find roots in given steps')
end

