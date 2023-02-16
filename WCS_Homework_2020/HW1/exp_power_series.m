function ans = exp_power_series(x,m) % The first m terms of the power series of e and take the value at x
temp(1)=1;
for i=2:m+1
    temp(i)=temp(i-1)*x/(i-1);
end
ans=sum(temp);
end

