function ans = TOT(m,P) % TOT : total offered traffic, P : Blocking Prob., m : number of channels
if m<=5
    k=[0:m];
    ans=secant_method(@(x) P*sum(x.^k./factorial(k))-(x^m)/factorial(m),m,m+1); % solve x , 1,2 are initial point
else %Taking ln scale since m! will become Inf as m large 
    %Fuethermore, approximation of ln(m!) is more precisely as m large
    ans=secant_method(@(x) log(P)+ln_factorial_approx(m)+log(exp_power_series(x,m))-m*log(x),1,2); % solve x , 1,2 are initial point
end

