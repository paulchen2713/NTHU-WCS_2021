function ans = ln_factorial_approx(n) % approximation of n! in ln scale
ans=n.*log(n)-n+log(2*pi*n)/2;
end


