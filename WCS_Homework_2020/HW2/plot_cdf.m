function  []=plot_cdf(counts,points,v,f_c)
total=sum(counts);
counts=counts/total;
for i=1:length(counts)
    cdf(i)=sum(counts(1:i));
end
figure
plot(points,cdf)
if length(v)==1
    titleStr=sprintf('cdf (v=%dkm/h, f_c=%dGHz)',v,f_c);
else
    titleStr=sprintf('cdf (V~U(20,90)km/h, f_c=%dGHz)',f_c);
end
title(titleStr)
end

