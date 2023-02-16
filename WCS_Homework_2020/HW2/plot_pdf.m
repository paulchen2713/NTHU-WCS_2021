function [counts,points] = plot_pdf(v,f_c,N) %compute the pdf of f_m*cos(X) where X~U(0,pi) where f_m=v*f_c/c where f is speed of light use N samples
%v:km/h f_c:GHz
X=pi*rand(1,N);
c=3*10^8;%(m/s)
f_m=(v*1000/3600)*(f_c*10^9)/c;
f_D=f_m.*cos(X);
[counts,points]=count_samples(f_D,10000,-f_m,f_m);
for i=1:length(points)-1
   points(i)=(points(i)+points(i+1))/2; %use midpoind to represent each interval
end
points(end)=[];
counts=counts/N;
figure
plot(points,counts)
if length(v)==1
    titleStr=sprintf('pdf (v=%dkm/h, f_c=%dGHz)',v,f_c);
else
    titleStr=sprintf('pdf (V~U(20,90)km/h, f_c=%dGHz)',f_c);
end
title(titleStr)

end

