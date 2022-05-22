% 
% Wireless Communication Systems HW2, 通訊所一年級 110064533 陳劭珩
%
% 1.Consider that a MS with a velocity v receives an unmodulated carrier
%   with a frequency f_c. The incidence angle θ(t) of the incoming wave is 
%   assumed to be uniformly distributed between –π and π.
%  (c) If f_c = 2 GHz and v is uniformly distributed between 20 km/hr and 
%      90 km/hr, find the cdf and the pdf of the observed Doppler shift 
%      via simulation.
% 
clear;
clc;
%
%
N = 1e6;              % number of samples N, say, 1 million samples
% 
% from $ doc rand, we can generate N random numbers in the interval (a, b) 
% in array r = a + (b-a).*rand(1, N), in our case incidence angle θ(t) are 
% uniformly distributed between –π and π
theta = -pi + (pi + pi).*rand(1, N);         % incidence angle θ(t)
v = (20 + (90-20).*rand(1, N))*1000/(60*60); % MS velocity v = 20~90 km/hr
% 
fc = 2*1e9;           % unmodulated carrier frequency f_c = 2 GHz
vc = 299792458;       % speed of light v_c = 3*10^8 m/sec
fm = v / (vc/fc);     % f_m = v / λ_c = v / (v_c/f_c)
% 
% because fm is vector, to perform elementwise multiplication, use ".*" 
fDn = fm.*cos(theta); % doppler frequency shift f_Dn = f_m*cos(θ(t))
% 
% [N, edges] = histcounts(X, edges, 'Normalization','pdf', ...
%                         'BinLimits', [BMIN, BMAX]);
% X - Data to distribute among bins, specified as a vector, or n-D array.
% edges - Number of bins, specified as a positive integer. If not specified
%         then histcountsautomatically calculates how many bins to use
%         based on the values in X.
% 'Normalization' - Type of normalization, 'count'(default) | 'probability'
%                   | 'countdensity' | 'pdf' | 'cumcount' | 'cdf'
% 'BinLimits' - Specified as a two-element vector, [bmin, bmax]
%
BMAX = 90*1000/(60*60) / (vc/fc); % v_max / λ_c = v_max / (v_c/f_c) = f_max
BMIN = (-1)*BMAX;                 % v_min / λ_c = v_min / (v_c/f_c) = f_min
[fDn_num, edges] = histcounts(fDn, 'normalization', 'pdf', ...
                              'BinLimits', [BMIN, BMAX]);
% 
% linspace(x1, x2, n) - Generates n points, and spacing between the points 
%                       is (x2 - x1) / (n - 1). Create a vector of numbers 
%                       with n evenly spaced points between x1 and x2.
% 
x = linspace(edges(1), edges(end), length(edges) - 1);  % x-coordinate
figure(1); % figure1 plot PDF
plot(x, fDn_num, 'Color', [0, 1, 0], 'LineWidth', 1.5, 'LineStyle', '-'); 
legend('PDF', 'Location', 'best');
xlabel('the observed Doppler shift (v = 20 km/hr, f_c = 2 GHz)');
ylabel('the Probability Density Function');
title('the simulated PDF of the observed Doppler shift');
% save the plot as 'HW2_c_pdf_110064533.png'
% print('HW2_c_pdf_110064533.png', '-dpng');
%
% 
figure(2); % figure2 plot CDF
CDF = cdfplot(fDn);
set(CDF, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', 'b');
legend('CDF', 'Location', 'best');
xlabel('the observed Doppler shift (v = 20 km/hr, f_c = 2 GHz)');
ylabel('the Cumulative Distribution Function');
title('the simulated CDF of the observed Doppler shift');
% save the plot as 'HW2_c_cdf_110064533.png'
% print('HW2_c_cdf_110064533.png', '-dpng');
% 
% 


