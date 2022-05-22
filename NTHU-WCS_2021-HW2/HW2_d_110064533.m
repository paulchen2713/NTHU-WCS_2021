% 
% Wireless Communication Systems HW2, 通訊所一年級 110064533 陳劭珩
%
% 1.Consider that a MS with a velocity v receives an unmodulated carrier
%   with a frequency f_c. The incidence angle θ(t) of the incoming wave is 
%   assumed to be uniformly distributed between –π and π.
%  (d) Derive the cdf and the pdf of the observed Doppler shift for fixed v
%      and fc. Compare the simulation results with the theoretical results.
% 
clear;
clc;
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% simulated result, settings are the same as before.
%
N = 1e6; % number of samples N, say, 1 million samples
% 
theta_simulated = -pi + (pi + pi).*rand(1, N); % incidence angle θ(t)
fDn_simulated = cos(theta_simulated); % normalized fDn = fDn/fm = cosθ(t)
%
[fDn_num, edges] = histcounts(fDn_simulated, 'normalization', 'pdf', ...
                              'BinLimits', [-1, 1]);
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% theoretical result, settings are base on the derived formula
% 
theta_theoretical = linspace((-1)*pi, pi, 200);
fDn_theoretical = cos(theta_theoretical);
pdf_theoretical = 1 ./ (pi * sqrt(1 - (fDn_theoretical).^2));
cdf_theoretical = acos((-1)*fDn_theoretical) / pi;
% cdf_theoretical = 1 - (acos(fDn_theoretical) / pi);
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%
% plot the comparison results
%
figure(1); % figure1 plot the comparison of PDFs
%
% first, the simulated PDF, in red solid line
x = linspace(edges(1), edges(end), length(edges) - 1);  % x-coordinate
plot(x, fDn_num, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '-'); 
hold on;
%
% then, the theoretical PDF, in green, dash-dot line
plot(cos(theta_theoretical), pdf_theoretical, 'g-.', 'LineWidth', 1.2);
%
xlabel('the observed Doppler shift (normalized)');
ylabel('the Probability Density Function');
title('the comparison between simulated PDF and theoretical PDF');
legend('simulated PDF', 'theoretical PDF', 'Location', 'best');
hold off;
% save the plot as 'HW2_d_pdf_110064533.png'
% print('HW2_d_pdf_110064533.png', '-dpng');
%
%
figure(2); % figure2 plot the comparison of CDFs
% 
% first, the simulated CDF, in blue solid line
simulated_CDF = cdfplot(fDn_simulated);
set(simulated_CDF, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', 'b');
hold on;
% 
% then, the theoretical CDF, in cyan, dash-dot line
plot(cos(theta_theoretical), cdf_theoretical, 'c-.', 'LineWidth', 1.2);
%
xlabel('the observed Doppler shift (normalized)');
ylabel('the Cumulative Distribution Function');
title('the comparison between simulated CDF and theoretical CDF');
legend('simulated CDF', 'theoretical CDF', 'Location', 'best');
hold off;
% save the plot as 'HW2_d_cdf_110064533.png'
% print('HW2_d_cdf_110064533.png', '-dpng');
%
%

