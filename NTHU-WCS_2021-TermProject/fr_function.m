% 
% f_r() function defined in Eq.(7.35)
%
% ref. Fundamentals of Wireless Communication 2005.
%
function  y = fr_function(omega, N, delta, L)
    %
    % f_r(£[_r) = (1/n_r) * exp(j£k*£G_r*£[_r(n_r - 1)) * (sin(£k*L_r*£[_r)) 
    %           * (sin(£k*L_r*£[_r / n_r)), where L_r = n_r * £G_r
    %
    y = (1/N) * exp(1i*pi*delta*omega*(N-1)) .* sin(pi*L*omega) ...
        ./ sin(pi*L*omega/N) ;
end
