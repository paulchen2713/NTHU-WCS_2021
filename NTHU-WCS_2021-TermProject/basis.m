%
% Compute unit spatial signatures e_r(£[), e_t(£[) basis vector
%
% ref. p.351~352, Fundamentals of Wireless Communication 2005.
%
function basis_vec_e = basis(omega, N, delta)
    %
    % unit spatial signature e_r(£[), e_t(£[), 
    % defined in Eq.(7.21), Eq.(7.25)
    %
    basis_vec_e = zeros(N, 1);
    for i = 1 : N
        basis_vec_e(i) = exp(-1*1i*2*pi*(i-1)*delta*omega);
    end
    basis_vec_e = basis_vec_e ./ sqrt(N);
end
