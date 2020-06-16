function residual=isotherm_regression(z,q_NL,c_NL)
% This function sets up the nonlinear Langmuir isotherm model. Note maximum
% membrane loading (Q) is denoted by z(2) and membrane equilibrium constant
% (K)is denoted by z(1).

NL_model = (z(1)*z(2)*c_NL) ./ (1+z(1)*c_NL);

r = q_NL-NL_model;
residual = r.'*r;