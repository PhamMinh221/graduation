delta = max;
q = 4;
T = 0.1;

f = @(lambda) delta + lambda.^2 / (4*(1 - q*T)) - ((-lambda + sqrt(lambda.^2 + 4*delta*q*T)) / (2*q*T))^2;

lambda0 = 0.17;

lambda_sol = fsolve(f, lambda0);

disp(['Nghiệm lambda: ', num2str(lambda_sol)])

lhs = delta + lambda_sol^2 / (4*(1 - q*T));
rhs = ((-lambda_sol + sqrt(lambda_sol^2 + 4*delta*q*T)) / (2*q*T))^2;

disp(['Vế trái: ', num2str(lhs)])
disp(['Vế phải: ', num2str(delta*(1+(1-q*T)/4))])
