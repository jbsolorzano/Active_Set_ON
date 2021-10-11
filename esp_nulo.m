function [x] = esp_nulo( G, A, c, b)
% Método del espacio nulo para
% Min (1/2)*x'*G*x + x'*c
% S.A. A*x = b
% donde G es una matrix simétrica positiva definida, A es mxn con m<= n
% rango(A) = m.

Z = null(A); % base ortonormal del espacio nulo de A.
xp = A\b;    % solución particular de A*x = b
%En las notas: resolver (Z'*G*Z)* xh=- Z'*(G*xp + c)
G1 = Z'*G*Z; % matriz simétrica positiva definida en R^(n-m)
c1 = Z'*(G*xp + c); % lado derecho
% resolver el sistema lineal: G1*xh = -c1 

% Cholesky
L = chol(G1)';
y = -L\c1;
xh = L'\y;

%solución final
x = xp + Z*xh;
end