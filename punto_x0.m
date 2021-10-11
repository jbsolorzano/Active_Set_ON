function xk = punto_x0(A, D, b, d)
% Encontrar un punto factible.
% Dado el problema: Min (1/2)x'Gx + c'x s.a. Ax = b Dx <= d
% Llevamos las restricciones a la forma estándar de un problema de
% programación lineal para encontrar un punto factible

% Llevando el problema a forma estándar
[m,n] = size(A);
[p,n] = size(D);
% agregando variables de holgura y separando a x en x+,x-
A1 = [A -A zeros(m, p); D -D eye(p)];
b1 = [b;d];
% agregamos una z auxiliar
A1 = [ A1 eye(m+p)];
c = [zeros(2*n + p,1); ones(m+p,1)];
% Resolvemos el problema min e'z sa A1x + z = b1
X = linprog(c,[],[],A1, b1, zeros(2*n + p + m + p),[]);
xk = X(1:n) - X((n+1):(2*n));

end