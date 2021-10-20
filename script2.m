% Problema de Klee-Minty
% n = 15
% Lo llevamos a la forma cuadrática
    % 1/2 * x'*G*x + x'*c
    % s.a. A*x = b
    %      D*x <= d
n = 15;
G = (10^(-4))*eye(15);
c = -ones(15,1);
A = []; % no hay restricciones de igualdad
b = []; % no hay restricciones de igualdad
D = [2*tril(ones(n),-1)+eye(n);-eye(n)];
d = [(2*ones(1,15)).^[1:n]-1 zeros(1,n)]';
x0 = zeros(15,1);
W0 = [21:30];
maxIter = []; %dejamos el máximo de iteraciones por default

%
[x_opt, iter,valor_min]= mActiveSet(G, c, A, b, D, d,x0,[],maxIter);

% Imprimimos la solución:
fprintf('Vector óptimo obtnenido es : [');
fprintf('%g, ', x_opt(1:end-1));
fprintf('%g]\n', x_opt(end));

% Imprimimos el valor mínimo:
fprintf('\nCon valor mínimo: %.4f\n', valor_min); 

% Valor óptimo del método quadprog de MatLab
[x_matlab,fval] = quadprog(G,c,D,d);
fprintf('Vector óptimo de MatLab: [');
fprintf('%g, ', x_matlab(1:end-1));
fprintf('%g]\n', x_matlab(end));

% Imprimimos el valor mínimo de MatLab:
fprintf('\nCon valor mínimo: %.4f\n', fval); 

msgbox({sprintf('Solución óptima obtnenida:  %.4f.',valor_min)
        sprintf('Número de iteraciones: %d.',iter)
        sprintf('\nSolución de matlab es: %.4f .',fval)},'Resultados');
