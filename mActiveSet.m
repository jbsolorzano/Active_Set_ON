function [xmin, iter] = mActiveSet(G, c, A, b, D, d)
% Queremos minimizar el problema
% 1/2 * x'*G*x + x'*c
% s.a. A*x = b
%      D*x <= d
% In:
% G.- matriz de nxn
% c.- vector de nx1
% A.- matriz de mxn (rango(A) = m, m<=n)
% b.- vector de mx1
% D.- matriz de rxn 
% d.- vector de rx1


% Definir el working set  W c CA
sigue = 0;
iter = 0;
[m, ~]  = size(A); % restricciones de igualdad
W = [1:m];
R = length(d); % restricciones de desigualdad
% restricciones activas de desigualdad en el working set

xk = linprog(zeros(length(A),1),D,d,A,b);
%xk = [3 0 0]';
%Resolvemos subproblema cuadrático
%min (1/2)p'Gp + p'g sa Ak*p 
CA = find(D*xk==d);
while ((sigue == 0) &&(iter < 10))
    xk
    W = sort(W);
    Ak = [A; D(CA,:)]; % matriz de restricciones activas
    Ak
    r  = length(CA);
    b1 = zeros(m+r,1);
    
    g  = G*xk + c;
    if (isempty(Ak))
        % Problema sin restricciones
        p = -G\g;
    else
        p  = esp_nulo( G, Ak, g, b1);
    end
    p
    %RAMA 2
    if (norm(p)<1.e-12) %le agregamos una tolerancia
        % calculamos los multiplicadores de Lagrange despejando de las CNPO
        lambda = -(Ak*Ak')\(Ak*g);
        lambda
        % encontrar A'
        if (lambda(m+1:end) >= 0)
            xmin = xk;
            sigue = 1;
        else
            %Todas las mu [ lambdas de desigualdad (m+1:end) menores a cero 
            %y queremos quitar solo una,la mas negativa]
            j = find(lambda == min(lambda(m+1:end)), 1, 'last');
            CA(j-m) = []; % quitamos la restriccion j
            %W(end+1)=j+m;
        end

    else  % p ~= 0
        jna = find(~ismember(1:R,CA)); % restricciones no activas
        D(jna,:)*p;
        %jna;

        jna1 = jna(D(jna,:)*p > 0);
        [tj,ind] = min((d(jna1)-D(jna1,:)*xk) ./ (D(jna1,:)*p))

        if (tj < 1)
            xk = xk + tj * p;
            CA = [CA jna1(ind)]; % agregamos la restriccion que se volvió activa
        else
            xk = xk + p;
        end              

    end
    
    iter = iter +1;  
end
end
