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
W = [];
sigue = 0;
iter = 0;
[m, ~]  = size(A); % restricciones de igualdad
R = length(d); % restricciones de desigualdad
% restricciones activas de desigualdad en el working set

xk  = punto_x0(A,D,b,d);

%Resolvemos subproblema cuadrático
%min (1/2)p'Gp + p'g sa Ak*p = phi

while (sigue == 0)
    CA = find(D*xk==d);
    Ak = [A; D(W,:)]; % matriz de restricciones activas
    
    g  = Q*xk + c;
    b1 = zeros(m+r,1);
    
    if (isempty(Ak))
        % Problema sin restricciones
        p = -G\g;
    else
        p  = esp_nulo( G, Ak, g, b1);
    end
    
    %RAMA 2
     if (norm(p)<1.e-12) %le agregamos una tolerancia
        % calculamos los multiplicadores de Lagrange
        mu;
        
        if (mu >= 0)
            xmin = xk;
            sigue = 1; 
        end
     end
iter = iter +1;  

end
end
