function [xmin, iter,valor_min] = mActiveSet(G, c, A, b, D, d, x0, W0, maxIter)
% Queremos minimizar el problema
% 1/2 * x'*G*x + x'*c
% s.a. A*x = b :indices
%      D*x <= d indices 
% In:
% G.- matriz de nxn
% c.- vector de nx1
% A.- matriz de mxn (rango(A) = m, m<=n)
% b.- vector de mx1
% D.- matriz de rxn 
% d.- vector de rx1
% x0.- vector de nx1, solución particular
% W0.- Conjunto de índices iniciales donde se da igualdad (respecto a x0);

%NOTA_1: Índices ={1,2,...,m,m+1,m+2,...,m+r}. Los primeros m corresponden a
    % las restricciones de igualdad; los últimos r índices corresponden a
    % las desigualdades.
%NOTA_2: Si el punto x0 se manda vacío, se encuentra una solución
%particular con ligprog.
%NOTA_3: Si el conjunto W0 se manda vacío, se toma por default W0= A(x0).



if (isempty(maxIter))
    maxIter = 100;
end
bandera = 0;
iter = 1;
[m, ~]  = size(A); % restricciones de igualdad
if (isempty(A))
    E = [];
    m = 0;
else
    [m, ~]  = size(A); % restricciones de igualdad
    E = [1:m]';
end
R = length(d); % restricciones de desigualdad

xk = x0;
%Encontramos una solución x0 en caso que no se proporcione
if(isempty(x0))
    num_var = max(size(A,2),size(D,2));
    xk = linprog(zeros(num_var,1),D,d,A,b);
end

%Resolvemos subproblema cuadrático
%min (1/2)p'Gp + p'g sa Ak*p 
CA = W0(m+1:end)'-m;
CA = sort(CA);

if(isempty(W0))
    CA = find(D*xk==d); %Conjunto activo de desigualdades; W0= A(x0)
end


while ((bandera == 0) &&(iter < maxIter))
    W= [E; (CA+m)]; % Conjunto activo
    W = sort(W);
%     W'
%     CA'
    %Resolvemos subproblema cuadrático
    %min (1/2)p'Gp + p'g sa Ak*p 
    Ak = [A; D(CA,:)]; % matriz de restricciones activas
    r  = length(CA);
    b1 = zeros(m+r,1);   
    g  = G*xk + c;
    if (isempty(Ak))
        % Problema sin restricciones
        p = -G\g;
    else
        p  = esp_nulo( G, Ak, g, b1);
    end
    if (sum(isnan(p))~=0)
        break
    end
    %RAMA 2
    if (max(abs(p))<10^(-9)) %le agregamos una tolerancia
        % calculamos los multiplicadores de Lagrange despejando de las CNPO
        lambda = -(Ak*Ak')\(Ak*g);
        
        if (lambda(m+1:end) >= 0)
            xmin = xk;
            bandera = 1;
            iter = iter-1;
        else
            %Todas las mu [lambdas de desigualdad (m+1:end) menores a cero 
                %y queremos quitar solo una,la mas negativa]
            j = find(lambda == min(lambda(m+1:end)), 1, 'last');
            length(CA)
            j
            sale = CA(j-m)+m;
            CA(j-m) = []; % quitamos la restriccion j
            %Imprimimos mu_j, j
            fprintf("\n%d-ésima iteración:\n",iter)
            fprintf("mu = %.4f ; j = %d \n",min(lambda(m+1:end)),sale)
        end
        
    %RAMA 1
    else  % p != 0
        jna = find(~ismember(1:R,CA)); % restricciones no activas
        
        jna1 = jna(D(jna,:)*p > 0);
        [alpha,ind] = min((d(jna1)-D(jna1,:)*xk) ./ (D(jna1,:)*p));

        if (alpha < 1)
            xk = xk + alpha * p;
            valor_xk = (xk'*G*xk)/2+xk'*c;
            
            CA = [CA; jna1(ind)]; % agregamos la restriccion que se volvió activa
            
            %Imprimimos ||d||_inf, alpha, índice que entra
            fprintf("\n%d-ésima iteración:\n",iter)
            fprintf("||d||_inf = %.4f ; q(x_k+1) = %.4f ; alpha = %.4f ; índice = %d\n",max(abs(p)),valor_xk,alpha,jna1(ind)+m)
        else
            xk = xk + p;
            valor_xk = (xk'*G*xk)/2+xk'*c;
            
            %Imprimimos ||d||_inf, alpha
            fprintf("\n%d-ésima iteración:\n",iter)
            fprintf("||d||_inf = %.4f ; q(x_k+1) = %.4f ; alpha = %d\n",max(abs(p)),valor_xk,alpha)
        end

    end   
    iter = iter + 1; 
    xmin=xk;

end
iter = iter-1;
valor_min = (xmin'*G*xmin)/2+xmin'*c;
end
