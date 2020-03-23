function [x, cout] = tv_o(y, a, lam, E)
% Entrées :
% - y : image à améliorer (type : double)
% - a : facteur numérique pour éviter l'annulation du terme gradient
% - lam : \lambda, multiplicateur de lagrange
% - E : couronne entourant la portion de l'image où on réalise l'inpainting
%       elle représente la zone bruitée entourant la région à réparer 

% Sorties :
% - x : image après une application de l'algo du gradient
% - coût : énergie J(x), évaluée pour constater sa décroissance au fil
%          des itérations

% Version de travail : il serait intéressant d'étudier une implémentation 
% en C++

[Nx, Ny] = size(y);
N = numel(y);
x = zeros(Nx,Ny);
grad_cout = 0;

for i=2:Nx-1
    for j=2:Ny-1
        u_O = y(i,j);
        lambda_e = lam*E(i,j);
        
        % Voisins
        u_E = y(i,j+1);
        u_N = y(i+1,j);
        u_W = y(i,j-1);
        u_S = y(i-1,j);
        
        u_NE = y(i-1,j+1);
        u_NW = y(i-1,j-1);
        u_SE = y(i+1,j+1);
        u_SW = y(i+1,j-1);
        
        % Gradients des voisins
        delta_e = sqrt((u_E-u_O)^2 + (u_NE + u_N - u_S - u_SE)^2/16);
        delta_w = sqrt((u_W-u_O)^2 + (u_NW + u_N - u_S - u_SW)^2/16);
        delta_n = sqrt((u_N-u_O)^2 + (u_E + u_NE - u_NW - u_W)^2/16);
        delta_s = sqrt((u_S-u_O)^2 + (u_E + u_SE - u_SW - u_W)^2/16);
        
        % Poids
        we = 1/(sqrt(delta_e^2 + a^2));
        ww = 1/(sqrt(delta_w^2 + a^2));
        wn = 1/(sqrt(delta_n^2 + a^2));
        ws = 1/(sqrt(delta_s^2 + a^2));
        
        % Calcul du gradient du point (pour le coût cf. infra)
        grad_cout = grad_cout + (1/we + 1/ww + 1/wn + 1/ws)/4;
        
        % Coefficients h
        h_OE = we/(we+ww+wn+ws + lambda_e);
        h_ON = wn/(we+ww+wn+ws + lambda_e);
        h_OW = ww/(we+ww+wn+ws + lambda_e);
        h_OS = ws/(we+ww+wn+ws + lambda_e);
        h_OO = lambda_e/(we + ww + wn + ws + lambda_e);
        
        % Calcul u_0
        u = h_OE*u_E + h_ON*u_N + h_OW*u_W + h_OS*u_S + h_OO*u_O;
        x(i,j) = u;
    end
end

% Fonctionnelle discrète qu'on évalue pour constater la décroissance 
% de l'énergie. sum(...,'all') ne fonctionne que pour Matlab > 2018b donc
% on pourra utiliser à la place :

% cout = grad_cout + lam/2*sum((y-x).^2, 'all');
cout = grad_cout + lam/2*sum(sum((y-x).^2));
cout = cout/N;

end

