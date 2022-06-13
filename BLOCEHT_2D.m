                         %% %%-%%-%% PROBLEME DE BLOCHET 2D %%-%%-%% %%
clear all
clc
close all 
%% PARAMETRES D'ENTREE :
H=2;             % Inclinaison des plans.
N=37;            % Nb de descritisation suivant x.
M=90;            % Nb de descritisation suivant z.
DX=1/N;          % pas suivant x : ∆x=1/N.
DZ=1/M;          % pas suivant z : ∆z=1/N.
B=20;            % Largeur du Blochet.
L=20;            % Longuer du Blochet.
omega=1.9;       % Poid de surrelaxation ω ∈ ]1 ÷ 2[.

%% EPAISSEUR DE FILM & CONSTANTS DE CALCUL

for i=2:N

    X(i)=1/N*(i-1);                     % Espacement,suivant x, entre deux noeuds successives. 
 
        h(i)=H-(H-1)*X(i);                  % h(i)
        hm(i)=H-(H-1)*(X(i)-0.5*DX);        % h(i-0.5)
        hp(i)=H-(H-1)*(X(i)+0.5*DX);        % h(i+0.5)
        hm1(i)=H-(H-1)*(X(i)-DX);           % h(i-1)
        hp1(i)=H-(H-1)*(X(i)+DX);           % h(i+1)

        hcub(i)=h(i).^3;                    % h(i)^3
        hmcub(i)=hm(i).^3;                  % h(i-0.5)^3
        hpcub(i)=hp(i).^3;                  % h(i+0.5)^3     

        Const1(i)=L^2/B^2;                  % (L/B)^2
        Const2(i)=(hpcub(i)+hmcub(i)+2*Const1(i)*hcub(i)*(DX/DZ)^2);
        AA(i)=hpcub(i)/Const2(i);                        % A(i)
        BB(i)=hmcub(i)/Const2(i);                        % B(i)
        CC(i)=(hcub(i).*Const1(i)*(DX/DZ)^2)/Const2(i);  % C(i)
        DD(i)=(3*DX*(hp1(i)-hm1(i)))./Const2(i);         % D(i)
end
%% CALCUL DE PRESSION :
tic
%-------------Inesialisation de pression--------------
p=zeros(N+1,M+1);
%-------------   Calcule de pression   ---------------

ITER=1000;              % Nombre max des itérations.
GAMA=1e-7;              % Marge d'erreur γ=10^-6
for K = 1:ITER
sumij=0.0;

    for i=2:N    
     for j=2:M 

        p(i,j)=[AA(i)*p(i+1,j)+BB(i)*p(i-1,j)+CC(i)*(p(i,j+1)+p(i,j-1))-DD(i)]*omega+(1-omega)*p(i,j);

        sumij=sumij+p(i,j);     
        sumijN=sumij+p(i,M-1);
     end
    end
    
%------------ Vérification de convergance -------------

sum(K+1)=sumij;
Pourcentage = abs(sum(K+1)-sum(K))/abs(sum(K+1));

if Pourcentage < GAMA
break
end

end

%-----------------     Résutats     -------------------

t=toc              % afficahge de temps 
Nb_it=K;           % Nombre des itérations
p_max=max(max(p))  % Maximume de champ pression
if not (B==L)
disp(['Comparaison avec la résultat 1D : p_1D = 0,2570                    B/L=1 et H=2'])
end
F=sumij*DZ*DZ      % Force de portance
if B==L
disp(['Comparaison avec la résultat analytique 2D : F_anal = 0,0689'])
else 
disp(['Comparaison avec la résultat numérique  1D : F_1D = 0,1588         B/L=1 et H=2'])
end 
%% CALCUL DES DEBITS SUR x ET z 

% ---------------Débit Qx---------------- 
sumQpx=0.0;
sumQcx=0.0;

for j=1:M+1

sumQpx=sumQpx + (DZ/DX)*((H-(H-1)*((1/N*(N-1))-0.5*DX))^3)*p(N-1,j)/12;
sumQcx=sumQcx + DZ*(H-(H-1)*((1/N*(N-1))-0.5*DX))/2;

end
Qx=sumQpx+sumQcx
if B==L
disp(['Comparaison avec la résultat analytique 2D : Qx_anal = 0.8473'])
else
disp(['Comparaison avec la résultat numérique  1D : Qx_1D = 0,6669     B/L=1 et H=2'])
end
% ---------------Débit Qz----------------

sumQpz=0.0;

for i=1:N+1

sumQpz=sumQpz+(DX/DZ)*((H-(H-1)*((1/M*(M-1))-0.5*DX))^3)*p(i,M-1)/12;

end
Qz=sumQpz
if B==L
disp(['Comparaison avec la résultat analytique 2D : Qz_anal = 0.2470'])
else
disp(['Comparaison avec la résultat numérique  1D : Qz_1D = 0             B/L=1 et H=2'])
end
%% CALCUL DE FORCE DE FROTTEMENT Ff
sumFf=0.0;
    for i=2:N   
     for j=2:M/2 


sumFf=sumFf+2*( 1/hp(i) + hp(i)/2*( (p(i+1,j+1)+p(i,j+1)) - (p(i+1,j)+p(i,j)) )/(4*DX) )*DX*DZ;

     end
    end
Ff=sumFf
Pf=Ff

%% PRESENTATION GRAPHIQUE DU CHAMP PRESSION  

figure
x=linspace(0,1,N+1);
y=linspace(0,1,M+1);
[Y,X]=meshgrid(x,y);
xx=transpose(X);
yy=transpose(Y);
surf(xx,yy,p);
title('Répartition de champ pression 2D pour un modèle de Blochet')
zlabel(' ̅p [-]' );
xlabel(' ̅z [-]');
ylabel(' ̅x [-]');
set(gca,'FontSize',14);
set(gca,'FontWeight','normal');
set(gca,'FontName','Times New Roman');
colormap jet;
colorbar;