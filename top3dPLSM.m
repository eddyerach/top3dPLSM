%Usage example top3dPLSM(12,6,6,0.3)
function top3dPLSM(nelx,nely,nelz,volfrac)
%clear;clc;clf;
%nelx = 24;nely = 16;nelz = 8;volfrac = 0.3;
%%  Inicializacion de la funcion level set
r = nely*0.1;%RADIUS OF INITIAL HOLES
hX = repmat(nelx*[repmat([1/6,5/6],1,3),repmat([0,1/3,2/3,1],1,2),1/2],[1,3]);
hY = repmat(nely*[kron([0,1/2,1],ones(1,2)),kron([1/4,3/4],ones(1,4)),1/2],[1,3]);
hZ = nelz*kron([0,1/2,1],ones(1,15));
[X,Y,Z] = meshgrid(-1:2/nelx:1,-1:2/nely:1,-1:2/nelz:1);%experimentar con el tamano de paso
[x,y,z] = meshgrid(-1:2/(nelx-1):1,-1:2/(nely-1):1,-1:2/(nelz-1):1);%experimentar con el tamano de paso
dX = bsxfun(@minus,repmat(X,[1,1,1,numel(hX)]),reshape(hX,1,1,1,numel(hX)));
dY = bsxfun(@minus,repmat(Y,[1,1,1,numel(hY)]),reshape(hY,1,1,1,numel(hY)));
Phi0 = zeros(nely+1,nelx+1,nelz+1);
Phi0(2:end-1,2:end-1,2:end-1) = repmat(meshgrid(0:1:nelx-2,0:1:nely-2),1,1,nelz-1);
%Phi0(2:end-1,2:end-1,2:end-1) = 2e-3*ones(nely-1,nelx-1,nelz-1);
Phi1 = max(-3,min(3,min(sqrt(dX.^2+dY.^2)-r,[],4)));
Phi = min(Phi0,Phi1);
Phi = Phi0;
%%  Inicializacion de la funcion de base radial
cRBF = 1e-5;
nNode = (nelx +1)*(nely+1)*(nelz+1);
% 
% Ax = abs(bsxfun(@minus, X(:),X(:)'));
% Ay = abs(bsxfun(@minus, Y(:),Y(:)'));
% Az = abs(bsxfun(@minus, Z(:),Z(:)'));

Ax = (bsxfun(@minus, X(:),X(:)'));
Ay = (bsxfun(@minus, Y(:),Y(:)'));
Az = (bsxfun(@minus, Z(:),Z(:)'));
% Ax = X(:) - X(:)';
% Ay = Y(:) - Y(:)';
% Az = Z(:) - Z(:)';
%cRBF = rand;
A = sqrt(Ax.^2+Ay.^2+Az.^2+cRBF^2); % Multi cuadratica %correcion (-) por (+)
%A = (Ax.^3+Ay.^3+Az.^3);
%r = (sqrt(Ax.^2+Ay.^2+Az.^2+cRBF^2));
%A = (Ax+Ay+Az).*log(Ax+Ay+Az);
%A = abs(Ax+Ay+Az);
%A = ((Ax+Ay+Az).^2).*log(abs(Ax+Ay+Az)); %Thin Plate Spline. No funciona.
%A = exp(-1*((Ax+Ay+Az).^2)); % Gaussiana No funciona
G = [A,ones(nNode,1),X(:),Y(:),Z(:);[ones(1,nNode);X(:)';Y(:)';Z(:)'],zeros(4,4)];
pGpX = [Ax./A,repmat([0,1,0,0],nNode,1);repmat([0;1;0;0],1,nNode),zeros(4,4)];
pGpY = [Ay./A,repmat([0,0,1,0],nNode,1);repmat([0;0;1;0],1,nNode),zeros(4,4)];
pGpZ = [Az./A,repmat([0,0,0,1],nNode,1);repmat([0;0;0;1],1,nNode),zeros(4,4)];
%pG = gradient(G,[X,Y,Z]);
Alpha = G\[Phi(:);0;0;0;0];
%Alpha = inv(G)*[Phi(:);0;0;0;0];
%%  Definicion de condiciones de borde Liu Tovar
%USER-DEFINED LOAD DOFs
il = nelx; jl = nely/2; kl = 0:nelz;                         % Coordinates*
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs*
loaddof = 3*loadnid(:) - 1;                             % DOFs* termino fuente
% USER-DEFINED SUPPORT FIXED DOFs
[jf,kf] = meshgrid(1:nely+1,1:nelz+1);                  % Coordinates*
fixednid = (kf-1)*(nely+1)*(nelx+1)+jf;                 % Node IDs*
fixeddofs = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs*
%% PREPARE FINITE ELEMENT ANALYSIS
E0 = 1; Emin = 1e-9; nu = 0.3; %MATERIAL PROPERTIES
KE = lk_H8(nu);
eleNode = auxiliar(nelx,nely,nelz);
edofMat = kron(eleNode,[3,3,3])+repmat([-2,-1,0],nelx*nely*nelz,8);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nelx*nely*nelz,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nelx*nely*nelz,1);
%  Defincion de condiciones de borde
% F = sparse(3*((nely+1)*(nelx+1)*((nelz+2)/2))-(3*(nely/2))-1,1,-100,3*nNode,1);%NODAL LOADS
% fixeddofs = retrieveFixeddofs(nelx,nely,nelz);
% freedofs = setdiff(1:3*nNode,fixeddofs);
% U = zeros(3*nNode,1);
nele = nelx*nely*nelz;%
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);%
F = sparse(loaddof,1,-1,ndof,1);%
U = zeros(ndof,1);%
freedofs = setdiff(1:ndof,fixeddofs);%
%%  Iteracion de optimizacion
nLoop = 160; nRelax = 100; %dt = 0.008
dt = 0.01; delta = 10; mu = 20; gamma = 0.05;
comp = zeros(nLoop,1); vol = zeros(nLoop,1); 
sumAlpha = zeros(nLoop,1);
sumAlphaX = zeros(nLoop,1);
sumAlphaY = zeros(nLoop,1);
sumAlphaZ = zeros(nLoop,1);
for iT = 1:nLoop
 %% FINITE ELEMENT ANALYSIS
  [s,t,u] = meshgrid(-1:0.1:1,-1:0.1:1,-1:0.1:1);
  tmpPhi = (1-s(:)).*(1-t(:)).*(1-u(:))/8*Phi(eleNode(:,1))'+...
        (1+s(:)).*(1-t(:)).*(1-u(:))/8*Phi(eleNode(:,2))'+...
        (1+s(:)).*(1+t(:)).*(1-u(:))/8*Phi(eleNode(:,3))'+...
        (1-s(:)).*(1+t(:)).*(1-u(:))/8*Phi(eleNode(:,4))'+...
        (1-s(:)).*(1-t(:)).*(1+u(:))/8*Phi(eleNode(:,5))'+...
        (1+s(:)).*(1-t(:)).*(1+u(:))/8*Phi(eleNode(:,6))'+...
        (1+s(:)).*(1+t(:)).*(1+u(:))/8*Phi(eleNode(:,7))'+...
        (1-s(:)).*(1+t(:)).*(1+u(:))/8*Phi(eleNode(:,8))';
  eleVol = sum(tmpPhi>=0,1)'/numel(s);
  %eleVol = sum(tmpPhi>=0,1)'/8;
  vol(iT) = sum(eleVol)/(nelx*nely*nelz);
  sK = reshape(KE(:)*(Emin+eleVol'*(E0-Emin)),24*24*nelx*nely*nelz,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs,1) = K(freedofs,freedofs)\F(freedofs,1);
  %aux1 = sum((U(edofMat)*KE).*U(edofMat),2);
  eleComp = sum((U(edofMat)*KE).*U(edofMat),2).*(Emin+eleVol*(E0-Emin));
  comp(iT) = sum(eleComp);
  %% RESULTS
  fprintf('No.%i, Obj:%f, Vol:%f\n',[iT,comp(iT),vol(iT)]);
  figure(1); subplot(2,1,1); plot(comp(1:iT),'-'); title('Compliance'); 
  subplot(2,1,2); plot(vol(1:iT),'-'); title('Volume fraction'); 
  figure(2);clf;isosurface(X,Y,Z,Phi,0);
  
 %% CONVERGENCE CHECK
  if iT>nRelax && abs(vol(iT)-volfrac)/volfrac<1e-3 && ...
    all(abs(comp(iT)-comp(iT-9:iT-1))/comp(iT)<1e-3)
    break;
  end
 %% LAGRANGE MULTIPLIER
  if iT<=nRelax
    lag = mu*(vol(iT)-vol(1)+(vol(1)-volfrac)*iT/nRelax);
  else
    lag = lag+gamma*(vol(iT)-volfrac);
    gamma = min(gamma+0.05,5);
  end
 %% LEVEL SET FUNCTION EVOLUTION
  gradG = (sqrt((pGpX*Alpha).^2+(pGpY*Alpha).^2+(pGpZ*Alpha).^2));
  indexDelta = (abs(Phi(:))<=delta);
  DeltaPhi = zeros(size(Phi));
  DeltaPhi(indexDelta) = (0.75/delta)*(1-(Phi(indexDelta).^2)/delta^2);
  
  eleComp = reshape(eleComp,nely,nelx,nelz);
  eleCompLR = ([eleComp(:,1,:), eleComp] + [eleComp, eleComp(:,end,:)]);
  eleCompLRUD = ([eleCompLR;eleCompLR(end,:,:)]+[eleCompLR(1,:,:);eleCompLR]);
  nodeComp = (cat(3,eleCompLRUD(:,:,1),eleCompLRUD)+cat(3,eleCompLRUD,eleCompLRUD(:,:,end)))/8;
  
  %eleComp = reshape(eleComp,nely,nelx,nelz);
  
%   eleCompLRUD = ([eleComp;eleComp(end,:,:)]+[eleComp(1,:,:);eleComp]);
%   eleCompLR = ([eleCompLRUD(:,1,:), eleCompLRUD] + [eleCompLRUD, eleCompLRUD(:,end,:)]);
%   nodeComp = (cat(3,eleCompLR(:,:,1),eleCompLR)+cat(3,eleCompLR,eleCompLR(:,:,end)))/8;

  
%M = scatteredInterpolant(x(:),y(:),z(:),eleComp(:),'linear','linear');  
%  nodeComp = M(X,Y,Z);
  %nodeComp = (nodeComp1 + nodeComp)./2;
  %nodeComp1 = reshape(nodeComp1,nely+1,nelx+1,nelz+1);
  B = (nodeComp(:)/median(nodeComp(:))-lag).*DeltaPhi(:)*delta/0.75; %original
  %B = (nodeComp(:)/median(nodeComp(:))-lag);
  B = [B;0;0;0;0];
  %B = B.*(abs(gradG.*Alpha));
  %B = (nodeComp(:)-lag).*DeltaPhi(:)*delta/0.75;
  %figure(3);clf;isosurface(X,Y,Z,nodeComp);
  Alpha = Alpha+dt*(G\B);
  sumAlpha(iT) = sum(Alpha);
  sumAlphaX(iT)=sum(Alpha.*[X(:);0;0;0;0]);
  sumAlphaY(iT)=sum(Alpha.*[Y(:);0;0;0;0]);
  sumAlphaZ(iT)=sum(Alpha.*[Z(:);0;0;0;0]);
  
  figure(3); subplot(4,1,1);plot(sumAlpha(1:iT),'-'); title('Sum(Alpha)');
  subplot(4,1,2); plot(sumAlphaX(1:iT),'-'); title('Sum(Alpha*X)'); 
  subplot(4,1,3); plot(sumAlphaY(1:iT),'-'); title('Sum(Alpha*Y)'); 
  subplot(4,1,4); plot(sumAlphaZ(1:iT),'-'); title('Sum(Alpha*Z)'); 
  
  %Alpha = Alpha/mean(gradG(unique(eleNode((eleVol<1 & eleVol>0),:))));
  if any(eleVol<1 & eleVol>0)
      Alpha = Alpha/mean(gradG(unique(eleNode((eleVol<1&eleVol>0),:))));
  end
  Phi = reshape(G(1:end-4,:)*Alpha,nely+1,nelx+1,nelz+1);
end
PhiExtended = -100.*ones(nely+3,nelx+3,nelz+3);
PhiExtended(2:end-1,2:end-1,2:end-1) = Phi;
figure(4);clf,isosurface(PhiExtended,0);
end
%plotSurface(X,Y,Z,Phi);
%clf;
%patch(fv);isosurface(X,Y,Z,Phi,0);
%%  Funciones auxiliares
function [cont] = auxiliar(nelx,nely,nelz)
cont = zeros(nelx*nely*nelz,8);
numElem = 1;
for z1 = 0:nelz-1
    for x1 = 0:nelx-1
        for y1 = nely-1:-1:0
            nidz = (nelx+1)*(nely+1);
            cont(numElem,4) = (z1*(nelx+1)*(nely+1))+(x1*(nely+1))+(nely+1-y1);
            cont(numElem,3) = cont(numElem,4) + (nely+1);
            cont(numElem,2) = cont(numElem,4) + nely;
            cont(numElem,1) = cont(numElem,4) - 1;
            cont(numElem,5) = cont(numElem,1) + nidz;
            cont(numElem,7) = cont(numElem,3) + nidz;
            cont(numElem,6) = cont(numElem,2) + nidz;
            cont(numElem,8) = cont(numElem,4) + nidz;
            numElem = numElem + 1;    
        end
    end
end
end
function [KE] = lk_H8(nu)
A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
    -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
k = 1/144*A'*[1; nu];

K1 = [k(1) k(2) k(2) k(3) k(5) k(5);
    k(2) k(1) k(2) k(4) k(6) k(7);
    k(2) k(2) k(1) k(4) k(7) k(6);
    k(3) k(4) k(4) k(1) k(8) k(8);
    k(5) k(6) k(7) k(8) k(1) k(2);
    k(5) k(7) k(6) k(8) k(2) k(1)];
K2 = [k(9)  k(8)  k(12) k(6)  k(4)  k(7);
    k(8)  k(9)  k(12) k(5)  k(3)  k(5);
    k(10) k(10) k(13) k(7)  k(4)  k(6);
    k(6)  k(5)  k(11) k(9)  k(2)  k(10);
    k(4)  k(3)  k(5)  k(2)  k(9)  k(12)
    k(11) k(4)  k(6)  k(12) k(10) k(13)];
K3 = [k(6)  k(7)  k(4)  k(9)  k(12) k(8);
    k(7)  k(6)  k(4)  k(10) k(13) k(10);
    k(5)  k(5)  k(3)  k(8)  k(12) k(9);
    k(9)  k(10) k(2)  k(6)  k(11) k(5);
    k(12) k(13) k(10) k(11) k(6)  k(4);
    k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
K4 = [k(14) k(11) k(11) k(13) k(10) k(10);
    k(11) k(14) k(11) k(12) k(9)  k(8);
    k(11) k(11) k(14) k(12) k(8)  k(9);
    k(13) k(12) k(12) k(14) k(7)  k(7);
    k(10) k(9)  k(8)  k(7)  k(14) k(11);
    k(10) k(8)  k(9)  k(7)  k(11) k(14)];
K5 = [k(1) k(2)  k(8)  k(3) k(5)  k(4);
    k(2) k(1)  k(8)  k(4) k(6)  k(11);
    k(8) k(8)  k(1)  k(5) k(11) k(6);
    k(3) k(4)  k(5)  k(1) k(8)  k(2);
    k(5) k(6)  k(11) k(8) k(1)  k(8);
    k(4) k(11) k(6)  k(2) k(8)  k(1)];
K6 = [k(14) k(11) k(7)  k(13) k(10) k(12);
    k(11) k(14) k(7)  k(12) k(9)  k(2);
    k(7)  k(7)  k(14) k(10) k(2)  k(9);
    k(13) k(12) k(10) k(14) k(7)  k(11);
    k(10) k(9)  k(2)  k(7)  k(14) k(7);
    k(12) k(2)  k(9)  k(11) k(7)  k(14)];
KE = 1/((nu+1)*(1-2*nu))*...
    [ K1  K2  K3  K4;
    K2'  K5  K6  K3';
    K3' K6  K5' K2';
    K4  K3  K2  K1'];
end
function [aux] = retrieveFixeddofs(nelx,nely,nelz)
acu = 1;
aux = zeros(1,(nely+1)*(nelz+1));
pivote = 3*(nely+1)*(nelx+1);
for z1 = 0:(nelz)
    for y1 = 1:3*(nely+1)
        aux(1,acu) = y1 + pivote * z1;
        acu = acu + 1;
    end
end
end

