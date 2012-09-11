function [sig,w,c]=BCVT(z,C,k,Pr,Ra,iBC1,iBCN,imode)
%
% USAGE: [sig,w]=BCVT(z,C,k,Pr,Ra,iBC1,iBCN,imode)
%
% Stability analysis for 2D bioconvection
% with a constant vertical velocity, U0=1 (non-dimensionalized)
% INPUTS:
% z = vertical coordinate vector (evenly spaced)
% C = concentration profile
% k = horizontal wavenumber 
% Pr = Prandtl number, nu/kappa
% Ra = Rayleigh number
% iBC1 = boundary conditions at z=z(0)
%    (1) velocity: 1=rigid (default), 0=frictionless
%    (2) buoyancy: 1=insulating, 0=fixed-buoyancy (default)
% iBCN = boundary conditions at z=z(N+1).
%        Definitions as for iBC1. Default: iBCN=iBC1.
% imode = mode choice (default imode=1)
%         imode=0: output all modes, sorted by growth rate
%
% OUTPUTS:
% sig = growth rate of FGM
% w = vertical velocity eigenfunction
% c = concentration eigenfunction
%
% CALLS:
% ddz, ddz2, ddz4
%
% Based on code by W. Smyth, OSU, Nov04
% Modified by J.R.Taylor, DAMTP, July12


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stage 1: Preliminaries
%

% check for equal spacing
if abs(std(diff(z))/mean(diff(z)))>.000001
    disp(['SSF: values not evenly spaced!'])
    sig=NaN;
    return
end

% defaults
if nargin<8;iBC1=[1 0];end
if nargin<9;iBCN=iBC1;end
if nargin<10;imode=1;end

% define constants
ii=complex(0.,1.);
del=mean(diff(z));N=length(z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stage 2: Derivative matrices and BCs
%
D1=ddz(z); %1st derivative matrix with 1-sided boundary terms
Cz=D1*C;

D2=ddz2(z); %2nd derivative matrix with 1-sided boundary terms

% Impermeable boundary
D2(1,:)=0;D2(1,1)=-2/del^2;D2(1,2)=1/del^2;
D2(N,:)=0;D2(N,N)=-2/del^2;D2(N,N-1)=1/del^2;

% Asymptotic boundary
% D2(1,:)=0;D2(1,1)=2*(-del*k-1)/del^2;D2(1,2)=2/del^2;
% D2(N,:)=0;D2(N,N)=2*(-del*k-1)/del^2;D2(N,N-1)=2/del^2;

% Fourth derivative
D4=ddz4(z);

% Rigid or frictionless BCs for 4th derivative
D4(1,:)=0;D4(1,1)=(5+2*iBC1(1))/del^4;D4(1,2)=-4/del^4;D4(1,3)=1/del^4;
D4(2,:)=0;D4(2,1)=-4/del^4;D4(2,2)=6/del^4;D4(2,3)=-4/del^4;D4(2,4)=1/del^4;
D4(N,:)=0;D4(N,N)=(5+2*iBCN(1))/del^4;D4(N,N-1)=-4/del^4;D4(N,N-2)=1/del^4;
D4(N-1,:)=0;D4(N-1,N)=-4/del^4;D4(N-1,N-1)=6/del^4;D4(N-1,N-2)=-4/del^4;D4(N-1,N-3)=1/del^4;

% Derivative matrix for buoyancy
D2c=ddz2(z); %2nd derivative matrix with 1-sided boundary terms
% Fixed-buoyancy boundary
D2c(1,:)=0;D2c(1,1)=-2/del^2;D2c(1,2)=1/del^2;
D2c(N,:)=0;D2c(N,N)=-2/del^2;D2c(N,N-1)=1/del^2;
% Insulating boundaries
if iBC1(2)==1;
    D2c(1,:)=0;D2c(1,1)=-2/(3*del^2);D2c(1,2)=2/(3*del^2);
end
if iBCN(2)==1;
    D2c(N,:)=0;D2c(N,N)=-2/(3*del^2);D2c(N,N-1)=2/(3*del^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stage 3: Assemble stability matrices
%
% Laplacian and squared Laplacian matrices
Id=eye(N);
L=D2-k^2*Id;
Lb=D2c-k^2*Id;
LL=D4-2*k^2*D2+k^4*Id;

N2=2*N;
NP=N+1;
A=zeros(N2,N2);B=zeros(N2,N2);

% assemble matrix A 
A=[L Id*0 ; Id*0 Id];

% Compute submatrices of B using Levi's syntax
b11=Pr*LL;
b21=-diag(Cz);
b12=k^2*Pr*Ra*Id;
b22=Lb-D1;

% assemble matrix B
B=[b11 b12 ; b21 b22];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stage 4: Solve eigenvalue problem and extract results
%
% Solve generalized eigenvalue problem
[v,e]=eig(B,A);sigma=diag(e);


% Sort eigvals
[sr,ind]=sort(real(sigma),1,'descend');
sigma=sigma(ind);
v=v(:,ind);


% Extract the selected mode(s)
if imode==0
    sig=sigma;
    w=v(1:N,:);
    c=v(NP:N2,:);
elseif imode>0
    sig=sigma(imode);
    w=v(1:N,imode);
    c=v(NP:N2,imode);
end

return
