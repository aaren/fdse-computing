% Use this script to add user-defined terms to the RHS
% of the momentum and scalar equations.
% Momentum forcing should be added to F1, F2, F3
% Scalar sources/sinks should be added to FTH

% ****** For Bioconvection ********
% This section adds a constant vertical velocity to the scalar equation
% The added term should be evaluated on the TH-grid (see grid diagram
% in docs folder), and should be added directly to FTH(ii,jj,n)
% For Bioconvection, uncomment the following lines
%U0=1; % Add a constant upward swimming 
%FTH(ii,jj,1)=FTH(ii,jj,1)-U0*((TH(ii,jp,1)-TH(ii,jj,1))./DY(ii,jp)+(TH(ii,jj,1)-TH(ii,jm,1))./DY(ii,jj))/2;

% ****** For internal wave generation *****
%amp=0.1; % Internal wave amplitude
%omega=1/sqrt(2);   % Forcing frequency (N=2*pi)
%width=0.01; % Size of forcing region (Gaussian width)
%for i=2:NX-1 
%for j=3:NY-1 % Loop over GY points within the physical domain
%  F2(i,j)=F2(i,j)+amp*exp(-((GXF(i)-LX/2).^2+(GY(j)-LY/2).^2)/(2*width^2))*sin(omega*TIME);
%end
%end
