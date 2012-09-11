% This file sets the boundary conditions for diablo
% set_params.m must be first (or NX,NY,N_TH) must be otherwise specified
% There are four walls where we need to apply boundary conditions to the velocity and scalar fields
% The boundary locations are designated XMIN, XMAX, and YMIN, YMAX
% For each boundary there are two parameters to set:
%	 The first, e.g. U_BC_XMIN sets the type of boundary condition:
%		0 = Dirichlet (value specified)
%		1 = Neumann (normal derivative specified)
%		3 = Periodic (cyclic)
%	The second parameter e.g. U_BC_XMIN_C1 specifies either the value or gradient
% 	If periodic boundary conditions are specified, the second parameter is not used
% 	This second parameter is an array allowing the value to vary along the wall
%	The size of U_BC_XMIN_C1, etc. is therefore the size of the transverse direction	
%		e.g. size(U_BC_XMIN_C1)=NY
% Important note: The boundary conditions on the pressure field cannot change in time
% since the pressure solve is currently implemented using an LU decomposition

% Check that NX,NY,N_TH are already set
if (not(exist('NX') & exist('NY') & exist('N_TH')))
  error('Error - NX, NY, and N_TH must be set before running set_params');
end
% Boundary conditions for X-velocity (U1)

% ******* User Input ********

% Boundaries at the edges of the x-coordinates
U1_BC_XMIN=3;	U1_BC_XMIN_C1(1:NY)=0; 
U2_BC_XMIN=3;	U2_BC_XMIN_C1(1:NY)=0;
U3_BC_XMIN=3;	U3_BC_XMIN_C1(1:NY)=0;
P_BC_XMIN=3;	P_BC_XMIN_C1(1:NY)=0;

U1_BC_XMAX=3;	U1_BC_XMAX_C1(1:NY)=0; 
U2_BC_XMAX=3;	U2_BC_XMAX_C1(1:NY)=0;
U3_BC_XMAX=3;	U3_BC_XMAX_C1(1:NY)=0;
P_BC_XMAX=3;	P_BC_XMAX_C1(1:NY)=0;


% Boundary conditions for each scalar
for n=1:N_TH
  TH_BC_XMIN(n)=3;	TH_BC_XMIN_C1(1:NY,n)=0; 
  TH_BC_XMAX(n)=3;	TH_BC_XMAX_C1(1:NY,n)=0; 
end

% Boundaries at the edges of the y-coordiantes
U1_BC_YMIN=1;	U1_BC_YMIN_C1(1:NX)=0; 
U2_BC_YMIN=0;	U2_BC_YMIN_C1(1:NX)=0;
U3_BC_YMIN=1;	U3_BC_YMIN_C1(1:NX)=0;   
P_BC_YMIN=1;	P_BC_YMIN_C1(1:NX)=0; 

U1_BC_YMAX=1;	U1_BC_YMAX_C1(1:NX)=0;
U2_BC_YMAX=0;	U2_BC_YMAX_C1(1:NX)=0; 
U3_BC_YMAX=1;	U3_BC_YMAX_C1(1:NX)=0;
P_BC_YMAX=1;	P_BC_YMAX_C1(1:NX)=0; 

% Boundary conditions for each scalar
for n=1:N_TH
TH_BC_YMIN(n)=1;	TH_BC_YMIN_C1(1:NX,n)=1; 
TH_BC_YMAX(n)=1;	TH_BC_YMAX_C1(1:NX,n)=1; 
end

% ******* End of User Input ********

% Check to make sure that unsupported boundary conditions aren't used
if ((U1_BC_XMIN==1)|(U1_BC_XMAX==1)|(U2_BC_YMIN==1)|(U2_BC_YMAX==1))
  error('Neumann boundary conditions on wall-normal velocity is not supported');
end
