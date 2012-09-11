% This file creates the X and Y coordinate arrays
% The grid arrays are produced as follows:
%   1 First, define the base grids (GX, GY) using on a specified function.
%   2 Then, define the fractional grids (GFX,GFY) halfway between neighboring
%     (GX,GY) grid locations
%   3 Stretch both fractional and base so that GF(2) and GF(N-1) now correspond
%     to the wall locations (e.g. GFX(2)=0, GFX(NY-1)=L)

% Check that NX,NY,NZ,N_TH,LX,LY are all specified:
if (not(exist('NX') & exist('NY') & exist('N_TH') & exist('LX') & exist('LY')))
  error('Error - NX, NY, N_TH, LX, and LY must be set before running create_grid');
end

% ****** User input *******
% Select the type of grid function to be used
% The two options are:
%	GRID_TYPE=1: High resolution at both ends
%	GRID_TYPE=2: High resolution at the center of the domain
GRID_TYPE_X=1;
GRID_TYPE_Y=1;
% Specify the grid stretching parameter, CS.  A value of zero indicates an unstretched grid
CSX=0.0;
CSY=0.0;
% ****** End of User input *******

% If we are using an unstretched grid, add a very small factor to avoid division by zero
CSX=max(CSX,1e-6);
CSY=max(CSY,1e-6);

if (GRID_TYPE_X==1)
  for i=1:NX
    GX(i)=(tanh(CSX*(2*i/NX-1))/tanh(CSX)+1)/2;
  end
elseif (GRID_TYPE_X==2)
  for i=1:NX
    GX(i)=((CSX*(2*i/NX-1)^3+(2*i/NX))/(CSX+1)+1)/2;
  end
else
  error('Unknown GRID_TYPE');
end

if (GRID_TYPE_Y==1)
  for j=1:NY
    GY(j)=(tanh(CSY*(2*j/NY-1))/tanh(CSY)+1)/2;
  end
elseif (GRID_TYPE_Y==2)
  for j=1:NY
    GY(j)=((CSY*(2*j/NY-1)^3+(2*j/NY-1))/(CSY+1)+1)/2;
  end
else
  error('Unknown GRID_TYPE');
end

% Now create the fractional grids, halfway between the GX and GY points
for i=2:NX-1
  GXF(i)=(GX(i)+GX(i+1))/2;
end
for j=2:NY-1
  GYF(j)=(GY(j)+GY(j+1))/2;
end

% Scale both grids so that GXF(2)=0 and GXF(NX-1)=LX (and the same for Y)
GX(1:NX)=GX(1:NX)*LX/(GXF(NX-1)-GXF(2));
GXF(2:NX-1)=GXF(2:NX-1)*LX/(GXF(NX-1)-GXF(2)); 
% Shift the grid to put GXF(2)=0
shift=GXF(2);
GX=GX-shift; GXF=GXF-shift;
GY(1:NY)=GY(1:NY)*LY/(GYF(NY-1)-GYF(2));
GYF(2:NY-1)=GYF(2:NY-1)*LY/(GYF(NY-1)-GYF(2)); 
shift=GYF(2);
% Shift the grid to put GYF(2)=0
GY=GY-shift; GYF=GYF-shift;

% Extrapolate the fractional grid to create ghost cells at i=1,NX, j=1,NY
GXF(1)=2*GXF(2)-GXF(3);
GXF(NX)=2*GXF(NX-1)-GXF(NX-2);
GYF(1)=2*GYF(2)-GYF(3);
GYF(NY)=2*GYF(NY-1)-GYF(NY-2);

% Create matrices for the grid spacing in each direction
% Note, DX,DY are defined at GX,GY gridpoints, and DXF,DYF are defined at GXF,GYF gridpoints
% The size of these matrices should match the state variables, ie U1
% It is assumed that the grid is cartesian, so that DX is constant over the y-direction, etc.
for i=2:NX
  DX(i,1:NY)=(GXF(i)-GXF(i-1));
end
DX(1,1:NY)=DX(2,1:NY); %Ghost cell
for j=2:NY
  DY(1:NX,j)=(GYF(j)-GYF(j-1));
end
DY(1:NX,1)=DY(1:NX,2); %Ghost cell

for i=2:NX-1
  DXF(i,1:NY)=(GX(i+1)-GX(i));
end
DXF(1,1:NY)=DXF(2,1:NY); DXF(NX,1:NY)=DXF(NX-1,1:NY);   % Ghost cells
for j=2:NY-1
  DYF(1:NX,j)=(GY(j+1)-GY(j));
end
DYF(1:NX,1)=DYF(1:NX,2); DYF(1:NX,NY)=DYF(1:NX,NY-1);   % Ghost cells



