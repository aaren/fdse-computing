% This script directly applies the boundary conditions to the velocity
% fields

if (U1_BC_XMIN==0) % Dirichlet boundary conditions
% Set the ghost cells to zero, since they won't be used
  U1(1,:)=0; 
  U1(2,:)=0;
% Set the value at the boundary (U1 @ GX(3))
  U1(3,:)=U1_BC_XMIN_C1(:)';
elseif (U1_BC_XMIN==1) % Neumann boundary conditions (d/dx(U1) @ GXF(2))
  U1(2,:)=U1(3,:)-DXF(2,:).*U1_BC_XMIN_C1(:)';
elseif (U1_BC_XMIN==3) % Periodic boundary conditions (U1(2)=U1(NX-2))
  U1(1,:)=U1(NX-2,:); % Note, this row isn't solved for
  U1(2,:)=U1(NX-1,:);
else
  error('Unknown boundary condition applied to U1 at XMIN');
end

if (U1_BC_XMAX==0) % Dirichlet boundary conditions
% Set the ghost cells to zero, since they won't be used
  U1(NX,:)=0;
% Set the value at the boundary (U1 @ GX(NX-1))
  U1(NX-1,:)=U1_BC_XMAX_C1(:)';
elseif (U1_BC_XMAX==1) % Neumann boundary conditions (d/dx(U1) @ GXF(NX-1))
 U1(NX,:)=U1(NX-1,:)+DXF(NX-1,:).*U1_BC_XMAX_C1(:)';
elseif (U1_BC_XMIN==3) % Periodic boundary conditions (U1(NX)=U1(3))
  U1(NX,:)=U1(3,:);
else
  error('Unknown boundary condition applied to U1 at XMAX');
end

if (U1_BC_YMIN==0) % Dirichlet boundary conditions
% Set the ghost cells to zero, since they won't be used
  U1(:,1)=0;
% Set the value at the boundary (U1 @ GYF(2))
  U1(:,2)=U1_BC_YMIN_C1(:)';
elseif (U1_BC_YMIN==1) % Neumann boundary conditions (d/dy(U1) @ GY(2))
  U1(:,1)=U1(:,2)-DY(1:NX,2).*U1_BC_YMIN_C1(1:NX)';
elseif (U1_BC_YMIN==3) % Periodic boundary conditions (U1(1)=U1(NY-2))
  U1(:,1)=U1(:,NY-2);
else
  error('Unknown boundary conditions applied to U1 at YMIN');
end

if (U1_BC_YMAX==0) % Dirichlet boundary conditions (U1 @ GYF(NY-1))
% Set the ghost cells to zero, since they won't be used
  U1(:,NY)=0;
% Set the value at the boundary
  U1(:,NY-1)=U1_BC_YMAX_C1(:)';
elseif (U1_BC_YMAX==1) % Neumann boundary conditions (d/dy(U1) @ GY(NY))
  U1(:,NY)=U1(:,NY-1)+DY(1:NX,NY).*U1_BC_YMAX_C1(1:NX)';
elseif (U1_BC_YMAX==3) % Periodic boundary conditions (U1(NY)=U1(3))
  U1(:,NY)=U1(:,3);
  U1(:,NY-1)=U1(:,2); % U1(NY-1)=U1(2)
else
  error('Unknown boundary conditions applied to U1 at YMAX');
end

% This cell (j=1) just here to keep the grids aligned properly
% Always set this value equal to zero
U2(:,1)=0;

if (U2_BC_XMIN==0) % Dirichlet boundary conditions (U2 @ GXF(2))
% Set the ghost cells to zero, since they won't be used
  U2(1,:)=zeros(NY,1); 
% Set the value at the boundary
  U2(2,:)=U2_BC_XMIN_C1(:)';
elseif (U2_BC_XMIN==1) % Neumann boundary conditions (d/dx(U2) @ GX(2))
  U2(1,:)=U2(2,:)-DX(2,:).*U2_BC_XMIN_C1(:)';
elseif (U2_BC_XMIN==3) % Periodic boundary conditions (U2(1)=U2(NX-2))
  U2(1,:)=U2(NX-2,:);
else
  error('Unknown boundary condition applied to U2 at XMIN');
end

if (U2_BC_XMAX==0) % Dirichlet boundary conditions (U2 @ GXF(NX-1))
% Set the ghost cells to zero, since they won't be used
  U2(NX,:)=0;
% Set the value at the boundary
  U2(NX-1,:)=U2_BC_XMAX_C1(:)';
elseif (U2_BC_XMAX==1) % Neumann boundary conditions (d/dx(U2) @ GX(NX))
  U2(NX,:)=U2(NX-1,:)+DX(NX,:).*U2_BC_XMAX_C1(:)';
elseif (U2_BC_XMIN==3) % Periodic boundary conditions (U2(NX)=U2(3))
  U2(NX,:)=U2(3,:);
  U2(NX-1,:)=U2(2,:);
else
  error('Unknown boundary condition applied to U2 at XMAX');
end

if (U2_BC_YMIN==0) % Dirichlet boundary conditions (U2 @ GY(3))
% Set the ghost cells to zero, since they won't be used
  U2(:,2)=zeros(NX,1);
% Set the value at the boundary
  U2(:,3)=U2_BC_YMIN_C1(:)';
elseif (U2_BC_YMIN==1) % Neumann boundary conditions (d/dy(U2) @ GYF(2))
  U2(:,2)=U2(:,3)-DYF(1:NX,2).*U2_BC_YMIN_C1(1:NX)';
elseif (U2_BC_YMIN==3) % Periodic boundary conditions (U2(2)=U2(NY-1))
  U2(:,1)=U2(:,NY-2); 
  U2(:,2)=U2(:,NY-1);
else
  error('Unknown boundary conditions applied to U2 at YMIN');
end
 
if (U2_BC_YMAX==0) % Dirichlet boundary conditions (U2 @ GY(NY-1))
% Set the ghost cells to zero, since they won't be used
  U2(:,NY)=0;
% Set the value at the boundary
  U2(:,NY-1)=U2_BC_YMAX_C1(:)';
elseif (U2_BC_YMAX==1) % Neumann boundary conditions (d/dy(U2) @ GYF(NY-1))
  U2(:,NY)=U2(:,NY-1)+DYF(1:NX,NY-1).*U2_BC_YMAX_C1(1:NX)';
elseif (U2_BC_YMAX==3) % Periodic boundary conditions (U2(NY)=U1(3))
  U2(:,NY)=U2(:,3);
else
  error('Unknown boundary conditions applied to U2 at YMAX');
end

if (U3_BC_XMIN==0) % Dirichlet boundary conditions (U3 @ GXF(2))
% Set the ghost cells to zero, since they won't be used
  U3(1,:)=0;
% Set the value at the boundary
  U3(2,:)=U3_BC_XMIN_C1(:)';
elseif (U3_BC_XMIN==1) % Neumann boundary conditions (d/dx(U3) @ GX(2))
  U3(1,:)=U3(2,:)-DX(2,:).*U3_BC_XMIN_C1(:)';
elseif (U3_BC_XMIN==3) % Periodic boundary conditions (U3(1)=U3(NX-2))
  U3(1,:)=U3(NX-2,:);
else
  error('Unknown boundary condition applied to U3 at XMIN');
end

if (U3_BC_XMAX==0) % Dirichlet boundary conditions (U3 @ GXF(NX-1))
% Set the ghost cells to zero, since they won't be used
  U3(NX,:)=0;
% Set the value at the boundary
  U3(NX-1,:)=U3_BC_XMAX_C1(:)';
elseif (U3_BC_XMAX==1) % Neumann boundary conditions (d/dx(U3) @ GX(NX))
  U3(NX,:)=U3(NX-1,:)+DX(NX,:).*U3_BC_XMAX_C1(:)';
elseif (U3_BC_XMIN==3) % Periodic boundary conditions (U3(NX)=U3(3))
  U3(NX,:)=U3(3,:);
  U3(NX-1,:)=U3(2,:);
else
  error('Unknown boundary condition applied to U3 at XMAX');
end

if (U3_BC_YMIN==0) % Dirichlet boundary conditions
% Set the ghost cells to zero, since they won't be used
  U3(:,1)=0;
% Set the value at the boundary (U3 @ GYF(2))
  U3(:,2)=U3_BC_YMIN_C1(:)';
elseif (U3_BC_YMIN==1) % Neumann boundary conditions (d/dy(U3) @ GY(2))
  U3(:,1)=U3(:,2)-DY(1:NX,2).*U3_BC_YMIN_C1(1:NX)';
elseif (U3_BC_YMIN==3) % Periodic boundary conditions (U3(1)=U3(NY-2))
  U3(:,1)=U3(:,NY-2);
else
  error('Unknown boundary conditions applied to U3 at YMIN');
end

if (U3_BC_YMAX==0) % Dirichlet boundary conditions (U3 @ GYF(NY-1))
% Set the ghost cells to zero, since they won't be used
  U3(:,NY)=0;
% Set the value at the boundary
  U3(:,NY-1)=U3_BC_YMAX_C1(:)';
elseif (U3_BC_YMAX==1) % Neumann boundary conditions (d/dy(U3) @ GY(NY))
  U3(:,NY)=U3(:,NY-1)+DY(1:NX,NY).*U3_BC_YMAX_C1(1:NX)';
elseif (U3_BC_YMAX==3) % Periodic boundary conditions (U3(NY)=U3(3))
  U3(:,NY)=U3(:,3);
  U3(:,NY-1)=U3(:,2); % U3(NY-1)=U3(2)
else
  error('Unknown boundary conditions applied to U3 at YMAX');
end




