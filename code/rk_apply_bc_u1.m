if (U1_BC_XMIN==0) % Dirichlet boundary conditions
% Set the ghost cells to zero, since they won't be used
  MAT(1:NX:NX*NY,1:NX:NX*NY)=eye(NY);
  VEC(1:NX:NX*NY)=zeros(NY,1); 
  MAT(2:NX:NX*NY,2:NX:NX*NY)=eye(NY);
  VEC(2:NX:NX*NY)=zeros(NY,1);
% Set the value at the boundary (U1 @ GX(3))
  MAT(3:NX:NX*NY,2:NX:NX*NY)=zeros(NY);
  MAT(3:NX:NX*NY,3:NX:NX*NY)=eye(NY);
  MAT(3:NX:NX*NY,4:NX:NX*NY)=zeros(NY);
  VEC(3:NX:NX*NY)=U1_BC_XMIN_C1;
elseif (U1_BC_XMIN==1) % Neumann boundary conditions (d/dx(U1) @ GXF(2))
% First, the coefficients of these rows to zero
  MAT(2:NX:NX*NY,:)=0;
  MAT(2:NX:NX*NY,2:NX:NX*NY)=-1*eye(NY);
  MAT(2:NX:NX*NY,3:NX:NX*NY)=eye(NY);
  VEC(2:NX:NX*NY)=DXF(2,:).*U1_BC_XMIN_C1(:)';
elseif (U1_BC_XMIN==3) % Periodic boundary conditions (U1(2)=U1(NX-1))
% U1(1)=U1(NX-2)
  MAT(1:NX:NX*NY,1:NX:NX*NY)=eye(NY);
  MAT(1:NX:NX*NY,NX-2:NX:NX*NY)=-1*eye(NY);
  VEC(1:NX:NX*NY)=zeros(NY,1);
% U1(2)=U1(NX-1)
% First, zero out these rows
  MAT(2:NX:NX*NY,:)=0;
  MAT(2:NX:NX*NY,2:NX:NX*NY)=eye(NY);
  MAT(2:NX:NX*NY,NX-1:NX:NX*NY)=-1*eye(NY);
  VEC(2:NX:NX*NY)=zeros(NY,1);
else
  error('Unknown boundary condition applied to U1 at XMIN');
end

if (U1_BC_XMAX==0) % Dirichlet boundary conditions
% Set the ghost cells to zero, since they won't be used
  MAT(NX:NX:NX*NY,NX:NX:NX*NY)=eye(NY);
  VEC(NX:NX:NX*NY)=zeros(NY,1);
% Set the value at the boundary (U1 @ GX(NX-1))
% First, set the coefficients of these rows to zero
  MAT(NX-1:NX:NX*NY,:)=0;
  MAT(NX-1:NX:NX*NY,NX-1:NX:NX*NY)=eye(NY);
  VEC(NX-1:NX:NX*NY)=U1_BC_XMAX_C1;
elseif (U1_BC_XMAX==1) % Neumann boundary conditions (d/dx(U1) @ GXF(NX-1))
  MAT(NX:NX:NX*NY,NX:NX:NX*NY)=eye(NY);
  MAT(NX:NX:NX*NY,NX-1:NX:NX*NY)=-1*eye(NY);
  VEC(NX:NX:NX*NY)=DXF(NX-1,:).*U1_BC_XMAX_C1(:)';
elseif (U1_BC_XMIN==3) % Periodic boundary conditions (U1(NX)=U1(3))
  MAT(NX:NX:NX*NY,NX:NX:NX*NY)=eye(NY);
  MAT(NX:NX:NX*NY,3:NX:NX*NY)=-1*eye(NY);
  VEC(NX:NX:NX*NY)=zeros(NY,1);
else
  error('Unknown boundary condition applied to U1 at XMAX');
end

if (U1_BC_YMIN==0) % Dirichlet boundary conditions
% Set the ghost cells to zero, since they won't be used
  MAT(1:NX,1:NX)=eye(NX);
  VEC(1:NX)=zeros(NX,1);
% Set the value at the boundary (U1 @ GYF(2))
% First, set the coefficients of these rows to zero
  MAT(NX+1:2*NX,:)=0;
  MAT(NX+1:2*NX,NX+1:2*NX)=eye(NX);
  VEC(NX+1:2*NX)=U1_BC_YMIN_C1(:);
elseif (U1_BC_YMIN==1) % Neumann boundary conditions (d/dy(U1) @ GY(2))
  MAT(1:NX,1:NX)=-1*eye(NX);
  MAT(1:NX,NX+1:2*NX)=eye(NX);
  VEC(1:NX)=DY(1:NX,2).*U1_BC_YMIN_C1(1:NX)';
elseif (U1_BC_YMIN==3) % Periodic boundary conditions (U1(1)=U1(NY-2))
  MAT(1:NX,1:NX)=-1*eye(NX);
  MAT(1:NX,NX*(NY-3)+1:NX*(NY-2))=eye(NX); 
  VEC(1:NX)=zeros(NX,1);
else
  error('Unknown boundary conditions applied to U1 at YMIN');
end

if (U1_BC_YMAX==0) % Dirichlet boundary conditions (U1 @ GYF(NY-1))
% Set the ghost cells to zero, since they won't be used
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-1)+1:NX*NY)=eye(NX);
  VEC(NX*(NY-1)+1:NX*NY)=zeros(NX,1);
% Set the value at the boundary
% First, set the coefficients of these rows to zero
  MAT(NX*(NY-2)+1:NX*(NY-1),:)=0;
  MAT(NX*(NY-2)+1:NX*(NY-1),NX*(NY-2)+1:NX*(NY-1))=eye(NX);
  VEC(NX*(NY-2)+1:NX*(NY-1))=U1_BC_YMAX_C1(:);
elseif (U1_BC_YMAX==1) % Neumann boundary conditions (d/dy(U1) @ GY(NY))
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-1)+1:NX*NY)=eye(NX);
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-2)+1:NX*(NY-1))=-1*eye(NX);
  VEC(NX*(NY-1)+1:NX*NY)=DY(1:NX,NY).*U1_BC_YMAX_C1(1:NX)';
elseif (U1_BC_YMAX==3) 
% Periodic boundary conditions (U1(NY)=U1(3))
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-1)+1:NX*NY)=eye(NX);
  MAT(NX*(NY-1)+1:NX*NY,2*NX+1:3*NX)=-eye(NX);     
  VEC(NX*(NY-1)+1:NX*NY)=zeros(NX,1);
% U1(NY-1)=U1(2)
% First, zero out these rows since they will have been set by time-stepping
  MAT(NX*(NY-2)+1:NX*(NY-1),:)=0;
  MAT(NX*(NY-2)+1:NX*(NY-1),NX*(NY-2)+1:NX*(NY-1))=eye(NX);
  MAT(NX*(NY-2)+1:NX*(NY-1),NX+1:2*NX)=-eye(NX);
  VEC(NX*(NY-2)+1:NX*(NY-1))=zeros(NX,1);
else
  error('Unknown boundary conditions applied to U1 at YMAX');
end


