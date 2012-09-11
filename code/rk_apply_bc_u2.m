% U2(j=1) is never used since U2(j=2) acts as the ghost cell
% This is just here to keep the grids aligned properly
% Always set this value equal to zero
MAT(1:NX,1:NX)=eye(NX);
VEC(1:NX)=zeros(NX,1);

if (U2_BC_XMIN==0) % Dirichlet boundary conditions (U2 @ GXF(2))
% Set the ghost cells to zero, since they won't be used
  MAT(1:NX:NX*NY,1:NX:NX*NY)=eye(NY);
  VEC(1:NX:NX*NY)=zeros(NY,1); 
% Set the value at the boundary
% First, set the coefficients of these rows to zero
  MAT(2:NX:NX*NY,:)=0;
  MAT(2:NX:NX*NY,2:NX:NX*NY)=eye(NY);
  VEC(2:NX:NX*NY)=U2_BC_XMIN_C1(:);
elseif (U2_BC_XMIN==1) % Neumann boundary conditions (d/dx(U2) @ GX(2))
  MAT(1:NX:NX*NY,1:NX:NX*NY)=-1*eye(NY);
  MAT(1:NX:NX*NY,2:NX:NX*NY)=eye(NY);
  VEC(1:NX:NX*NY)=DX(2,:).*U2_BC_XMIN_C1(:)';
elseif (U2_BC_XMIN==3) % Periodic boundary conditions (U2(1)=U2(NX-2))
  MAT(1:NX:NX*NY,1:NX:NX*NY)=eye(NY);
  MAT(1:NX:NX*NY,NX-2:NX:NX*NY)=-1*eye(NY);
  VEC(1:NX:NX*NY)=zeros(NY,1);
else
  error('Unknown boundary condition applied to U2 at XMIN');
end

if (U2_BC_XMAX==0) % Dirichlet boundary conditions (U2 @ GXF(NX-1))
% Set the ghost cells to zero, since they won't be used
  MAT(NX:NX:NX*NY,NX:NX:NX*NY)=eye(NY);
  VEC(NX:NX:NX*NY)=zeros(NY,1);
% Set the value at the boundary
% First, set the coefficients of these rows to zero
  MAT(NX-1:NX:NX*NY,:)=0;
  MAT(NX-1:NX:NX*NY,NX-1:NX:NX*NY)=eye(NY);
  VEC(NX-1:NX:NX*NY)=U2_BC_XMAX_C1;
elseif (U2_BC_XMAX==1) % Neumann boundary conditions (d/dx(U2) @ GX(NX))
  MAT(NX:NX:NX*NY,NX:NX:NX*NY)=eye(NY);
  MAT(NX:NX:NX*NY,NX-1:NX:NX*NY)=-1*eye(NY);
  VEC(NX:NX:NX*NY)=DX(NX,:).*U2_BC_XMAX_C1(:)';
elseif (U2_BC_XMIN==3) % Periodic boundary conditions (U2(NX)=U2(3))
  MAT(NX:NX:NX*NY,NX:NX:NX*NY)=eye(NY);
  MAT(NX:NX:NX*NY,3:NX:NX*NY)=-1*eye(NY);
  VEC(NX:NX:NX*NY)=zeros(NY,1);
% U2(NX-1)=U2(2);
% First, set the coefficients of these rows to zero
  MAT(NX-1:NX:NX*NY,:)=0;
  MAT(NX-1:NX:NX*NY,NX-1:NX:NX*NY)=eye(NY);
  MAT(NX-1:NX:NX*NY,2:NX:NX*NY)=-1*eye(NY);
  VEC(NX-1:NX:NX*NY)=zeros(NY,1);
else
  error('Unknown boundary condition applied to U2 at XMAX');
end

if (U2_BC_YMIN==0) % Dirichlet boundary conditions (U2 @ GY(3))
% Set the ghost cells to zero, since they won't be used
  MAT(NX+1:2*NX,:)=0;
  MAT(NX+1:2*NX,NX+1:2*NX)=eye(NX);
  VEC(NX+1:2*NX)=zeros(NX,1);
% Set the value at the boundary
  MAT(2*NX+1:3*NX,:)=0;
  MAT(2*NX+1:3*NX,2*NX+1:3*NX)=eye(NX);
  VEC(2*NX+1:3*NX)=U2_BC_YMIN_C1(:);
elseif (U2_BC_YMIN==1) % Neumann boundary conditions (d/dy(U2) @ GYF(2))
  MAT(2*NX+1:3*NX,:)=0;
  MAT(2*NX+1:3*NX,2*NX+1:3*NX)=-1*eye(NX);
  MAT(2*NX+1:3*NX,3*NX+1:4*NX)=eye(NX);
  VEC(2*NX+1:3*NX)=DYF(1:NX,2).*U2_BC_YMIN_C1(1:NX)';
elseif (U2_BC_YMIN==3) % Periodic boundary conditions (U2(2)=U2(NY-1))
% U2(1)=U2(NY-2)
  MAT(1:NX,1:NX)=eye(NX);
  MAT(1:NX,NX*(NY-3)+1:NX*(NY-2))=-1*eye(NX);
  VEC(1:NX)=zeros(NX,1);
% U2(2)=U2(NY-1)
% First, zero out these rows in the coefficient matrix
  MAT(NX+1:2*NX,:)=0;
  MAT(NX+1:2*NX,NX+1:2*NX)=eye(NX);
  MAT(NX+1:2*NX,NX*(NY-2)+1:NX*(NY-1))=-1*eye(NX); 
  VEC(NX+1:2*NX)=zeros(NX,1);
else
  error('Unknown boundary conditions applied to U2 at YMIN');
end

if (U2_BC_YMAX==0) % Dirichlet boundary conditions (U2 @ GY(NY-1))
% Set the ghost cells to zero, since they won't be used
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-1)+1:NX*NY)=eye(NX);
  VEC(NX*(NY-1)+1:NX*NY)=zeros(NX,1);
% Set the value at the boundary
  MAT(NX*(NY-2)+1:NX*(NY-1),:)=0;
  MAT(NX*(NY-2)+1:NX*(NY-1),NX*(NY-2)+1:NX*(NY-1))=eye(NX);
  VEC(NX*(NY-2)+1:NX*(NY-1))=U2_BC_YMAX_C1(:);
elseif (U2_BC_YMAX==1) % Neumann boundary conditions (d/dy(U2) @ GYF(NY-1))
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-1)+1:NX*NY)=eye(NX);
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-2)+1:NX*(NY-1))=-1*eye(NX);
  VEC(NX*(NY-1)+1:NX*NY)=DYF(1:NX,NY-1).*U2_BC_YMAX_C1(1:NX)';
elseif (U2_BC_YMAX==3) % Periodic boundary conditions (U2(NY)=U2(3))
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-1)+1:NX*NY)=eye(NX);
  MAT(NX*(NY-1)+1:NX*NY,2*NX+1:3*NX)=-eye(NX);     
  VEC(NX*(NY-1)+1:NX*NY)=zeros(NX,1);
else
  error('Unknown boundary conditions applied to U2 at YMAX');
end


