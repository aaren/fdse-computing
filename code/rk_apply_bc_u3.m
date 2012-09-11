if (U3_BC_XMIN==0) % Dirichlet boundary conditions (U3 @ GXF(2))
% Set the ghost cells to zero, since they won't be used
  MAT(1:NX:NX*NY,1:NX:NX*NY)=eye(NY);
  VEC(1:NX:NX*NY)=zeros(NY,1);
% Set the value at the boundary
% First, set the coefficents of these rows to zero
  MAT(2:NX:NX*NY,:)=0;
  MAT(2:NX:NX*NY,2:NX:NX*NY)=eye(NY);
  VEC(2:NX:NX*NY)=U3_BC_XMIN_C1(:)';
elseif (U3_BC_XMIN==1) % Neumann boundary conditions (d/dx(U3) @ GX(2))
  MAT(1:NX:NX*NY,1:NX:NX*NY)=-1*eye(NY);
  MAT(1:NX:NX*NY,2:NX:NX*NY)=eye(NY);
  VEC(1:NX:NX*NY)=DX(2,:).*U3_BC_XMIN_C1(:)';
elseif (U3_BC_XMIN==3) % Periodic boundary conditions (U3(1)=U3(NX-2))
  MAT(1:NX:NX*NY,1:NX:NX*NY)=eye(NY);
  MAT(1:NX:NX*NY,NX-2:NX:NX*NY)=-1*eye(NY);
  VEC(1:NX:NX*NY)=zeros(NY,1);
else
  error('Unknown boundary condition applied to U3 at XMIN');
end

if (U3_BC_XMAX==0) % Dirichlet boundary conditions (U3 @ GXF(NX-1))
% Set the ghost cells to zero, since they won't be used
  MAT(NX:NX:NX*NY,NX:NX:NX*NY)=eye(NY);
  VEC(NX:NX:NX*NY)=zeros(NY,1);
% Set the value at the boundary
% First, set the coefficients of these rows to zero
  MAT(NX-1:NX:NX*NY,:)=0;
  MAT(NX-1:NX:NX*NY,NX-1:NX:NX*NY)=eye(NY);
  VEC(NX-1:NX:NX*NY)=U3_BC_XMAX_C1(:)';
elseif (U3_BC_XMAX==1) % Neumann boundary conditions (d/dx(U3) @ GX(NX))
  MAT(NX:NX:NX*NY,NX:NX:NX*NY)=eye(NY);
  MAT(NX:NX:NX*NY,NX-1:NX:NX*NY)=-1*eye(NY);
  VEC(NX:NX:NX*NY)=DX(NX,:).*U3_BC_XMAX_C1(:)';
elseif (U3_BC_XMAX==3) % Periodic boundary conditions (U3(NX)=U3(3))
  MAT(NX:NX:NX*NY,NX:NX:NX*NY)=eye(NY);
  MAT(NX:NX:NX*NY,3:NX:NX*NY)=-1*eye(NY);
  VEC(NX:NX:NX*NY)=zeros(NY,1);
% U3(NX-1)=U3(2);
% First, set these rows to zero
  MAT(NX-1:NX:NX*NY,:)=0;
  MAT(NX-1:NX:NX*NY,NX-1:NX:NX*NY)=eye(NY);
  MAT(NX-1:NX:NX*NY,2:NX:NX*NY)=-1*eye(NY);
  VEC(NX-1:NX:NX*NY)=zeros(NY,1);
else
  error('Unknown boundary condition applied to U3 at XMAX');
end

if (U3_BC_YMIN==0) % Dirichlet boundary conditions
% Set the ghost cells to zero, since they won't be used
  MAT(1:NX,1:NX)=eye(NX);
  VEC(1:NX)=zeros(NX,1);
% Set the value at the boundary (U3 @ GYF(2))
% First, set the coefficients of these rows to zero
  MAT(NX+1:2*NX,:)=0;
  MAT(NX+1:2*NX,NX+1:2*NX)=eye(NX);
  VEC(NX+1:2*NX)=U3_BC_YMIN_C1(:);
elseif (U3_BC_YMIN==1) % Neumann boundary conditions (d/dy(U3) @ GY(2))
  MAT(1:NX,1:NX)=-1*eye(NX);
  MAT(1:NX,NX+1:2*NX)=eye(NX);
  VEC(1:NX)=DY(1:NX,2).*U3_BC_YMIN_C1(1:NX)';
elseif (U3_BC_YMIN==3) % Periodic boundary conditions (U3(1)=U3(NY-2))
  MAT(1:NX,1:NX)=-1*eye(NX);
  MAT(1:NX,NX*(NY-3)+1:NX*(NY-2))=eye(NX); 
  VEC(1:NX)=zeros(NX,1);
else
  error('Unknown boundary conditions applied to U3 at YMIN');
end

if (U3_BC_YMAX==0) % Dirichlet boundary conditions (U3 @ GYF(NY-1))
% Set the ghost cells to zero, since they won't be used
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-1)+1:NX*NY)=eye(NX);
  VEC(NX*(NY-1)+1:NX*NY)=zeros(NX,1);
% Set the value at the boundary
% First, set the coefficients of these rows to zero
  MAT(NX*(NY-2)+1:NX*(NY-1),:)=0;
  MAT(NX*(NY-2)+1:NX*(NY-1),NX*(NY-2)+1:NX*(NY-1))=eye(NX);
  VEC(NX*(NY-2)+1:NX*(NY-1))=U3_BC_YMAX_C1(:)';
elseif (U3_BC_YMAX==1) % Neumann boundary conditions (d/dy(U3) @ GY(NY))
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-1)+1:NX*NY)=eye(NX);
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-2)+1:NX*(NY-1))=-1*eye(NX);
  VEC(NX*(NY-1)+1:NX*NY)=DY(1:NX,NY).*U3_BC_YMAX_C1(1:NX)';
elseif (U3_BC_YMAX==3) % Periodic boundary conditions (U3(NY)=U3(3))
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-1)+1:NX*NY)=eye(NX);
  MAT(NX*(NY-1)+1:NX*NY,2*NX+1:3*NX)=-eye(NX);     
  VEC(NX*(NY-1)+1:NX*NY)=zeros(NX,1);
% U3(NY-1)=U3(2)
% First, set the coefficients of these rows to zero
  MAT(NX*(NY-2)+1:NX*(NY-1),:)=0;
  MAT(NX*(NY-2)+1:NX*(NY-1),NX*(NY-2)+1:NX*(NY-1))=eye(NX);
  MAT(NX*(NY-2)+1:NX*(NY-1),NX+1:2*NX)=-eye(NX);
  VEC(NX*(NY-2)+1:NX*(NY-1))=zeros(NX,1);
else
  error('Unknown boundary conditions applied to U3 at YMAX');
end


