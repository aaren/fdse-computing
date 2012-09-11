% Here, we apply boundary conditions to the pressure update
% We will use homogeneous Neumann boundary conditions at all walls
% However, this results in a singular matrix because the Boussinesq equations
% don't depend on the value of P.  To avoid this, set P(2,2)=0

if (P_BC_XMIN==0) % Dirichlet boundary conditions (P @ GXF(2))
% Set the ghost cells to zero, since they won't be used
  MAT(1:NX:NX*NY,1:NX:NX*NY)=eye(NY);
  VEC(1:NX:NX*NY)=zeros(NY,1);
% Set the value at the boundary
% First, set the coefficients of these rows to zero
  MAT(2:NX:NX*NY,:)=0;
  MAT(2:NX:NX*NY,2:NX:NX*NY)=eye(NY);
  VEC(2:NX:NX*NY)=P_BC_XMIN_C1;
elseif (P_BC_XMIN==1) % Neumann boundary conditions (d/dx(P) @ GX(2))
  MAT(1:NX:NX*NY,1:NX:NX*NY)=-1*eye(NY);
  MAT(1:NX:NX*NY,2:NX:NX*NY)=eye(NY);
  VEC(1:NX:NX*NY)=DX(2,:).*P_BC_XMIN_C1(:)';
elseif (P_BC_XMIN==3) % Periodic boundary conditions (P(1)=P(NX-2))
  MAT(1:NX:NX*NY,1:NX:NX*NY)=eye(NY);
  MAT(1:NX:NX*NY,NX-2:NX:NX*NY)=-1*eye(NY);
  VEC(1:NX:NX*NY)=zeros(NY,1);
else
  error('Unknown boundary condition applied to P at XMIN');
end

if (P_BC_XMAX==0) % Dirichlet boundary conditions (P @ GXF(NX-1))
% Set the ghost cells to zero, since they won't be used
  MAT(NX:NX:NX*NY,NX:NX:NX*NY)=eye(NY);
  VEC(NX:NX:NX*NY)=zeros(NY,1);
% Set the value at the boundary
% First, set the coefficients of these rows to zero
  MAT(NX-1:NX:NX*NY,:)=0;
  MAT(NX-1:NX:NX*NY,NX-1:NX:NX*NY)=eye(NY);
  VEC(NX-1:NX:NX*NY)=P_BC_XMAX_C1;
elseif (P_BC_XMAX==1) % Neumann boundary conditions (d/dx(P) @ GX(NX))
  MAT(NX:NX:NX*NY,NX:NX:NX*NY)=eye(NY);
  MAT(NX:NX:NX*NY,NX-1:NX:NX*NY)=-1*eye(NY);
  VEC(NX:NX:NX*NY)=DX(NX,:).*P_BC_XMAX_C1(:)';
elseif (P_BC_XMAX==3) % Periodic boundary conditions (P(NX)=P(3))
  MAT(NX:NX:NX*NY,NX:NX:NX*NY)=eye(NY);
  MAT(NX:NX:NX*NY,3:NX:NX*NY)=-1*eye(NY);
  VEC(NX:NX:NX*NY)=zeros(NY,1);
% P(NX-1)=P(2);
% First, set the coefficients of these rows to zero 
  MAT(NX-1:NX:NX*NY,:)=0;
  MAT(NX-1:NX:NX*NY,NX-1:NX:NX*NY)=eye(NY);
  MAT(NX-1:NX:NX*NY,2:NX:NX*NY)=-1*eye(NY);
  VEC(NX-1:NX:NX*NY)=zeros(NY,1);
else
  error('Unknown boundary condition applied to P at XMAX');
end

if (P_BC_YMIN==0) % Dirichlet boundary conditions
% Set the ghost cells to zero, since they won't be used
  MAT(1:NX,1:NX)=eye(NX);
  VEC(1:NX)=zeros(NX,1);
% Set the value at the boundary (P @ GYF(2))
% First, set the coefficients of these rows to zero
  MAT(NX+1:2*NX,:)=0;
  MAT(NX+1:2*NX,NX+1:2*NX)=eye(NX);
  VEC(NX+1:2*NX)=P_BC_YMIN_C1(:);
elseif (P_BC_YMIN==1) % Neumann boundary conditions (d/dy(P) @ GY(2))
  MAT(1:NX,1:NX)=-1*eye(NX);
  MAT(1:NX,NX+1:2*NX)=eye(NX);
  VEC(1:NX)=DY(1:NX,2).*P_BC_YMIN_C1(1:NX)';
elseif (P_BC_YMIN==3) % Periodic boundary conditions (P(1)=P(NY-2))
  MAT(1:NX,1:NX)=-1*eye(NX);
  MAT(1:NX,NX*(NY-3)+1:NX*(NY-2))=eye(NX);
  VEC(1:NX)=zeros(NX,1);
else
  error('Unknown boundary conditions applied to P at YMIN');
end

if (P_BC_YMAX==0) % Dirichlet boundary conditions (P @ GYF(NY-1))
% Set the ghost cells to zero, since they won't be used
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-1)+1:NX*NY)=eye(NX);
  VEC(NX*(NY-1)+1:NX*NY)=zeros(NX,1);
% Set the value at the boundary
% First, set the coefficients of these rows to zero
  MAT(NX*(NY-2)+1:NX*(NY-1),:)=0;
  MAT(NX*(NY-2)+1:NX*(NY-1),NX*(NY-2)+1:NX*(NY-1))=eye(NX);
  VEC(NX*(NY-2)+1:NX*(NY-1))=P_BC_YMAX_C1(:);
elseif (P_BC_YMAX==1) % Neumann boundary conditions (d/dy(P) @ GY(NY))
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-1)+1:NX*NY)=eye(NX);
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-2)+1:NX*(NY-1))=-1*eye(NX);
  VEC(NX*(NY-1)+1:NX*NY)=DY(1:NX,NY).*P_BC_YMAX_C1(1:NX)';
elseif (P_BC_YMAX==3) % Periodic boundary conditions (P(NY)=P(3))
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-1)+1:NX*NY)=eye(NX);
  MAT(NX*(NY-1)+1:NX*NY,2*NX+1:3*NX)=-eye(NX);
  VEC(NX*(NY-1)+1:NX*NY)=zeros(NX,1);
% P(NY-1)=P(2)
% First, set the coefficients of these rows to zero
  MAT(NX*(NY-2)+1:NX*(NY-1),:)=0;
  MAT(NX*(NY-2)+1:NX*(NY-1),NX*(NY-2)+1:NX*(NY-1))=eye(NX);
  MAT(NX*(NY-2)+1:NX*(NY-1),NX+1:2*NX)=-eye(NX);
  VEC(NX*(NY-2)+1:NX*(NY-1))=zeros(NX,1);
else
  error('Unknown boundary conditions applied to P at YMAX');
end

