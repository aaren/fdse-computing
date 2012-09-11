if (TH_BC_XMIN(n)==0) % Dirichlet boundary conditions (TH @ GXF(2))
% Set the ghost cells to zero, since they won't be used
  MAT(1:NX:NX*NY,1:NX:NX*NY)=eye(NY);
  VEC(1:NX:NX*NY)=zeros(NY,1);
% Set the value at the boundary
  MAT(2:NX:NX*NY,1:NX:NX*NY)=zeros(NY);
  MAT(2:NX:NX*NY,2:NX:NX*NY)=eye(NY);
  MAT(2:NX:NX*NY,3:NX:NX*NY)=zeros(NY);
  VEC(2:NX:NX*NY)=TH_BC_XMIN_C1(:,n)';
elseif (TH_BC_XMIN(n)==1) % Neumann boundary conditions (d/dx(TH) @ GX(2))
  MAT(1:NX:NX*NY,1:NX:NX*NY)=-1*eye(NY);
  MAT(1:NX:NX*NY,2:NX:NX*NY)=eye(NY);
  VEC(1:NX:NX*NY)=DX(2,:).*TH_BC_XMIN_C1(:,n)';
elseif (TH_BC_XMIN(n)==3) % Periodic boundary conditions (TH(1)=TH(NX-2))
  MAT(1:NX:NX*NY,1:NX:NX*NY)=eye(NY);
  MAT(1:NX:NX*NY,NX-2:NX:NX*NY)=-1*eye(NY);
  VEC(1:NX:NX*NY)=zeros(NY,1);
else
  error('Unknown boundary condition applied to TH at XMIN');
end

if (TH_BC_XMAX(n)==0) % Dirichlet boundary conditions (TH @ GXF(NX-1))
% Set the ghost cells to zero, since they won't be used
  MAT(NX:NX:NX*NY,NX:NX:NX*NY)=eye(NY);
  VEC(NX:NX:NX*NY)=zeros(NY,1);
% Set the value at the boundary
  MAT(NX-1:NX:NX*NY,NX-2:NX:NX*NY)=zeros(NY);
  MAT(NX-1:NX:NX*NY,NX-1:NX:NX*NY)=eye(NY);
  MAT(NX-1:NX:NX*NY,NX:NX:NX*NY)=zeros(NY);
  VEC(NX-1:NX:NX*NY)=TH_BC_XMAX_C1(:,n)';
elseif (TH_BC_XMAX(n)==1) % Neumann boundary conditions (d/dx(TH) @ GX(NX))
  MAT(NX:NX:NX*NY,NX:NX:NX*NY)=eye(NY);
  MAT(NX:NX:NX*NY,NX-1:NX:NX*NY)=-1*eye(NY);
  VEC(NX:NX:NX*NY)=DX(NX,:).*TH_BC_XMAX_C1(:,n)';
elseif (TH_BC_XMAX(n)==3) % Periodic boundary conditions (TH(NX)=TH(3))
  MAT(NX:NX:NX*NY,NX:NX:NX*NY)=eye(NY);
  MAT(NX:NX:NX*NY,3:NX:NX*NY)=-1*eye(NY);
  VEC(NX:NX:NX*NY)=zeros(NY,1);
% TH(NX-1)=TH(2);
% First, set these rows to zero
  MAT(NX-1:NX:NX*NY,:)=0;

  MAT(NX-1:NX:NX*NY,NX-1:NX:NX*NY)=eye(NY);
  MAT(NX-1:NX:NX*NY,2:NX:NX*NY)=-1*eye(NY);
  VEC(NX-1:NX:NX*NY)=zeros(NY,1);
else
  error('Unknown boundary condition applied to TH at XMAX');
end

if (TH_BC_YMIN(n)==0) % Dirichlet boundary conditions
% Set the ghost cells to zero, since they won't be used
  MAT(1:NX,1:NX)=eye(NX);
  VEC(1:NX)=zeros(NX,1);
% Set the value at the boundary (TH @ GYF(2))
  MAT(NX+1:2*NX,NX+1:2*NX)=eye(NX);
  MAT(NX+1:2*NX,1:NX)=zeros(NX);
  MAT(NX+1:2*NX,2*NX+1:3*NX)=zeros(NX);
  VEC(NX+1:2*NX)=TH_BC_YMIN_C1(:,n)';
elseif (TH_BC_YMIN(n)==1) % Neumann boundary conditions (d/dy(TH) @ GY(2))
  MAT(1:NX,1:NX)=-1*eye(NX);
  MAT(1:NX,NX+1:2*NX)=eye(NX);
  VEC(1:NX)=DY(1:NX,2)'.*TH_BC_YMIN_C1(1:NX,n)';
elseif (TH_BC_YMIN(n)==3) % Periodic boundary conditions (TH(1)=TH(NY-2))
  MAT(1:NX,1:NX)=-1*eye(NX);
  MAT(1:NX,NX*(NY-3)+1:NX*(NY-2))=eye(NX); 
  VEC(1:NX)=zeros(NX,1);
else
  error('Unknown boundary conditions applied to TH at YMIN');
end

if (TH_BC_YMAX(n)==0) % Dirichlet boundary conditions (TH @ GYF(NY-1))
% Set the ghost cells to zero, since they won't be used
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-1)+1:NX*NY)=eye(NX);
  VEC(NX*(NY-1)+1:NX*NY)=zeros(NX,1);
% Set the value at the boundary
  MAT(NX*(NY-2)+1:NX*(NY-1),NX*(NY-2)+1:NX*(NY-1))=eye(NX);
  MAT(NX*(NY-2)+1:NX*(NY-1),NX*(NY-1)+1:NX*NY)=zeros(NX);
  MAT(NX*(NY-2)+1:NX*(NY-1),NX*(NY-3)+1:NX*(NY-2))=zeros(NX);
  VEC(NX*(NY-2)+1:NX*(NY-1))=TH_BC_YMAX_C1(:,n)';
elseif (TH_BC_YMAX(n)==1) % Neumann boundary conditions (d/dy(TH) @ GY(NY))
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-1)+1:NX*NY)=eye(NX);
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-2)+1:NX*(NY-1))=-1*eye(NX);
  VEC(NX*(NY-1)+1:NX*NY)=DY(1:NX,NY)'.*TH_BC_YMAX_C1(1:NX,n)';
elseif (TH_BC_YMAX(n)==3) % Periodic boundary conditions (TH(NY)=TH(3))
  MAT(NX*(NY-1)+1:NX*NY,NX*(NY-1)+1:NX*NY)=eye(NX);
  MAT(NX*(NY-1)+1:NX*NY,2*NX+1:3*NX)=-eye(NX);     
  VEC(NX*(NY-1)+1:NX*NY)=zeros(NX,1);
% TH(NY-1)=TH(2)
% First, zero these rows:
  MAT(NX*(NY-2)+1:NX*(NY-1),:)=0;

  MAT(NX*(NY-2)+1:NX*(NY-1),NX*(NY-2)+1:NX*(NY-1))=eye(NX);
  MAT(NX*(NY-2)+1:NX*(NY-1),NX+1:2*NX)=-eye(NX);
  VEC(NX*(NY-2)+1:NX*(NY-1))=zeros(NX,1);
else
  error('Unknown boundary conditions applied to TH at YMAX');
end


