
disp('Creating preconditioner matrices...');

TEMP1=0; % It generally works to create preconditioner matrices using DELTA_T=0
%If this doesn't work well, it might be worth trying a nonzero H_BAR in the
%next line
%TEMP1=(1/Re)*H_BAR(1)/2;

setup.type='ilutp';

ri=zeros((NX-2)*(NY-2)*5,1);
ci=zeros((NX-2)*(NY-2)*5,1);
vi=zeros((NX-2)*(NY-2)*5,1);
MAT=zeros(NX,NY);

% Create some arrays that will be useful for constructing the implicit
% matrices
ii2=2:NXM;
jj2=2:NYM;
ip2=3:NX;
im2=1:NX-2;
jp2=3:NY;
jm2=1:NY-2;
row2=floor((0:(NX-2)*(NY-2)*5-1)/5);
row2=row2+NX+1; % Skip the j=1 rows, which will be handeled with boundary conditions
row2=row2+floor((0:(NX-2)*(NY-2)*5-1)/(5*(NX-2)))+1; % Skip the i=1 rows
row2=row2+floor((0:(NX-2)*(NY-2)*5-1)/(5*(NX-2)))+0; % Skip the i=0 rows

ri(1:(NX-2)*(NY-2)*5,1)=row2';

ci(1:5:(NX-2)*(NY-2)*5,1)=row2(1:5:(NX-2)*(NY-2)*5)';
ci(2:5:(NX-2)*(NY-2)*5,1)=row2(1:5:(NX-2)*(NY-2)*5)'-1;
ci(3:5:(NX-2)*(NY-2)*5,1)=row2(1:5:(NX-2)*(NY-2)*5)'+1;
ci(4:5:(NX-2)*(NY-2)*5,1)=row2(1:5:(NX-2)*(NY-2)*5)'-NX;
ci(5:5:(NX-2)*(NY-2)*5,1)=row2(1:5:(NX-2)*(NY-2)*5)'+NX;

if (N_TH>0)
vi(1:5:(NX-2)*(NY-2)*5,1)=reshape(1-TEMP1*(-1./DX(ii2,jj2)./DXF(ii2,jj2)-1./DX(ip2,jj2)./DXF(ii2,jj2) ...
                          -1./DY(ii2,jp2)./DYF(ii2,jj2)-1./DY(ii2,jj2)./DYF(ii2,jj2)),(NX-2)*(NY-2),1);
vi(2:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DXF(ii2,jj2)./DX(ii2,jj2),(NX-2)*(NY-2),1);
vi(3:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DXF(ii2,jj2)./DX(ip2,jj2),(NX-2)*(NY-2),1);
vi(4:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DY(ii2,jj2)./DYF(ii2,jj2),(NX-2)*(NY-2),1);
vi(5:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DY(ii2,jp2)./DYF(ii2,jj2),(NX-2)*(NY-2),1);
MAT=sparse(ri,ci,vi,NX*NY,NX*NY);
n=1; rk_apply_bc_th;  % Apply boundary conditions (for one of the scalars)
[L_TH,U_TH]=ilu(MAT,setup);
end

vi(1:5:(NX-2)*(NY-2)*5,1)=reshape(1-TEMP1*(-1./DX(ii2,jj2)./DXF(ii2,jj2)-1./DX(ii2,jj2)./DXF(im2,jj2) ...
                          -1./DY(ii2,jp2)./DYF(ii2,jj2)-1./DY(ii2,jj2)./DYF(ii2,jj2)),(NX-2)*(NY-2),1);
vi(2:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DXF(im2,jj2)./DX(ii2,jj2),(NX-2)*(NY-2),1);
vi(3:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DXF(ii2,jj2)./DX(ii2,jj2),(NX-2)*(NY-2),1);
vi(4:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DY(ii2,jj2)./DYF(ii2,jj2),(NX-2)*(NY-2),1);
vi(5:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DY(ii2,jp2)./DYF(ii2,jj2),(NX-2)*(NY-2),1);
MAT=sparse(ri,ci,vi,NX*NY,NX*NY);
rk_apply_bc_u1;
[L_U1,U_U1]=ilu(MAT,setup);


vi(1:5:(NX-2)*(NY-2)*5,1)=reshape(1-TEMP1*(-1./DX(ii2,jj2)./DXF(ii2,jj2)-1./DX(ip2,jj2)./DXF(ii2,jj2) ...
                          -1./DY(ii2,jj2)./DYF(ii2,jm2)-1./DY(ii2,jj2)./DYF(ii2,jj2)),(NX-2)*(NY-2),1);
vi(2:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DXF(ii2,jj2)./DX(ii2,jj2),(NX-2)*(NY-2),1);
vi(3:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DXF(ii2,jj2)./DX(ip2,jj2),(NX-2)*(NY-2),1);
vi(4:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DY(ii2,jj2)./DYF(ii2,jm2),(NX-2)*(NY-2),1);
vi(5:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DY(ii2,jj2)./DYF(ii2,jj2),(NX-2)*(NY-2),1);
MAT=sparse(ri,ci,vi,NX*NY,NX*NY);
rk_apply_bc_u2;
[L_U2,U_U2]=ilu(MAT,setup);

vi(1:5:(NX-2)*(NY-2)*5,1)=reshape(1-TEMP1*(-1./DX(ii2,jj2)./DXF(ii2,jj2)-1./DX(ip2,jj2)./DXF(ii2,jj2) ...
                          -1./DY(ii2,jp2)./DYF(ii2,jj2)-1./DY(ii2,jj2)./DYF(ii2,jj2)),(NX-2)*(NY-2),1);
vi(2:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DXF(ii2,jj2)./DX(ii2,jj2),(NX-2)*(NY-2),1);
vi(3:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DXF(ii2,jj2)./DX(ip2,jj2),(NX-2)*(NY-2),1);
vi(4:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DY(ii2,jj2)./DYF(ii2,jj2),(NX-2)*(NY-2),1);
vi(5:5:(NX-2)*(NY-2)*5,1)=reshape(-TEMP1./DY(ii2,jp2)./DYF(ii2,jj2),(NX-2)*(NY-2),1);
MAT=sparse(ri,ci,vi,NX*NY,NX*NY);
rk_apply_bc_u3;
[L_U3,U_U3]=ilu(MAT,setup);

disp('Performing LU decomposition for pressure Poisson solve...');
vi(1:5:(NX-2)*(NY-2)*5,1)=reshape(-1./DX(ii2,jj2)./DXF(ii2,jj2)-1./DX(ip2,jj2)./DXF(ii2,jj2) ...
                          -1./DY(ii2,jp2)./DYF(ii2,jj2)-1./DY(ii2,jj2)./DYF(ii2,jj2),(NX-2)*(NY-2),1);
vi(2:5:(NX-2)*(NY-2)*5,1)=reshape(1./DXF(ii2,jj2)./DX(ii2,jj2),(NX-2)*(NY-2),1);
vi(3:5:(NX-2)*(NY-2)*5,1)=reshape(1./DXF(ii2,jj2)./DX(ip2,jj2),(NX-2)*(NY-2),1);
vi(4:5:(NX-2)*(NY-2)*5,1)=reshape(1./DY(ii2,jj2)./DYF(ii2,jj2),(NX-2)*(NY-2),1);
vi(5:5:(NX-2)*(NY-2)*5,1)=reshape(1./DY(ii2,jp2)./DYF(ii2,jj2),(NX-2)*(NY-2),1);
MAT=sparse(ri,ci,vi,NX*NY,NX*NY);
rk_apply_bc_p;
[L_P,U_P,p_P]=lu(MAT,'vector'); % For the pressure, use an LU decomposition to solve directly since the coefficient matrix doesn't change
