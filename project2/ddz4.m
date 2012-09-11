
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=ddz4(z)
% Compute fourth derivative matrix for independent variable vector z.

% check for equal spacing
if abs(std(diff(z))/mean(diff(z)))>.000001
    disp(['ddz2: values not evenly spaced!'])
    d=NaN;
    return
end

del=z(2)-z(1);N=length(z);

d=zeros(N,N);
for n=3:N-2
    d(n,n-2)=1.;
    d(n,n-1)=-4.;
    d(n,n)=6.;
    d(n,n+1)=-4.;
    d(n,n+2)=1.;
end
% assume f(0)=f(N+1)=f''(0)=f''(N+1)=0 (1-sided BC's not needed)
d(1,1)=5;d(1,2)=-4;d(1,3)=1.;
d(2,1)=-4.;d(2,2)=6.;d(2,3)=-4.;d(2,4)=1.;
d(N-1,N)=-4.;d(N-1,N-1)=6.;d(N-1,N-2)=-4.;d(N-1,N-3)=1.;
d(N,N)=5;d(N,N-1)=-4;d(N,N-2)=1;
d=d/del^4;
return
end