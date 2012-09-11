function d=ddz2(z)
% Second derivative matrix for independent variable z.
% 2nd order centered differences.

% check for equal spacing
if abs(std(diff(z))/mean(diff(z)))>.000001
    disp(['ddz2: values not evenly spaced!'])
    d=NaN;
    return
end

del=z(2)-z(1);N=length(z);

d=zeros(N,N);
for n=2:N-1
    d(n,n-1)=1.;
    d(n,n)=-2.;
    d(n,n+1)=1.;
end
d(1,1)=2;d(1,2)=-5;d(1,3)=4;d(1,4)=-1;
d(N,N)=2;d(N,N-1)=-5;d(N,N-2)=4;d(N,N-3)=-1;
d=d/del^2;
return
end
