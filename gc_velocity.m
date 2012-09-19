%calculates the velocity of the gravity current

ind=zeros(size(TIME_save));

for t=1:length(TIME_save);
    
    th2=TH_save(ii+1,1,t);
    th1=TH_save(ii,1,t);
    gradth=th2-th1;
    
    [diff,ind(t)]=max(gradth);

    %figure; plot(gradth);
end
x=GXF(ind);
figure; plot(TIME_save,x); xlabel('time','fontsize',14); ylabel('x_n(t)','fontsize',14);

ind1=find(x>=1,1);
ind2=find(x>=3.5,1);

p=polyfit(TIME_save(ind1:ind2),x(ind1:ind2),1);

vel=p(1);


t=ind2;

thH=mean(TH_save(50:75,:,t),1);
ref=thH(80);
indh=find(thH>=(0.8*ref),1);

H=GYF(indh);
