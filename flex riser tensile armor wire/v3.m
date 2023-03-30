
clc
close all
import com.comsol.model.*
import com.comsol.model.util.*
ModelUtil.setServerBusyHandler(ServerBusyHandler(2000));
%%
l=1.21/2;


figure
hold on
for gg=10:100;
clc
f=gg*1e3

%%
model = ModelUtil.create('Model');
model.modelPath('C:\Users\Eduardo Becker Groth\Desktop\trabalhos\helice');
model.modelNode.create('comp1');
geometria

j=401;
for i=1:j;
ii=i-1;
z=(l/j)*ii;
model.geom('geom1').feature.create(['pt',sprintf('%01d',i)], 'Point');
model.geom('geom1').feature(['pt',sprintf('%01d',i)]).setIndex('p',0.0775*cos(2*pi*0.83*z)-.0775, 0);
model.geom('geom1').feature(['pt',sprintf('%01d',i)]).setIndex('p',z, 2);
model.geom('geom1').feature(['pt',sprintf('%01d',i)]).setIndex('p',0.0775*sin(2*pi*0.83*z), 1);
model.geom('geom1').feature(['pt',sprintf('%01d',i)]).set('createselection', 'on');
model.geom('geom1').feature(['pt',sprintf('%01d',i)]).set('contributeto', 'csel1');
model.geom('geom1').run(['pt',sprintf('%01d',i)]);
zz(i)=z;
zx(i)=0.0775*cos(2*pi*0.83*z)-.0775;
zy(i)=0.0775*sin(2*pi*0.83*z);

teta(i)=atan((zy(i))/(zx(i)))*2+pi;


end
model.geom('geom1').run;
%plot(teta,'.')

teta2=1.18;



resto

a=mphtable(model,'tbl3');
u=mphtable(model,'tbl_u');
v=mphtable(model,'tbl_v');
w=mphtable(model,'tbl_w');




% for g = 1:10;
% figure
% model.result('pg1').feature('surf1').set('data', 'dset1');
% model.result('pg1').feature('surf1').set('looplevel', sprintf('%01d',g));
% model.result('pg1').feature('surf1').set('expr', 'solid.disp');
% model.result('pg1').run;
% subplot(4,2,[1 3 5 7])
% mphplot(model,'pg1')
% subplot(4,2,2)
% plot(a.data(g,2:i),'.k');
% subplot(4,2,4)
% plot(u.data(g,2:i),'.y');
% subplot(4,2,6)
% plot(v.data(g,2:i),'.g');
% subplot(4,2,8)
% plot(w.data(g,2:i),'.');
% end
% 
% % 
% for g=1:10;
% norm(g,:) = sqrt((w.data(g,2:i).*cos(teta(2:i))).^2+(v.data(g,2:i).*sin(teta(2:i)).^2));
% end
% 
% clear norm
% g=1;
% norm = (u.data(g,2:i).*sin(teta(2:i))+v.data(g,2:i).*cos(teta(2:i)));%.*cos(teta(2:i))+v.data(g,2:i).*sin(teta(2:i)));
% %plot(norm)
for g=1:40;


h1=u.data(g,2:i).*cos(teta(2:i))-v.data(g,2:i).*sin(teta(2:i));%normal
V=u.data(g,2:i).*sin(teta(2:i))+v.data(g,2:i).*cos(teta(2:i));
W=w.data(g,2:i);
h2=W.*cos(pi/6)+V.*cos(pi/3);
h3=V.*sin(pi/3)-W.*sin(pi/6);

if max(h1)>max(h2)
h=h1;
else
h=h2;    
end

if max(h)>max(h3)
    h=h;
else
    h=h3;
end



padx =5000;

ka = 1/((l/401));
%D = fft(fft(L1_U1,padt,1),padx,2);
%D = fft(u.data(g,2:i).*cos(teta(2:i)),padx);

k = 0:ka/padx:ka-ka/padx;
%figure
% %plot(k*2*pi,abs(D))
% 
% figure
% hold on


D = fft((h),padx);
[A(g),P(g)]=max(abs(D(1:2000)));

plot(P(g),u.data(g,1),'o')



end

end
% 
% 
% 
% 
% 
% fff=8;
% for g = fff;
% figure
% model.result('pg1').feature('surf1').set('data', 'dset1');
% model.result('pg1').feature('surf1').set('looplevel', sprintf('%01d',g));
% model.result('pg1').feature('surf1').set('expr', 'solid.disp');
% model.result('pg1').run;
% subplot(4,2,[1 3 5 7])
% mphplot(model,'pg1')
% subplot(4,2,2)
% plot(a.data(g,2:i),'.k');
% subplot(4,2,4)
% plot(u.data(g,2:i),'.y');
% subplot(4,2,6)
% plot(v.data(g,2:i),'.g');
% subplot(4,2,8)
% plot(w.data(g,2:i),'.');
% 
% end
% 
%close all
 for g = 10%:20;
h1=u.data(g,2:i).*cos(teta(2:i))-v.data(g,2:i).*sin(teta(2:i));%normal
V=u.data(g,2:i).*sin(teta(2:i))+v.data(g,2:i).*cos(teta(2:i));
W=w.data(g,2:i);
h2=W.*cos(pi/6);%+V.*cos(pi/3);
h3=V.*sin(pi/3)-W.*sin(pi/6);
% % 
% % % 
% figure 
% hold on
% plot(u.data(g,2:i))
% plot(v.data(g,2:i))
% plot(w.data(g,2:i))
% legend('u','v','w')

D = fft((h3),padx);
%for g = 1:10;
figure
model.result('pg1').feature('surf1').set('data', 'dset1');
model.result('pg1').feature('surf1').set('looplevel', sprintf('%01d',g));
model.result('pg1').feature('surf1').set('expr', 'solid.disp');
model.result('pg1').run;
subplot(4,2,[1 3 5 7])
mphplot(model,'pg1')
subplot(4,2,2)
plot(k,abs(D),'.k');
subplot(4,2,4)
plot(h1,'m');
subplot(4,2,6)
plot(h2,'.g');
subplot(4,2,8)
plot(h3,'.');
end
% 
% % norm = (u.data(g,2:i).*sin(teta(2:i))+v.data(g,2:i).*cos(teta(2:i)));%.*cos(teta(2:i))+v.data(g,2:i).*sin(teta(2:i)));
% % %plot(norm)
% % 
% % padx =5000;
% % hold on
% % ka = 1/((l/201));
% % %D = fft(fft(L1_U1,padt,1),padx,2);
% % D = fft(h2,padx);
% % 
% % k = 0:ka/padx:ka-ka/padx;
% % figure
% % plot(k*2*pi,(D))
% 
% 
% 
% 
