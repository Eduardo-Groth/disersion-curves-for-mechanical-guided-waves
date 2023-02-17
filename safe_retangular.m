%% Modelo criado para rotina safe - Eduardo Becker Groth
clc
clear all
close all

import com.comsol.model.*
import com.comsol.model.util.*

model=ModelUtil.create('Model');
model.modelPath('C:\Program Files\COMSOL\COMSOL44');
model.name('SAFE_arame_de_riser.mph');

%% Prâmetros

model.param.set('C11', '281.8e9[Pa]'); % módulo de Young
model.param.set('C22', 'C11');
model.param.set('C33', 'C11');
model.param.set('C44', 'C66');
model.param.set('C55', 'C66');
model.param.set('C66', '84.3e9[Pa]'); % módulo de cisalhamento
model.param.set('C12', 'C11-2*C66');
model.param.set('C13', 'C12');
model.param.set('C23', 'C12');
model.param.set('rho', '7.932e3[kg/m^3]'); % densidade
model.param.set('VLsteel', 'sqrt(C11/rho)');
model.param.set('VTsteel', 'sqrt(C66/rho)');
model.param.set('f', '100e3 [Hz]');
model.param.set('W', '15e-3'); % base da seção
model.param.set('H', '5e-3'); % altura da seção
model.param.set('FL', '0.0000001'); % filete
%%

model.modelNode.create('comp1');

model.func.create('an1', 'Analytic');
model.func('an1').model('comp1');

%% Geometria
model.geom.create('geom1', 2);
model.geom('geom1').lengthUnit('m');
model.geom('geom1').feature.create('r1', 'Rectangle');
model.geom('geom1').feature.create('fil2', 'Fillet');
model.geom('geom1').feature('r1').set('size', {'W' 'H'});
model.geom('geom1').feature('r1').set('base', 'center');
model.geom('geom1').feature('fil2').set('radius', 'FL');
model.geom('geom1').feature('fil2').selection('point').set('r1(1)', [1 2 3 4]);
model.geom('geom1').run;

%%

model.variable.create('var3');
model.variable('var3').model('comp1');
model.variable('var3').set('freq', 'f');
model.variable('var3').set('omg', '2*pi*freq*1[rad]');
model.variable('var3').set('kx', 'lambda');
model.variable('var3').set('Ta', 'abs(imag(kx)/real(kx))');
model.variable('var3').set('Vph', '2*pi*freq/real(kx)');

model.variable.create('var4');
model.variable('var4').model('comp1');
model.variable('var4').set('Sxy', 'C66*(uy+vx)');
model.variable('var4').set('Syy', 'C22*vy+C12*ux+C23*i*w1');
model.variable('var4').set('Syz', 'C44*(wy+i*v1)');
model.variable('var4').set('Sxx', 'C11*ux+C12*vy+C13*i*w1');
model.variable('var4').set('Szz', 'C13*ux+C23*vy+C33*i*w1');
model.variable('var4').set('Sxz', 'C55*(i*u1+wx)');
model.variable('var4').set('Px', '0.5*real(Sxx*conj(i*omg*u)+Sxy*conj(i*omg*v)+Sxz*conj(i*omg*w))');
model.variable('var4').set('Py', '0.5*real(Sxy*conj(i*omg*u)+Syy*conj(i*omg*v)+Syz*conj(i*omg*w))');
model.variable('var4').set('Pz', '0.5*real(Sxz*conj(i*omg*u)+Syz*conj(i*omg*v)+Szz*conj(i*omg*w))');

%%
model.physics.create('c', 'CoefficientFormPDE', 'geom1');
model.physics('c').identifier('riser');
model.physics('c').field('dimensionless').component({'u' 'v' 'w' 'u1' 'v1' 'w1'});

%%
model.mesh.create('mesh1', 'geom1');
model.mesh('mesh1').feature.create('ftri1', 'FreeTri');
model.mesh('mesh1').feature.create('ref1', 'Refine');

%%
model.view('view1').axis.set('xmin', '-37.203163146972656');
model.view('view1').axis.set('ymin', '-4.766697883605957');
model.view('view1').axis.set('xmax', '41.597747802734375');
model.view('view1').axis.set('ymax', '4.572669506072998');

%%
model.physics('c').feature('cfeq1').set('a',{'-rho*omg^2';'0';'0';'0';'0';'0';'0';'-rho*omg^2';'0';'0';'0';'0';'0';'0';'-rho*omg^2';'0';'0';'0';'0';...
                                             '0';'0';'-rho*omg^2';'0';'0';'0';'0';'0';'0';'-rho*omg^2';'0';'0';'0';'0';'0';'0';'-rho*omg^2'});
model.physics('c').feature('cfeq1').set('c',{'C11' '0' '0' 'C66';
                                             '0' 'C12' 'C66' '0';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             '0' 'C66' 'C12' '0';
                                             'C66' '0' '0' 'C22';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             'C55' '0' '0' 'C44';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             '1' '0' '0' '1';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0';
                                             '1' '0' '0' '1';
                                             '0' '0' '0' '0';
                                             '0' '0' '0' '0'; 
                                             '0' '0' '0' '0'; 
                                             '0' '0' '0' '0'; 
                                             '0' '0' '0' '0'; 
                                             '0' '0' '0' '0';
                                             '1' '0' '0' '1'});
model.physics('c').feature('cfeq1').set('be', {'0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0';  ...
'0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0';  ...
'-i*C13' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '-i*C23'; '0' '0'; '0' '0'; '0' '0';  ...
'-i*C55' '0'; '0' '-i*C44'; '0' '0'; '0' '0'; '0' '0'; '0' '0'});
model.physics('c').feature('cfeq1').set('da', {'1'; '0'; '0'; '-rho*omg^2'; '0'; '0'; '0'; '1'; '0'; '0';  ...
'-rho*omg^2'; '0'; '0'; '0'; '1'; '0'; '0'; '-rho*omg^2'; '-C55'; '0';  ...
'0'; '1'; '0'; '0'; '0'; '-C44'; '0'; '0'; '1'; '0';  ...
'0'; '0'; '-C33'; '0'; '0'; '1'});
model.physics('c').feature('cfeq1').set('al', {'0' '0'; '0' '0';'0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0';  ...
'0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0';  ...
'i*C55' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' '0'; '0' 'i*C44'; '0' '0'; '0' '0'; '0' '0';  ...
'i*C13' '0'; '0' 'i*C23'; '0' '0'; '0' '0'; '0' '0'; '0' '0'});

%%
model.mesh('mesh1').feature('size').set('hauto', 3);
model.mesh('mesh1').run;
%%
model.study.create('std1');
model.study('std1').feature.create('eigv', 'Eigenvalue');
%%
model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature.create('e1', 'Eigenvalue');
%%
model.study('std1').feature('eigv').set('initstudyhide', 'on');
model.study('std1').feature('eigv').set('initsolhide', 'on');
model.study('std1').feature('eigv').set('notstudyhide', 'on');
model.study('std1').feature('eigv').set('notsolhide', 'on');
%%
model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').name('Compile Equations: Eigenvalue');
model.sol('sol1').feature('st1').set('studystep', 'eigv');
model.sol('sol1').feature('v1').set('control', 'eigv');

autovalores=40;
PAR.no_pts=0;
limite=100;

matriz_lambda=zeros(autovalores,limite);
vetor_frequencia=zeros(autovalores,limite);
matriz_lambda_imaginario=zeros(autovalores,limite);
matriz_vph=zeros(autovalores,limite);

for fx = 1:limite
    clc
    disp (['iteracao:' num2str(fx) ' de ' num2str(limite) '.' ])
    incremento_frequencia(fx)=fx*1000;
    model.study('std1').feature('eigv').set('neigs', num2str(sprintf('%1.0f',autovalores)));
    model.study('std1').feature('eigv').set('shift', num2str(PAR.no_pts));
    model.param.set('f', num2str(incremento_frequencia(fx)));
    model.sol('sol1').feature('e1').set('shift', num2str(PAR.no_pts));
    model.sol('sol1').feature('e1').set('neigs', num2str(sprintf('%1.0f',autovalores)));
    %model.sol('sol1').feature('e1').set('neigs', num2str( sprintf('%1.0f',PAR.eig_no)));
    model.sol('sol1').runAll;
 
    %%% READ ALL SOLUTIONS

    
    lambda=mphglobal(model, 'lambda' ,'dataset','dset1', 'outersolnum',1);
    lambda_imaginario=mphglobal(model, 'imag(lambda)' ,'dataset','dset1', 'outersolnum',1);
    SOL.Ta=mphglobal(model, 'Ta' ,'dataset','dset1', 'outersolnum',1);
    vph=mphglobal(model, 'Vph' ,'dataset','dset1', 'outersolnum',1);
    SOL.pd= mpheval(model,{'u','v','w', 'Px','Py','Pz'} ,'dataset','dset1', 'outersolnum',1, 'dataonly', 'on');
    frequencia=mphglobal(model, 'f' ,'dataset','dset1', 'outersolnum',1);
    matriz_lambda(:,fx)=lambda;
    matriz_lambda_imaginario(:,fx)=lambda_imaginario;
    vetor_frequencia(:,fx)=frequencia;
    matriz_vph(:,fx)=vph;
    modos(fx).deslocamentos=SOL.pd;
end
clc
ii=0;
for i=1:limite
    for j=1:autovalores
         if abs(matriz_lambda_imaginario(j,i))>1e-3
            matriz_lambda(j,i)=NaN;
         end
         if matriz_lambda(j,i)<=0
           matriz_lambda(j,i)=NaN;
         end
        if matriz_lambda(j,i)>0
        ii=ii+1;
        a(ii,1)=i;
        a(ii,2)=j;
        end
    end
end


%%%%%%% dispersion curves
h=figure; 
hold on
title('Curvas de dispersão','fontsize',14)
set(gca,'fontsize', 14)
for i=1:limite
    plot(matriz_lambda(:,i),vetor_frequencia(:,i),'o')
    %axis([-0.1 0.5 0 10])
xlabel('k [1/m]','fontsize',14)
ylabel('f [Hz]','fontsize',14)
end
axis([0 350 0 10e4])
%saveas(h,'dispersao_quadrada.png');
%%%%%% fase velocity
vph=(vetor_frequencia./matriz_lambda).*2*pi;
figure, hold on
title('Velocidade de fase','fontsize',14)
set(gca,'fontsize', 14)
for i=1:limite
plot(vetor_frequencia(:,i),vph(:,i),'*k')
ylabel('Velocidade de fase [m/s]','fontsize',14)
xlabel('Frequência [Hz]','fontsize',14)
axis([0 1e5 0 9000])
end
clc

%% modal displacement visualization 
% the plot will do four modes of the middle of the range that you set - for see mor modes chage the ii and explore the 'a' matrix  
s=25;
dat=mpheval(model,'u');
tri= dat.t' + 1; 
x= dat.p( 1, :);
y= dat.p( 2, :);

for ii = size(a)/2:size(a)/2+4
fh = figure();
fh.WindowState = 'maximized';
subplot(3,1,1)
trimesh(tri,x,y,real(modos(a(ii,1)).deslocamentos.d3(a(ii,2),:)));title('u_x') 
axis off
subplot(3,1,2)
set(gca,'fontsize',1)
trimesh(tri,x,y,real(modos(a(ii,1)).deslocamentos.d2(a(ii,2),:)));title('u_y')
str={[' frequency ' num2str(real(vetor_frequencia(a(ii,2),a(ii,1)))/1000)  ' khs '], [' wavenumber ' num2str(real(matriz_lambda(a(ii,2),a(ii,1)))) ' 2pi/m ']}
annotation('textbox', [0.02,0.5,0,0],'String',str,'FitBoxToText','on');
axis off
subplot(3,1,3)
set(gca,'fontsize',1)
trimesh(tri,x,y,real(modos(a(ii,1)).deslocamentos.d1(a(ii,2),:)));title('u_z')
axis off
end

