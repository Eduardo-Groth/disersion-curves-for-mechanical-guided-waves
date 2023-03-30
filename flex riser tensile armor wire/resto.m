model.material.create('mat1');
model.material('mat1').propertyGroup('def').set('density', {'7860'});
model.material('mat1').propertyGroup('def').set('youngsmodulus', {'200e9'});
model.material('mat1').propertyGroup('def').set('poissonsratio', {'0.3'});

model.mesh.create('mesh1', 'geom1');
model.mesh('mesh1').autoMeshSize(3);
model.mesh('mesh1').run;

model.physics.create('solid', 'SolidMechanics', 'geom1');
model.physics('solid').feature.create('pc1', 'PeriodicCondition',2);
model.physics('solid').feature('pc1').selection.set([1 10]);
model.physics('solid').feature.create('fix1', 'Fixed', 2);
model.physics('solid').feature('fix1').selection.set([1 10]);

model.study.create('std1');
model.study('std1').feature.create('eig', 'Eigenfrequency');
model.study('std1').feature('eig').activate('solid', true);
model.study('std1').feature('eig').set('neigs', 40);
model.study('std1').feature('eig').set('shift', sprintf('%04d',f));

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'eig');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'eig');
model.sol('sol1').feature.create('e1', 'Eigenvalue');
model.sol('sol1').feature('e1').set('shift',  sprintf('%04d',f));
model.sol('sol1').feature('e1').set('neigs',10);
model.sol('sol1').feature('e1').set('transform', 'eigenfrequency');
model.sol('sol1').feature('e1').set('control', 'eig');
model.sol('sol1').attach('std1');

model.result.create('pg1', 3);
model.result('pg1').set('data', 'dset1');
model.result('pg1').feature.create('surf1', 'Surface');
%model.result('pg1').feature('surf1').set('expr', {'solid.disp'});
model.result('pg1').name('desplazamiento');
model.result('pg1').feature('surf1').feature.create('def', 'Deform');
% %model.result('pg1').feature('surf1').feature('def').set('expr', {'u' 'v' 'w'});
% model.result('pg1').feature('surf1').feature('def').set('descr', 'Displacement field (Material)');

% model.result.create('pg2', 3);
% model.result('pg2').set('data', 'dset1');
% model.result('pg2').feature.create('surf1', 'Surface');
% %model.result('pg1').feature('surf1').set('expr', {'solid.disp'});
% model.result('pg2').name('desplazamiento2');
% model.result('pg2').feature('surf1').feature.create('def', 'Deform');
% model.result('pg2').feature('surf1').feature('def').set('expr', {'u' 'v' 'w'});
% %model.result('pg1').feature('surf1').feature('def').set('descr', 'Displacement field (Material)');

model.sol('sol1').runAll;

model.result('pg1').run;
% for i=1:2;
% 
% model.result.numerical.create(['pev',sprintf('%01d',i)], 'EvalPoint');
% model.result.numerical(['pev',sprintf('%01d',i)]).selection.set([i]);
% model.result.table.create(['tbl',sprintf('%01d',i)], 'Table');
% model.result.numerical(['pev',sprintf('%01d',i)]).set('expr', 'u');
% model.result.table(['tbl',sprintf('%01d',i)]).comments('mae natureza');
% model.result.numerical(['pev',sprintf('%01d',i)]).set('table',['tbl',sprintf('%01d',i)]);
% model.result.numerical(['pev',sprintf('%01d',i)]).setResult;
% 
% end
%%
model.result.numerical.create('pev3', 'EvalPoint');
model.result.numerical('pev3').selection.named('geom1_csel1_pnt');
model.result.table.create('tbl3', 'Table');
model.result.numerical('pev3').set('expr', 'solid.disp');
model.result.table('tbl3').comments('Point Evaluation 3 (solid.disp)');
model.result.numerical('pev3').set('table', 'tbl3');
model.result.numerical('pev3').setResult;
model.result.table.create('tbl_u', 'Table');
model.result.numerical('pev3').set('expr', 'u');
model.result.numerical('pev3').set('descr', 'Displacement field, X component');
model.result.numerical('pev3').set('table', 'tbl_u');
model.result.numerical('pev3').setResult;

model.result.table.create('tbl_v', 'Table');
model.result.numerical('pev3').set('expr', 'v');
model.result.numerical('pev3').set('descr', 'Displacement field, Y component');
model.result.numerical('pev3').set('table', 'tbl_v');
model.result.numerical('pev3').setResult;

model.result.table.create('tbl_w', 'Table');
model.result.numerical('pev3').set('expr', 'w');
model.result.numerical('pev3').set('descr', 'Displacement field, Z component');
model.result.numerical('pev3').set('table', 'tbl_w');
model.result.numerical('pev3').setResult;