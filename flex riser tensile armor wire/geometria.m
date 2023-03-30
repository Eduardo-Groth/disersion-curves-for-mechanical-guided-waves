model.geom.create('geom1', 3);
model.modelNode('comp1').name('wire_tensile_armor');
model.geom('geom1').feature.create('wp1', 'WorkPlane');
model.geom('geom1').feature('wp1').set('unite', true);
model.geom('geom1').feature.create('pc1', 'ParametricCurve');
model.geom('geom1').feature('wp1').geom.feature.create('r1', 'Rectangle');
model.geom('geom1').feature('wp1').geom.feature('r1').set('base', 'center');
model.geom('geom1').feature('pc1').setIndex('coord', 's', 2);
model.geom('geom1').feature('wp1').geom.feature('r1').setIndex('pos', '0', 0);
model.geom('geom1').feature('wp1').geom.feature('r1').setIndex('pos', '0', 1);
model.geom('geom1').feature('pc1').setIndex('coord', '0.0775*cos(2*pi*0.83*s)-.0775', 0);
model.geom('geom1').feature('pc1').setIndex('coord', '0.0775*sin(2*pi*0.83*s)', 1);
model.geom('geom1').feature('pc1').set('axistype', 'cartesian');
model.geom('geom1').run('pc1');
model.geom('geom1').feature('pc1').set('parmax', l);
model.geom('geom1').feature.create('swe1', 'Sweep');
model.geom('geom1').feature('swe1').set('twistcomp', 'off');
model.geom('geom1').feature('swe1').set('includefinal', false);
model.geom('geom1').feature('swe1').set('maxknots', '1000000');
model.geom('geom1').feature('swe1').selection('face').set('wp1', [1]);
model.geom('geom1').feature('swe1').selection('edge').set('pc1', [1]);


model.geom('geom1').feature('wp1').geom.feature.create('fil1', 'Fillet');
model.geom('geom1').feature('wp1').geom.feature('fil1').selection('point').set('r1', [1 2 3 4]);
model.geom('geom1').feature('wp1').geom.feature('fil1').set('radius', '1.5e-3');
model.geom('geom1').feature('wp1').set('planetype', 'coordinates');
model.geom('geom1').feature('wp1').geom.feature('r1').setIndex('pos', '0', 0);
model.geom('geom1').feature('pc1').setIndex('pos', '0', 0);
model.geom('geom1').feature('pc1').setIndex('pos', '0', 1);
model.geom('geom1').feature('wp1').setIndex('genpoints', '1', 1, 0);
model.geom('geom1').feature('wp1').setIndex('genpoints', '0', 1, 1);

model.geom('geom1').feature('wp1').setIndex('genpoints', '-cos(1.2)', 2, 2);


model.geom('geom1').feature('wp1').setIndex('genpoints', '0', 2, 0);
model.geom('geom1').feature('wp1').geom.feature('r1').setIndex('size', '5e-3', 0);
model.geom('geom1').feature('wp1').geom.feature('r1').setIndex('size', '15e-3', 1);
model.geom('geom1').feature('wp1').setIndex('genpoints', 'sin(pi/3)', 2, 1);



model.geom('geom1').selection.create('csel1', 'CumulativeSelection');
model.geom('geom1').selection('csel1').name('Cumulative Selection 1');

