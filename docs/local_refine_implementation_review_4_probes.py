from pathlib import Path
import ast
import runpy
import json
import numpy as np
import pytest
import shapely
from fvcom_mesh_tools import patch
from fvcom_mesh_tools.refine import limit_rfactor
from fvcom_mesh_tools.io.fvcom_native import (export_fvcom_case, read_fvcom_case, read_obc_types, apply_obc_depth_control)
from fvcom_mesh_tools.qa import QACheck, run_qa
from fvcom_mesh_tools.sizing import _geometry_from_file
H = runpy.run_path('tests/test_patch.py')
M = runpy.run_path('tests/test_local_refine_driver.py')['square_mesh']

@pytest.mark.parametrize('ident', ['0', 0.5, -1])
def test_malformed_element_is_not_inherited(ident):
    c = QACheck('probe','quality',True,False,'test','test',1,
                offenders=[{'kind':'element','id':ident}])
    assert patch.introduced_violations([c], 2, M().elements)


def test_inherited_real_qa_truncation():
    m = M()
    # Both 45-degree triangles are unchanged and retained, but fail 46 degrees.
    qa = run_qa(m, min_angle_deg=46, max_offenders=1)
    c = next(c for c in qa.checks if c.check_id == 'c1_min_angle')
    assert c.n_violations == 2 and len(c.offenders) == 1
    assert patch.introduced_violations([c], 2, m.elements) == []


def test_conflicts_matches_actual_field():
    xy, tri = H['grid_mesh'](15,15)
    regions = [(shapely.box(490,690,510,710),5,1000),
               (shapely.box(990,690,1010,710),12,100)]
    p = np.array([[1000.,700.]])
    alone = patch.patch_sizing(xy,tri,[regions[1]],distmesh_scale=1)(p)
    joint = patch.patch_sizing(xy,tri,regions,distmesh_scale=1)(p)
    assert np.array_equal(alone,joint) and alone[0] == 12
    assert 'coarse' not in patch.region_conflicts(regions,['fine','coarse'])['finer_than_declared']


def test_thin_gradation_has_no_false_zero():
    rep = patch.field_gradation(lambda p: 100+np.asarray(p)[:,0],
                               shapely.box(0,0,100,20),spacing=10)
    assert rep['n_samples'] == 9
    assert rep['max_slope'] == pytest.approx(1.)


def test_obc_type_survives_depth_control(tmp_path):
    m=M(); m.obc_type=2
    controlled,_=apply_obc_depth_control(m)
    files=export_fvcom_case(controlled,tmp_path,'b',obc_depth_control=False)
    assert read_obc_types(files['obc']) == [2]


def test_grid_nan_rejected(tmp_path):
    files=export_fvcom_case(M(),tmp_path,'b',obc_depth_control=False)
    rows=files['grd'].read_text().splitlines()
    fields=rows[4].split(); fields[1]='nan'; rows[4]=' '.join(fields)
    files['grd'].write_text('\n'.join(rows)+'\n')
    with pytest.raises(ValueError):
        read_fvcom_case(files['grd'],files['dep'],files['obc'])


def test_grid_node_labels_follow_fvcom_row_order(tmp_path):
    files=export_fvcom_case(M(),tmp_path,'b',obc_depth_control=False)
    rows=files['grd'].read_text().splitlines()
    fields=rows[4].split(); fields[0]='2'; rows[4]=' '.join(fields)
    files['grd'].write_text('\n'.join(rows)+'\n')
    back=read_fvcom_case(files['grd'],files['dep'],files['obc'])
    assert np.array_equal(back.nodes,M().nodes)

@pytest.mark.parametrize('bad', [0.,-1.,np.nan,np.inf])
def test_limiter_invalid_depth_refused(bad):
    with pytest.raises(ValueError):
        limit_rfactor(np.array([[0,1,2]]), np.full(3,bad), np.ones(3,bool), .2)

# Controls that should pass.
def test_mapping_direction_permuted_and_deleted():
    source=Path('notebooks/420_local_refine.py').read_text()
    fn=next(n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name=='attempt')
    block=next(n for n in fn.body if isinstance(n,ast.If) and "cfg['rfactor_limit']" in ast.unparse(n.test))
    # Exercise real driver mapping with deleted base ID 1 and nontrivial permutation.
    from types import SimpleNamespace
    env=dict(np=np,cfg={'rfactor_limit':.2},base_rmax=.2,
             base=SimpleNamespace(elements=np.array([[0,2,3]]),depths=np.ones(5)),
             nodes=np.zeros((4,2)),node_map=np.array([2,-1,0,3,1]),
             elements=np.array([[0,1,2]]),depths=np.ones(4),is_new=np.zeros(4,bool),out={},say=lambda *a:None)
    def fake(*a,**kw):
        return np.ones(4),dict(frozen_pair_edges_over_rmax=[[2,0],[0,1]],n_depths_changed=0,
                              max_depth_change_m=0,converged=False,rounds=1,max_frozen_depth_change_m=0)
    env['limit_rfactor']=fake
    stop=next(i for i,n in enumerate(block.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Name) and n.value.func.id=='say')
    exec(compile(ast.Module(body=block.body[:stop],type_ignores=[]),'<driver rfactor>','exec'),env)
    assert env['out']['rfactor']['n_new_frozen_pair_over_rmax']==1

@pytest.mark.parametrize('rmax',[0.,.2,.999])
def test_positive_limiter_bounds(rmax):
    out,rep=limit_rfactor(M().elements,np.array([8.,80.,.8,8.]),
                          np.array([False,True,True,False]),rmax,depth_min=.8,depth_max=80,rounds=500)
    assert rep['converged']
    assert np.array_equal(out[[0,3]],[8,8])
    assert np.isfinite(out).all() and out.min()>=.8 and out.max()<=80
    e=np.array([[0,1],[1,2],[0,2],[2,3],[0,3]])
    assert np.max(abs(out[e[:,0]]-out[e[:,1]])/out[e].sum(axis=1))<=rmax+1e-9

@pytest.mark.parametrize('angle',[0,.2,.8,1.7])
def test_gradation_affine_rotations(angle):
    f=lambda p: 100+np.asarray(p)@np.array([np.cos(angle),np.sin(angle)])
    rep=patch.field_gradation(f,shapely.box(0,0,100,100),spacing=10)
    assert rep['max_slope']==pytest.approx(1.)


def test_geometry_projected_crs_filter_buffer(tmp_path):
    import geopandas as gpd
    from pyproj import Transformer
    x,y=Transformer.from_crs(4326,32654,always_xy=True).transform(139.8,35.3)
    poly=shapely.box(x,y,x+100,y+100)
    path=tmp_path/'regions.gpkg'
    gpd.GeoDataFrame({'name':['a','b']},geometry=[poly,shapely.box(x+500,y,x+600,y+100)],crs=32654).to_file(path)
    geom,_=_geometry_from_file({'file':str(path),'where':{'name':'a'},'buffer_m':20})
    back=shapely.ops.transform(Transformer.from_crs(4326,32654,always_xy=True).transform,geom)
    assert abs(back.area-(10000+8000+np.pi*400))<150
    with pytest.raises(ValueError): _geometry_from_file({'file':str(path)})
    with pytest.raises(ValueError): _geometry_from_file({'file':str(path),'index':2})


def test_invalid_hole_before_buffer(tmp_path):
    p=shapely.Polygon([(139,35),(139.01,35),(139.01,35.01),(139,35.01)],
                      holes=[[(140,36),(140.01,36),(140.01,36.01),(140,36.01)]])
    path=tmp_path/'bad.geojson'; path.write_text(json.dumps(shapely.geometry.mapping(p)))
    with pytest.raises(ValueError): _geometry_from_file({'file':str(path),'buffer_m':100})

@pytest.mark.parametrize('value',['nan','inf','-inf'])
def test_dep_depth_nonfinite(tmp_path,value):
    files=export_fvcom_case(M(),tmp_path,'b',obc_depth_control=False)
    rows=files['dep'].read_text().splitlines(); fields=rows[1].split(); fields[2]=value; rows[1]=' '.join(fields)
    files['dep'].write_text('\n'.join(rows)+'\n')
    with pytest.raises(ValueError): read_fvcom_case(files['grd'],files['dep'],files['obc'])


def test_mixed_types_refused_and_default_preserved(tmp_path):
    files=export_fvcom_case(M(),tmp_path/'a','a',obc_type=3,obc_depth_control=False)
    m=read_fvcom_case(files['grd'],files['dep'],files['obc'])
    out=export_fvcom_case(m,tmp_path/'b','b')
    assert read_obc_types(out['obc'])==[3]
    rows=files['obc'].read_text().splitlines(); rows[-1]=rows[-1][:-1]+'2'
    files['obc'].write_text('\n'.join(rows)+'\n')
    with pytest.raises(ValueError): read_fvcom_case(files['grd'],files['dep'],files['obc'])


def test_ceiling_holes_zero_width_and_order():
    xy,tri=H['grid_mesh'](15,15)
    g=shapely.Polygon([(200,200),(1200,200),(1200,1200),(200,1200)],
                      holes=[[(500,500),(900,500),(900,900),(500,900)]])
    rs=[(g,20,300,0),(shapely.box(400,400,1000,1000),80,0,999)]
    pts=np.array([[300,300],[700,700],[600,700],[100,100]])
    h=patch.patch_sizing(xy,tri,rs,distmesh_scale=1)(pts)
    pieces=[patch.patch_sizing(xy,tri,[r],distmesh_scale=1)(pts) for r in rs]
    assert np.isfinite(h).all() and np.allclose(h,np.min(pieces,axis=0))
    assert np.allclose(h,patch.patch_sizing(xy,tri,rs[::-1],distmesh_scale=1)(pts))

@pytest.mark.parametrize('dte',[1.7,2.3,5.0])
def test_ramp_rounding(dte):
    import re
    env=runpy.run_path('notebooks/383_m2_case_prep.py'); fn=env['namelist']; fn.__globals__['DTE']=dte
    nml=fn(Path('/tmp/input'),Path('/tmp/output'))
    ramp=int(re.search(r'(?m)^\s*IRAMP\s*=\s*(\d+)',nml).group(1))
    assert abs(ramp*dte*env['ISPLIT']-env['RAMP_SECONDS'])<=dte*env['ISPLIT']/2


def test_driver_cap_rejects_large_unchanged_base():
    m=M(); m.nodes,m.elements=H['grid_mesh'](73,73)
    m.nodes[:,1]*=.2; m.depths=np.full(len(m.nodes),8.)
    m.open_boundaries=[]; m.land_boundaries=[]
    qa=run_qa(m,max_offenders=10000)
    c=next(c for c in qa.checks if c.check_id=='c1_min_angle')
    assert c.n_violations==10368 and len(c.offenders)==10000
    assert patch.introduced_violations([c],len(m.elements),m.elements)==[]
