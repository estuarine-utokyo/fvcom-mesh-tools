"""Policy parsing, generation wiring, and synthetic throat regression."""
import ast
from pathlib import Path

import numpy as np
import pytest

from fvcom_mesh_tools.channel_policy import resolve_narrow_channels
from fvcom_mesh_tools.one_wide import (
    configured_one_wide,
    generation_options,
    parse_one_wide,
)
from fvcom_mesh_tools.sizing import load_sizing
from tests.test_channel_policy_strict import _strip_with_one_wide_tail


@pytest.mark.parametrize('value', ['', 'ALLOW', None, True, 'allow_if_harmless'])
def test_invalid(value):
    with pytest.raises(ValueError, match='one_wide'):
        parse_one_wide(value)


def test_sources_and_knobs():
    assert configured_one_wide(environ={}) == 'forbid'
    assert configured_one_wide({'one_wide': 'allow'}, environ={}) == 'allow'
    assert configured_one_wide({'one_wide': 'allow'},
                               environ={'SR_ONE_WIDE': 'forbid'}) == 'forbid'
    assert generation_options('forbid', environ={}) == dict(
        feature_rows=3., min_rows=2, widen_factor=.875,
        attain_bar_h=1.5, force_two_rows=False)
    assert generation_options('allow', environ={'SR_FS': '9', 'SR_FORCE2ROWS': 'on',
                              'SR_ATTAIN_BAR': '9', 'SR_WIDEN_FACTOR': '9'}) == dict(
        feature_rows=1., min_rows=1, widen_factor=1., attain_bar_h=0., force_two_rows=False)


def test_recipe(tmp_path):
    source = Path('recipes/sizing/tokyo_bay.yaml').read_text()
    p = tmp_path / 'recipe.yaml'
    for mode in ('allow', 'forbid', 'invalid'):
        p.write_text(source.replace('one_wide: forbid', f'one_wide: {mode}'))
        if mode == 'invalid':
            with pytest.raises(ValueError, match='one_wide'):
                load_sizing(p)
        else:
            assert load_sizing(p)['one_wide'] == mode


def test_throat_kept_or_pruned_and_default_regression():
    mesh, nmain = _strip_with_one_wide_tail()
    kw = dict(min_basin_elements=25, apply_widen=False,
              strict_boundary_flag=True, max_rounds=8)
    default, di = resolve_narrow_channels(mesh, **kw)
    forbid, fi = resolve_narrow_channels(mesh, one_wide='forbid', **kw)
    allow, ai = resolve_narrow_channels(mesh, one_wide='allow', **kw)
    np.testing.assert_array_equal(default.elements, forbid.elements)
    np.testing.assert_array_equal(default.nodes, forbid.nodes)
    assert di['n_deleted_elements'] == fi['n_deleted_elements'] == 3
    assert len(forbid.elements) == nmain
    assert allow is mesh
    assert ai['n_deleted_elements'] == ai['n_widened'] == 0
    assert ai['n_flagged'] >= 2
    assert all(c['action'] == 'keep' for c in ai['clusters'])


def test_scripts_thread_every_call():
    for script in ('325_sample_repro.py', '331_finish2.py'):
        tree = ast.parse(Path('notebooks', script).read_text())
        names = {'detect_waterways', 'apply_waterway_policy', 'normalize_unresolved_water',
                 'resolve_narrow_channels', 'finish_obc_mesh'}
        calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
                 and isinstance(n.func, ast.Name) and n.func.id in names]
        assert calls
        assert all(any(k.arg == 'one_wide' for k in c.keywords) for c in calls)


def test_allow_geometry_and_normalization():
    from shapely.geometry import box
    from shapely.ops import unary_union

    from fvcom_mesh_tools.waterways import (
        apply_waterway_policy,
        detect_waterways,
        normalize_unresolved_water,
    )
    domain = box(0, 0, 20000, 10000)
    land = unary_union([box(8000, 2000, 12000, 4850),
                        box(8000, 5150, 12000, 8000)])
    kw = dict(h_mesh_m=350., metric_scale=(1., 1.))
    recs = detect_waterways(land, domain, obc_point=(1000, 5000),
                           one_wide='allow', min_resolve_width_frac=10, **kw)
    assert any(r['action'] == 'keep' for r in recs)
    allowed, ai = apply_waterway_policy(land, domain, recs, one_wide='allow',
                                      force_two_rows=True, attain_bar_h=10, **kw)
    assert ai['kept'] and not ai['band_n']
    recs = detect_waterways(land, domain, obc_point=(1000, 5000), **kw)
    default, di = apply_waterway_policy(land, domain, recs, **kw)
    recs = detect_waterways(land, domain, obc_point=(1000, 5000), one_wide='forbid', **kw)
    forbid, fi = apply_waterway_policy(land, domain, recs, one_wide='forbid', **kw)
    assert default.equals_exact(forbid, 0)
    assert allowed.area > forbid.area  # one-row target removes less land
    fills, info = normalize_unresolved_water(
        allowed, domain, obc_point=(1000, 5000), one_wide='allow', **kw)
    assert not fills and info['area_filled_ha'] == 0


def test_finish_reuses_generation_policy(tmp_path):
    from fvcom_mesh_tools.one_wide import finishing_one_wide
    p = tmp_path / 'channel_policy.json'
    assert finishing_one_wide(p, environ={}) == 'forbid'
    p.write_text('{"one_wide": "allow"}')
    assert finishing_one_wide(p, environ={}) == 'allow'
    with pytest.raises(ValueError, match='differs'):
        finishing_one_wide(p, environ={'SR_ONE_WIDE': 'forbid'})


@pytest.mark.parametrize('mode,expected', [('allow', []), ('forbid', ['split', 'widen'])])
def test_finish_width_operators(monkeypatch, mode, expected):
    from fvcom_mesh_tools import mesh_clean
    from fvcom_mesh_tools import mesh_clean_phase_h as h
    from fvcom_mesh_tools.algorithms import obc_finish as f
    from fvcom_mesh_tools.algorithms import perp_local
    mesh, _ = _strip_with_one_wide_tail()
    calls = []
    def noop(m, *a, **kw):
        return m, {}
    monkeypatch.setattr(perp_local, 'align_open_boundary_local', noop)
    monkeypatch.setattr(h, 'phase_h_finish', noop)
    monkeypatch.setattr(h, '_stochastic_local_fix_round', lambda *a, **kw: {})
    monkeypatch.setattr(mesh_clean, 'compact_nodes', noop)
    for name in ('flip_for_obc_perp', 'flip_c4_edges'):
        monkeypatch.setattr(f, name, lambda *a, **kw: {})
    monkeypatch.setattr(f, 'fix_r4', lambda *a, **kw: {'unfixed': []})
    for name in ('split_c4_edges', 'collapse_short_boundary_edges', 'cfl_polish'):
        monkeypatch.setattr(f, name, noop)

    def split(m):
        calls.append('split')
        return m, {}

    def widen(m, land):
        calls.append('widen')
        return m, {'widened': 0}

    monkeypatch.setattr(f, 'split_choke_edges', split)
    monkeypatch.setattr(f, 'widen_choke_sections', widen)
    result, _ = f.finish_obc_mesh(mesh, one_wide=mode, land_union=object())
    assert result is mesh and calls == expected
