import ast
from pathlib import Path

import numpy as np
import pytest

from fvcom_mesh_tools.sizing import apply_sizing_regions, load_sizing


def lattice():
    x, y = np.meshgrid(np.arange(5)/111000, np.arange(5)/111000, indexing="ij")
    return np.full((5, 5), 100.), x, y


def region(target=30, **kw):
    return dict(name="port", geometry={"bbox": [0, 0, 1.1/111000, 4.1/111000]},
                target_h_m=target, **kw)


def apply(regions, **kwargs):
    v, x, y = lattice()
    return apply_sizing_regions(v, x, y, regions, gradation=kwargs.pop("gradation", 1000),
                                **kwargs)


def test_apply_pure_and_hmin():
    v, x, y = lattice()
    out, report = apply_sizing_regions(v, x, y, [region()], gradation=1000, hmin_m=40)
    assert np.all(v == 100)
    assert np.all(out[:2] == 40)
    assert np.all(out[2:] == 100)
    r = report["regions"][0]
    assert r["hmin_yield_fraction"] == 1
    assert r["area_m2"] == pytest.approx(4.51)
    assert r["achieved_median_m"] == 40


def test_overlap():
    a, b = region(), dict(region(60), name="other")
    first, _ = apply([a, b])
    second, _ = apply([b, a])
    np.testing.assert_array_equal(first, second)
    assert first[0, 0] == 30
    b["priority"] = 1
    assert apply([a, b])[0][0, 0] == 60


def test_transition():
    out, _ = apply([region(transition_m=2)])
    assert out[2, 0] == pytest.approx(30+70*.9/2)
    assert out[3, 0] == pytest.approx(30+70*1.9/2)
    assert out[4, 0] == 100


def test_gradation():
    out, _ = apply([region()], gradation=2)
    assert out[4, 0] == pytest.approx(36)
    assert np.max(abs(np.diff(out, axis=0))) <= 2+1e-10
    assert out.min() == 30


def test_cfl_diagnostic_only():
    out, report = apply([region()], cfl={"dt_s": 15, "cr": .45,
                                        "depth_m": np.full((5, 5), 10.)})
    r = report["regions"][0]
    assert out[0, 0] == 30
    assert r["target_below_cfl_fraction"] == 1
    assert r["cfl_floor_max_m"] == pytest.approx(np.sqrt(98.1)*15/.45)
    assert report["implied_dt_s"] == pytest.approx(.45*30/np.sqrt(98.1))
    assert r["implied_dt_cr1_s"] == pytest.approx(30/np.sqrt(98.1))


def test_empty_dry_unsampled():
    r = dict(region(), geometry={"bbox": [1, 1, 2, 2]})
    _, report = apply([r], cfl={"dt_s": 15, "cr": .45, "depth_m": np.zeros((5, 5))})
    assert report["implied_dt_s"] is None
    assert report["regions"][0]["achieved_min_m"] is None


def test_polygon_hole():
    r = region()
    r["geometry"] = {"type": "Polygon", "coordinates": [
        [[-1,-1],[5,-1],[5,5],[-1,5],[-1,-1]],
        [[1,1],[3,1],[3,3],[1,3],[1,1]]]}
    r["geometry"]["coordinates"] = (np.array(r["geometry"]["coordinates"])/111000).tolist()
    out, _ = apply([r])
    assert out[0, 0] == 30
    assert out[2, 2] == 100


def test_loader(tmp_path):
    cfg = load_sizing(Path(__file__).resolve().parents[1]/"recipes/sizing/tokyo_bay.yaml")
    assert cfg["coastal_target_m"] == 290
    assert cfg["regions"] == []
    import yaml
    for mutate in (lambda c: c.update(typo=1), lambda c: c.update(gradation=-1),
                   lambda c: c.update(regions=[region(target=-30)]),
                   lambda c: c.update(cfl={"dt_s": 15}),
                   lambda c: c.update(regions=[region(), region()]),
                   lambda c: c.update(regions=[dict(region(), geometry={"bbox": [2,1,0,3]})])):
        bad = dict(cfg)
        mutate(bad)
        p = tmp_path/"bad.yaml"
        p.write_text(yaml.safe_dump(bad))
        with pytest.raises(ValueError):
            load_sizing(p)


@pytest.mark.parametrize("change", ["depth", "shape", "nan"])
def test_invalid_arrays(change):
    v, x, y = lattice()
    kwargs = {}
    if change == "depth":
        kwargs["cfl"] = {"dt_s": 15, "cr": .45, "depth_m": np.ones(2)}
    elif change == "shape":
        x = x[:-1]
    else:
        v[0, 0] = np.nan
    with pytest.raises(ValueError):
        apply_sizing_regions(v, x, y, [], gradation=.165, **kwargs)


def test_default_hook_dry_run(monkeypatch):
    """Execute actual notebook opt-in hooks with SR_SIZING absent, without inputs."""
    import os
    monkeypatch.delenv("SR_SIZING", raising=False)
    source = (Path(__file__).resolve().parents[1]/"notebooks/325_sample_repro.py").read_text()
    tree = ast.parse(source)
    hooks = [n for n in tree.body if isinstance(n, ast.If) and
             ("SR_SIZING" in ast.unparse(n.test) or "_sizing_recipe" in ast.unparse(n.test))]
    assert len(hooks) == 2
    class Grid:
        values = np.arange(25., dtype=float).reshape(5, 5)
    g = Grid()
    before = g.values.tobytes()
    env = dict(os=os, _sizing_recipe=None, g=g)
    exec(compile(ast.Module(body=hooks, type_ignores=[]), "hooks", "exec"), env)
    assert g.values.tobytes() == before


def test_recipe_empty_hooks_preserve_values(monkeypatch):
    import os
    recipe = Path(__file__).resolve().parents[1]/"recipes/sizing/tokyo_bay.yaml"
    monkeypatch.setenv("SR_SIZING", str(recipe))
    source = (recipe.parents[2]/"notebooks/325_sample_repro.py").read_text()
    tree = ast.parse(source)
    hooks = [n for n in tree.body if isinstance(n, ast.If) and
             ("SR_SIZING" in ast.unparse(n.test) or "_sizing_recipe" in ast.unparse(n.test))]
    class Grid:
        values = np.arange(25., dtype=float).reshape(5, 5)
    g = Grid()
    before = g.values.tobytes()
    env = dict(os=os, _sizing_recipe=None, g=g)
    exec(compile(ast.Module(body=hooks, type_ignores=[]), "hooks", "exec"), env)
    assert g.values.tobytes() == before
    assert [env[k] for k in ("H0", "MAXEL", "GRADE", "DT", "CRMIN")] == [
        290, 1400, .165, 15, .45]


def test_measurement_script(tmp_path):
    import subprocess
    import sys
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        [sys.executable, str(root/"notebooks/391_sizing_report.py"),
         str(root/"tests/fixtures/tiny.fort14"), str(root/"recipes/sizing/tokyo_bay.yaml")],
        check=True, capture_output=True, text=True)
    assert "not historical provenance" in result.stderr
    assert "measured_dt_cr1_s" in result.stdout
    assert len(result.stdout.splitlines()) == 3


@pytest.mark.parametrize("content", ["regions: []\nregions: []", "- bad", "regions: ["])
def test_malformed_yaml(tmp_path, content):
    import yaml
    path = tmp_path/"malformed.yaml"
    path.write_text(content)
    with pytest.raises((ValueError, yaml.YAMLError)):
        load_sizing(path)
