"""Extension recipes (recipes/extend/*.yaml)."""

from pathlib import Path

import pytest
import yaml

from fvcom_mesh_tools.extend_recipe import REQUIRED_SETTINGS, load_extend_recipe

REPO = Path(__file__).resolve().parents[1]
RECIPE = REPO / "recipes/extend/tokyo_bay_enshu.yaml"


def _write(tmp_path, **change):
    raw = yaml.safe_load(RECIPE.read_text())
    base = tmp_path / "base"
    base.mkdir(exist_ok=True)
    for kind in ("grd", "dep", "obc"):
        (base / f"{raw['base_case']}_{kind}.dat").write_text("x\n")
    (tmp_path / "obc.csv").write_text("lon,lat\n139.0,34.0\n139.5,34.0\n")
    raw.update(base="base", open_boundary="obc.csv")
    for k, v in change.items():
        if v is None:
            raw.pop(k, None)
        else:
            raw[k] = v
    p = tmp_path / "r.yaml"
    p.write_text(yaml.safe_dump(raw))
    return p


def test_the_enshu_recipe_writes_out_every_setting():
    raw = yaml.safe_load(RECIPE.read_text())
    assert set(raw["settings"]) == set(REQUIRED_SETTINGS)


def test_a_valid_recipe_loads_with_absolute_paths(tmp_path):
    r = load_extend_recipe(_write(tmp_path))
    assert Path(r["base"]).is_absolute() and Path(r["open_boundary"]).is_file()


@pytest.mark.parametrize("change, match", [
    ({"case": None}, "missing"),
    ({"extra": 1}, "unknown"),
    ({"bathymetry": {"sizing": ["nope"], "depths": ["srtm15plus"]}}, "unknown or repeated"),
    ({"bathymetry": {"sizing": ["m7001", "m7001"], "depths": ["m7001"]}}, "unknown or repeated"),
    ({"bathymetry": {"sizing": ["m7001"]}}, "exactly 'sizing' and 'depths'"),
    ({"depths": {"min_m": 3, "max_m": 1, "rfactor": 0.2}}, "min_m < max_m"),
    ({"depths": {"min_m": 3, "max_m": None, "rfactor": 1.5}}, "rfactor"),
    ({"land": {"bbox": [1, 2, 3]}}, "land"),
])
def test_a_bad_recipe_is_refused(tmp_path, change, match):
    with pytest.raises(ValueError, match=match):
        load_extend_recipe(_write(tmp_path, **change))


def test_settings_must_be_complete_and_positive(tmp_path):
    s = dict(yaml.safe_load(RECIPE.read_text())["settings"])
    s.pop("dm_scale")
    with pytest.raises(ValueError, match="dm_scale"):
        load_extend_recipe(_write(tmp_path, settings=s))
    s = dict(yaml.safe_load(RECIPE.read_text())["settings"], max_edge_m=-1)
    with pytest.raises(ValueError, match="positive"):
        load_extend_recipe(_write(tmp_path, settings=s))
    s = dict(yaml.safe_load(RECIPE.read_text())["settings"], coast_h_m=9000)
    with pytest.raises(ValueError, match="exceeds"):
        load_extend_recipe(_write(tmp_path, settings=s))


def test_a_missing_base_file_is_refused(tmp_path):
    p = _write(tmp_path)
    next((tmp_path / "base").glob("*_dep.dat")).unlink()
    with pytest.raises(ValueError, match="no .*_dep.dat"):
        load_extend_recipe(p)


@pytest.mark.parametrize("key, value, match", [
    ("cfl_dt_s", float("nan"), "finite"),
    ("max_iter", 0.5, "integer"),
    ("gen_seed", 0.5, "integer"),
    ("fin_seed", -1, r"\[0, 2\*\*32\)"),
])
def test_settings_must_be_finite_and_integral_where_counted(tmp_path, key, value, match):
    """Review F18: NaN steps and fractional seeds/iterations were accepted."""
    s = dict(yaml.safe_load(RECIPE.read_text())["settings"], **{key: value})
    with pytest.raises(ValueError, match=match):
        load_extend_recipe(_write(tmp_path, settings=s))


@pytest.mark.parametrize("case", ["../escape", "/abs/path", "a/b", "", "."])
def test_case_names_are_one_file_name_component(tmp_path, case):
    """Review F28: a path in the case name could leave the output directory."""
    with pytest.raises(ValueError, match="file-name component"):
        load_extend_recipe(_write(tmp_path, case=case))


@pytest.mark.parametrize("bbox", [[float("nan"), 33, 141, 36], [141, 36, 137, 33],
                                  [137, -95, 141, 36]])
def test_land_bbox_must_be_finite_and_ordered(tmp_path, bbox):
    """Review round 2 F11."""
    with pytest.raises(ValueError, match="land.bbox"):
        load_extend_recipe(_write(tmp_path, land={"bbox": bbox}))


def test_a_boundary_whose_report_hash_differs_is_refused(tmp_path):
    """Review round 3 F5: a half-published boundary (new report, old CSV)."""
    import json

    p = _write(tmp_path)
    obc = Path(load_extend_recipe(p)["open_boundary"])
    obc.with_suffix(".json").write_text(json.dumps({"csv_sha256": "0" * 64}))
    with pytest.raises(ValueError, match="hash mismatch"):
        load_extend_recipe(p)


def test_the_recipe_digest_is_of_the_bytes_parsed(tmp_path, monkeypatch):
    """Review round 5 F3: the driver and its stages must read one recipe."""
    import hashlib

    from fvcom_mesh_tools.extend_recipe import EXPECT_ENV, check_expected

    p = _write(tmp_path)
    r = load_extend_recipe(p)
    assert r["recipe_sha256"] == hashlib.sha256(p.read_bytes()).hexdigest()
    monkeypatch.setenv(EXPECT_ENV["recipe_sha256"], r["recipe_sha256"])
    monkeypatch.setenv(EXPECT_ENV["open_boundary_sha256"], r["open_boundary_sha256"])
    check_expected(r)
    p.write_text(p.read_text() + "# edited\n")
    with pytest.raises(ValueError, match="recipe changed"):
        check_expected(load_extend_recipe(p))


def test_a_boundary_left_for_recovery_is_refused(tmp_path):
    """Review round 5 F5: a failed rollback leaves a marker."""
    p = _write(tmp_path)
    (tmp_path / "obc.csv.RECOVER").write_text("{}")
    with pytest.raises(ValueError, match="RECOVER"):
        load_extend_recipe(p)


def test_the_recipe_name_is_one_path_component(tmp_path):
    """Review round 15 F7."""
    with pytest.raises(ValueError):
        load_extend_recipe(_write(tmp_path, name="x/../../escaped"))


def test_a_repeated_key_is_refused(tmp_path):
    """Review round 21 F7: the last of two values won silently."""
    p = _write(tmp_path)
    p.write_text(p.read_text() + "case: TokyoBayEnshu\n")
    with pytest.raises(ValueError, match="duplicate key 'case'"):
        load_extend_recipe(p)


def test_merge_keys_still_work():
    """Review round 22 F5."""
    from fvcom_mesh_tools.yaml_strict import load_unique

    assert load_unique("d:\n  <<: {a: 1, b: 2}\n  b: 3\n") == {"d": {"a": 1, "b": 3}}


def test_duplicates_inside_merged_mappings_are_refused():
    """Review round 23 F1."""
    from fvcom_mesh_tools.yaml_strict import load_unique

    for text in ("d:\n  <<: {a: 1, a: 2}\n", "x: &x {p: 1, p: 2}\ny:\n  <<: *x\n"):
        with pytest.raises(ValueError, match="duplicate key"):
            load_unique(text)
