"""Base-mesh recipes (recipes/base/*.yaml)."""

from pathlib import Path

import numpy as np
import pytest
import yaml

from fvcom_mesh_tools.base_recipe import (
    REQUIRED_SETTINGS,
    compare_to_reference,
    load_base_recipe,
    read_open_boundary,
)

REPO = Path(__file__).resolve().parents[1]
RECIPE = REPO / "recipes/base/tokyo_bay_tool.yaml"

# the arc notebook 325 carried in its code before it became an input
GOTO2023_ARC = np.array([
    [139.6713, 35.1396], [139.6737, 35.1288], [139.6772, 35.1168],
    [139.6816, 35.1031], [139.6871, 35.0877], [139.6946, 35.0705],
    [139.7000, 35.0576], [139.7069, 35.0445], [139.7134, 35.0327],
    [139.7216, 35.0184], [139.7289, 35.0047], [139.7373, 34.9916],
    [139.7497, 34.9750]])


def test_the_tokyo_bay_recipe_loads_and_writes_out_every_setting():
    r = load_base_recipe(RECIPE)
    assert set(REQUIRED_SETTINGS) <= set(r["settings"])
    assert Path(r["open_boundary"]).is_file() and Path(r["edits"]).is_dir()
    assert r["reference"]["n_nodes"] == 4734


def test_the_open_boundary_input_is_the_arc_the_code_carried_bit_for_bit():
    arc = np.array(read_open_boundary(REPO / "recipes/base/tokyo_bay_obc.csv"))
    assert np.array_equal(arc, GOTO2023_ARC)


def _write(tmp_path, **change):
    raw = yaml.safe_load(RECIPE.read_text())
    for k, v in change.items():
        if v is None:
            raw.pop(k, None)
        else:
            raw[k] = v
    for f in ("tokyo_bay_obc.csv", "tokyo_bay_domain.json"):
        (tmp_path / f).write_text((REPO / "recipes/base" / f).read_text())
    raw["edits"] = str(REPO / "recipes/edits/sample_repro")
    p = tmp_path / "r.yaml"
    p.write_text(yaml.safe_dump(raw))
    return p


def test_a_missing_setting_is_refused(tmp_path):
    raw = yaml.safe_load(RECIPE.read_text())
    s = dict(raw["settings"])
    s.pop("SR_FIN_SEED")
    with pytest.raises(ValueError, match="SR_FIN_SEED"):
        load_base_recipe(_write(tmp_path, settings=s))


def test_an_unknown_or_path_setting_is_refused(tmp_path):
    raw = yaml.safe_load(RECIPE.read_text())
    for bad in ("SR_H00", "SR_OUT"):
        s = dict(raw["settings"], **{bad: "1"})
        with pytest.raises(ValueError, match="unknown settings"):
            load_base_recipe(_write(tmp_path, settings=s))


def test_a_missing_key_or_bad_depths_is_refused(tmp_path):
    with pytest.raises(ValueError, match="missing"):
        load_base_recipe(_write(tmp_path, open_boundary=None))
    with pytest.raises(ValueError, match="depths"):
        load_base_recipe(_write(tmp_path, depths="srtm"))


@pytest.mark.parametrize("text", ["lon,lat\n139.7,35.1\n", "lon,lat\n139.7,abc\n1,2\n",
                                  "lon,lat\n139.7,95\n139.8,35\n"])
def test_a_bad_open_boundary_is_refused(tmp_path, text):
    p = tmp_path / "obc.csv"
    p.write_text(text)
    with pytest.raises(ValueError):
        read_open_boundary(p)


def test_a_malformed_reference_is_refused(tmp_path):
    ref = yaml.safe_load(RECIPE.read_text())["reference"]
    for bad in ({k: v for k, v in ref.items() if k != "wet_area_km2"},
                dict(ref, fort14_sha256="abc"),
                dict(ref, tolerance={"n_nodes": 0.05}),
                dict(ref, tolerance=dict(ref["tolerance"], n_nodes=1.5)),
                dict(ref, n_nodes=0)):
        with pytest.raises(ValueError, match="reference"):
            load_base_recipe(_write(tmp_path, reference=bad))


def _summary(ref, **kw):
    s = {k: ref[k] for k in ("fort14_sha256", "n_nodes", "n_elements", "n_obc_nodes",
                             "wet_area_km2")}
    s.update(kw)
    return s


def test_the_reference_is_reproduced_without_byte_identity():
    """A slightly different mesh reproduces; byte identity is reported apart."""
    ref = load_base_recipe(RECIPE)["reference"]
    same = compare_to_reference(_summary(ref), ref, qa_passed=True)
    assert same["byte_identical"] and same["reproduces"]
    near = compare_to_reference(
        _summary(ref, fort14_sha256="0" * 64, n_nodes=ref["n_nodes"] + 100,
                 wet_area_km2=ref["wet_area_km2"] * 1.005), ref, qa_passed=True)
    assert near["reproduces"] and not near["byte_identical"]


@pytest.mark.parametrize("kw, qa", [({}, False),
                                    ({"n_obc_nodes": 14}, True),
                                    ({"n_elements": 9000}, True),
                                    ({"wet_area_km2": 1400.0}, True)])
def test_the_reference_is_not_reproduced(kw, qa):
    ref = load_base_recipe(RECIPE)["reference"]
    assert not compare_to_reference(_summary(ref, **kw), ref, qa_passed=qa)["reproduces"]
