"""fmesh-refine: a thin, careful wrapper round the GPL-side generator."""

from __future__ import annotations

from fvcom_mesh_tools.cli.plot_views import main as plot_main
from fvcom_mesh_tools.cli.refine_run import main, repo_root


def test_the_generator_is_found_next_to_the_package():
    assert (repo_root() / "notebooks" / "420_local_refine.py").exists()


def test_dry_run_prints_the_command_and_runs_nothing(tmp_path, capsys):
    recipe = tmp_path / "r.yaml"
    recipe.write_text("refine: []\n")
    land = tmp_path / "land.shp"
    land.write_text("")
    assert main([str(recipe), "--out", str(tmp_path / "o"), "--land", str(land),
                 "--seeds", "0:1", "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "420_local_refine.py" in out and "LR_SEEDS=0,1" in out
    assert not (tmp_path / "o").exists()


def test_it_refuses_a_missing_recipe_missing_land_and_an_existing_result(tmp_path):
    land = tmp_path / "land.shp"
    land.write_text("")
    assert main([str(tmp_path / "none.yaml"), "--land", str(land)]) == 2
    recipe = tmp_path / "r.yaml"
    recipe.write_text("refine: []\n")
    assert main([str(recipe), "--land", str(tmp_path / "no.shp")]) == 2
    done = tmp_path / "done"
    done.mkdir()
    (done / "report.json").write_text("{}")
    assert main([str(recipe), "--out", str(done), "--land", str(land)]) == 2


def test_plot_views_refuses_a_directory_without_a_mesh(tmp_path):
    assert plot_main([str(tmp_path)]) == 2
