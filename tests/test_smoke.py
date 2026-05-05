import fvcom_mesh_tools
from fvcom_mesh_tools.backend import Backend


def test_version_string() -> None:
    assert isinstance(fvcom_mesh_tools.__version__, str)
    assert fvcom_mesh_tools.__version__


def test_backend_protocol_runtime_check() -> None:
    class Dummy:
        name = "dummy"

        def is_available(self) -> bool:
            return False

    assert isinstance(Dummy(), Backend)
