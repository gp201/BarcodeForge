import toml
from barcodeforge import __version__ as init_version


def test_version_consistency():
    """Tests if the version in pyproject.toml and barcodeforge/__init__.py are consistent."""
    try:
        with open("pyproject.toml", "r") as f:
            data = toml.load(f)
        project_version = data["project"]["version"]
    except FileNotFoundError:
        assert False, "pyproject.toml not found"
    except KeyError:
        assert False, "Version not found in pyproject.toml"

    assert (
        project_version == init_version
    ), f"Version mismatch: pyproject.toml has {project_version}, __init__.py has {init_version}"
