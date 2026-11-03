"""The MCP tool must not override a capability's registered invoke_path.

`FactoryTools.invoke_capability` resolves `path or cap.invoke_path or "/"`. The MCP tool
declared `path: str = "/"`, and "/" is truthy, so every call made through an assistant went to
the service root instead of the registered route — `chemprop-admet` was POSTed to `:8572/`
rather than `/admet`. The library default (None) routed correctly, which is why nothing in the
composer path showed the bug, and no test covered `build_server`.
"""
import inspect

from hcls_common import mcp_server


def _invoke_tool_signature():
    """Pull invoke_capability's signature out of build_server without needing `mcp`.

    build_server imports FastMCP lazily and the package is an optional extra, so the source
    is read rather than executed.
    """
    src = inspect.getsource(mcp_server.build_server)
    marker = "def invoke_capability("
    start = src.index(marker)
    end = src.index(")", start)
    return src[start:end]


def test_mcp_path_parameter_defaults_to_none_not_slash():
    sig = _invoke_tool_signature()
    assert "path: str | None = None" in sig, (
        'invoke_capability must default path to None; a default of "/" is truthy and '
        "overrides each capability's invoke_path"
    )
    assert 'path: str = "/"' not in sig


def test_factory_tools_resolves_registered_path_when_none_is_passed():
    """None must fall through to the capability's own invoke_path."""
    src = inspect.getsource(mcp_server.FactoryTools.invoke_capability)
    assert 'path or getattr(cap, "invoke_path", None) or "/"' in src


def test_build_server_names_the_missing_package():
    """A documented `python -m hcls_common.mcp_server` failed with a bare import error."""
    src = inspect.getsource(mcp_server.build_server)
    assert "hcls-common[mcp]" in src
