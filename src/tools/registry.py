import importlib
import pkgutil
import os
from typing import Callable, Dict, Any

# Global tool registry
TOOLS: Dict[str, Dict[str, Any]] = {}

def register_tools(name: str, description: str, args: Dict[str, str], returns: str):
    """
    Register utility functions to TOOLS.

    Args:
        name (str): Tool name
        description (str): Tool description
        args (dict): Parameter description, format {Parameter name: type string}
        returns (str): Return value description
    """
    def decorator(func: Callable):
        tool_name = name or func.__name__
        TOOLS[tool_name] = {
            "func": func,
            "description": description or func.__doc__ or "",
            "args": args or {},
            "returns": returns or ""
        }
        return func
    return decorator


def auto_import_tools(package: str = "src.tools"):
    """
    Automatically import all modules under the tools directory.

    Args:
        path: The path to the tools path.
    """
    package_path = package.replace(".", "/")
    if not os.path.isdir(package_path):
        raise RuntimeError(f"Package path not found: {package_path}")

    for _, module_name, is_pkg in pkgutil.iter_modules([package_path]):
        if is_pkg:
            continue
        if module_name == "registry":
            continue
        importlib.import_module(f"{package}.{module_name}")


def export_tools_doc(format: str = "markdown"):
    """
    Export information of all the tools, including name, description, and args.

    Args:
        format (str): "markdown" or "json"
    """
    if format == "markdown":
        lines = []
        for name, meta in TOOLS.items():
            lines.append(f"## {name}")
            lines.append(f"- **Description**: {meta['description']}")
            if meta["args"]:
                lines.append(f"- **Arguments**:")
                for arg, typ in meta["args"].items():
                    lines.append(f"  - `{arg}` ({typ})")
            if meta["returns"]:
                lines.append(f"- **Returns**: {meta['returns']}")
            lines.append("")
        return "\n".join(lines)

    elif format == "json":
        import json
        return json.dumps(
            {k: {kk: vv for kk, vv in v.items() if kk != "func"} for k, v in TOOLS.items()},
            indent=2
        )

    else:
        raise ValueError("Unsupported format. Use 'markdown' or 'json'.")
