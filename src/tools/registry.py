import importlib
import pkgutil
import os
from typing import Callable, Dict, Any

TOOLS: Dict[str, Dict[str, Any]] = {}

def register_tools(name: str, description: str, args: Dict[str, str], returns: str):
    """
    Register a utility function in TOOLS.

    The first line of `description` is the tool's summary: the planning step is
    shown only that line, and the code-generation step gets the whole entry for
    the tools the plan selected. Write it as one sentence saying what the tool
    is for.

    Args:
        name (str): Tool name, which must match the decorated function's name.
        description (str): Summary line, then any detail below it.
        args (dict): Parameter description, format {parameter name: description}
        returns (str): Return value description
    """
    def decorator(func: Callable):
        tool_name = name or func.__name__
        if tool_name != func.__name__:
            raise RuntimeError(
                f"register_tools(name='{tool_name}') is decorating '{func.__name__}'. The decorator binds "
                f"to the function defined immediately below it; move any helper above the decorator, or "
                f"correct the registered name."
            )
        TOOLS[tool_name] = {
            "func": func,
            "description": description or func.__doc__ or "",
            "summary": (description or func.__doc__ or "").strip().split("\n")[0],
            "args": args or {},
            "returns": returns or ""
        }
        return func
    return decorator


def auto_import_tools(package: str = "src.tools"):
    """
    Automatically import all modules under the tools directory.

    Args:
        package: Dotted package name holding the tool modules.
    """
    # Resolve against the repository root rather than the process working directory, so importing the agent from anywhere still finds the tools.
    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    package_path = os.path.join(repo_root, *package.split("."))
    if not os.path.isdir(package_path):
        raise RuntimeError(f"Package path not found: {package_path}")

    for _, module_name, is_pkg in pkgutil.iter_modules([package_path]):
        if is_pkg:
            continue
        if module_name == "registry":
            continue
        importlib.import_module(f"{package}.{module_name}")


def tool_summaries() -> Dict[str, str]:
    """
    One line per tool: what the planning step is shown.

    Sending every full description costs thousands of tokens per call and asks
    the planner to read detail it cannot act on yet. It only has to rule the
    irrelevant tools out; the code-generation step receives the full entry for
    whatever the plan selected.
    """
    return {name: meta.get("summary") or meta["description"].strip().split("\n")[0]
            for name, meta in TOOLS.items()}


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
