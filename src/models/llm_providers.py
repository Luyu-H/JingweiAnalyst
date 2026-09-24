"""
Which language model the agent talks to, and how it is built.

The agent uses a chat model in one way only - `invoke([SystemMessage, HumanMessage])`
returning a message whose `.content` is text - so any LangChain chat model can drive
it. What differs between vendors is the client class, the endpoint, the environment
variable holding the key, and a handful of argument names. That knowledge lives in
the table below rather than in the agent, so adding a vendor is a table entry and
switching vendor is one line of configs/agent_with_tools.yaml.

Three kinds of client cover the field:

  openai_compatible  one endpoint speaking the OpenAI chat-completions protocol,
                     driven by ChatOpenAI with a base_url. DeepSeek, OpenAI itself,
                     Qwen's DashScope compatible mode, and every local or hosted
                     server that mimics the protocol (vLLM, Ollama, SiliconFlow,
                     Moonshot, ...) are all this kind, so they need no extra package.
  anthropic          Claude, through langchain-anthropic.
  google             Gemini, through langchain-google-genai.

The two native clients are imported only when they are actually asked for, so the
project keeps working with langchain-openai alone.
"""

import os
from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Tuple

# Client arguments accepted from the configuration for every provider. Anything
# else a vendor's client understands can still be passed through `extra_params`.
COMMON_PARAMS = ("max_tokens", "timeout", "max_retries", "top_p")

# Every key the `llm` section understands. Named here so a caller overriding one
# setting can tell a model setting from an agent setting without guessing.
LLM_KEYS = ("provider", "model_name", "temperature", "base_url", "api_key_env",
            "extra_params") + COMMON_PARAMS


@dataclass(frozen=True)
class ProviderSpec:
    """One vendor: how to reach it, and what to assume when the config is silent."""

    kind: str                                   # openai_compatible | anthropic | google
    default_model: str
    api_key_env: str
    base_url: Optional[str] = None              # None: the client's own default endpoint
    package: Optional[str] = None               # pip name, when a client must be installed
    import_name: Optional[str] = None
    aliases: Tuple[str, ...] = ()
    defaults: Dict[str, Any] = field(default_factory=dict)
    note: str = ""


PROVIDERS: Dict[str, ProviderSpec] = {
    "deepseek": ProviderSpec(
        kind="openai_compatible",
        default_model="deepseek-reasoner",
        api_key_env="DEEPSEEK_API_KEY",
        base_url="https://api.deepseek.com",
        note="deepseek-reasoner (thinking) or deepseek-chat; temperature is ignored by the reasoner.",
    ),
    "openai": ProviderSpec(
        kind="openai_compatible",
        default_model="gpt-5",
        api_key_env="OPENAI_API_KEY",
        aliases=("gpt",),
        note="The reasoning models accept only their default temperature; set temperature: null for those.",
    ),
    "anthropic": ProviderSpec(
        kind="anthropic",
        default_model="claude-sonnet-5",
        api_key_env="ANTHROPIC_API_KEY",
        package="langchain-anthropic",
        import_name="langchain_anthropic",
        aliases=("claude",),
        # The client's own default is 1024 tokens, which truncates generated code
        # mid-function; the agent then repairs a syntax error that was never the
        # model's mistake.
        defaults={"max_tokens": 8192},
    ),
    "google": ProviderSpec(
        kind="google",
        default_model="gemini-2.5-pro",
        api_key_env="GOOGLE_API_KEY",
        package="langchain-google-genai",
        import_name="langchain_google_genai",
        aliases=("gemini",),
    ),
    "qwen": ProviderSpec(
        kind="openai_compatible",
        default_model="qwen-plus",
        api_key_env="DASHSCOPE_API_KEY",
        # Alibaba's OpenAI-compatible endpoint. Mainland China accounts use
        # https://dashscope.aliyuncs.com/compatible-mode/v1 instead; set base_url.
        base_url="https://dashscope-intl.aliyuncs.com/compatible-mode/v1",
        aliases=("dashscope", "tongyi"),
        note="qwen-max / qwen-plus / qwen-turbo, or a qwen3-* model id.",
    ),
    "openai_compatible": ProviderSpec(
        kind="openai_compatible",
        default_model="",
        api_key_env="LLM_API_KEY",
        aliases=("custom", "compatible"),
        note="Anything else speaking the OpenAI protocol; base_url and model_name are required.",
    ),
}

_ALIASES = {alias: name for name, spec in PROVIDERS.items() for alias in spec.aliases}


def list_providers() -> Tuple[str, ...]:
    """Provider names the configuration accepts, aliases aside."""
    return tuple(PROVIDERS)


def resolve_provider(name: str) -> Tuple[str, ProviderSpec]:
    """Map a configured name, or one of its aliases, onto its spec."""
    key = (name or "").strip().lower().replace("-", "_")
    key = _ALIASES.get(key, key)
    if key not in PROVIDERS:
        raise ValueError(
            f"Unknown llm.provider {name!r}. Known providers: "
            f"{', '.join(sorted(PROVIDERS))} "
            f"(aliases: {', '.join(sorted(_ALIASES))}). "
            f"For a vendor not listed, use provider: openai_compatible with its base_url."
        )
    return key, PROVIDERS[key]


def resolve_settings(llm_cfg: Dict[str, Any]) -> Dict[str, Any]:
    """
    Turn the `llm` section into exactly what will reach the client.

    Whatever the configuration leaves out falls back to the provider's entry above,
    so a section naming only a provider and a model is enough. The result doubles as
    the record written into each run's manifest.
    """
    llm_cfg = dict(llm_cfg or {})
    name, spec = resolve_provider(llm_cfg.get("provider", "deepseek"))

    model = llm_cfg.get("model_name") or spec.default_model
    base_url = llm_cfg.get("base_url") or spec.base_url
    if not model:
        raise ValueError(f"llm.model_name is required for provider {name!r}.")
    if spec.kind == "openai_compatible" and name == "openai_compatible" and not base_url:
        raise ValueError("llm.base_url is required for provider 'openai_compatible'.")

    params = dict(spec.defaults)
    for key in COMMON_PARAMS:
        if llm_cfg.get(key) is not None:
            params[key] = llm_cfg[key]
    params.update(llm_cfg.get("extra_params") or {})

    return {
        "provider": name,
        "kind": spec.kind,
        "model": model,
        # None means "do not send it": the reasoning models reject any temperature
        # but their own, and a default of 0.0 would make them fail on every call.
        "temperature": llm_cfg.get("temperature"),
        "base_url": base_url,
        "api_key_env": llm_cfg.get("api_key_env") or spec.api_key_env,
        "params": params,
    }


def resolve_api_key(settings: Dict[str, Any], api_key: Optional[str] = None,
                    hint: Optional[str] = None) -> str:
    """
    The key comes from the caller or from the environment variable the settings name.

    It is deliberately never read out of the YAML file, which is committed.
    """
    if api_key:
        return api_key
    variable = settings["api_key_env"]
    key = os.environ.get(variable, "")
    if not key:
        where = f" in {hint}" if hint else ""
        raise SystemExit(
            f"{variable} is not set in the environment, and provider "
            f"'{settings['provider']}' needs it. Export it, or point llm.api_key_env"
            f"{where} at the variable that holds the key."
        )
    return key


def _import_client(spec: ProviderSpec, attribute: str):
    """Import a vendor's client late, and say what to install when it is missing."""
    import importlib

    try:
        module = importlib.import_module(spec.import_name)
    except ImportError as exc:
        raise ImportError(
            f"This provider needs the {spec.package} package: pip install {spec.package}"
        ) from exc
    return getattr(module, attribute)


def build_chat_model(settings: Dict[str, Any], api_key: Optional[str] = None,
                     hint: Optional[str] = None):
    """Construct the chat model described by `resolve_settings`."""
    key = resolve_api_key(settings, api_key, hint)
    params = dict(settings["params"])
    temperature = settings.get("temperature")
    if temperature is not None:
        params["temperature"] = temperature

    kind = settings["kind"]
    if kind == "openai_compatible":
        from langchain_openai import ChatOpenAI

        if settings.get("base_url"):
            params["base_url"] = settings["base_url"]
        return ChatOpenAI(model=settings["model"], api_key=key, **params)

    if kind == "anthropic":
        spec = PROVIDERS[settings["provider"]]
        chat_anthropic = _import_client(spec, "ChatAnthropic")
        if settings.get("base_url"):
            params["base_url"] = settings["base_url"]
        return chat_anthropic(model=settings["model"], api_key=key, **params)

    if kind == "google":
        spec = PROVIDERS[settings["provider"]]
        chat_google = _import_client(spec, "ChatGoogleGenerativeAI")
        # The Gemini client spells this one differently from everyone else.
        if "max_tokens" in params:
            params["max_output_tokens"] = params.pop("max_tokens")
        return chat_google(model=settings["model"], google_api_key=key, **params)

    raise ValueError(f"Unsupported provider kind: {kind!r}")


def build_llm(llm_cfg: Dict[str, Any], api_key: Optional[str] = None,
              hint: Optional[str] = None) -> Tuple[Any, Dict[str, Any]]:
    """The chat model for this `llm` section, and the settings it was built from."""
    settings = resolve_settings(llm_cfg)
    return build_chat_model(settings, api_key, hint), settings


if __name__ == "__main__":
    # python -m src.models.llm_providers - what can be configured, and whether the
    # key and client for each provider are actually present on this machine.
    import importlib.util

    header = (f"{'provider':<18} {'default model':<20} {'key variable':<20} "
              f"{'client package':<24} status")
    print(header)
    print("-" * len(header))
    for name, spec in PROVIDERS.items():
        installed = (spec.package is None
                     or importlib.util.find_spec(spec.import_name) is not None)
        if not installed:
            status = f"pip install {spec.package}"
        elif not os.environ.get(spec.api_key_env):
            status = f"export {spec.api_key_env}=..."
        else:
            status = "ready"
        print(f"{name:<18} {spec.default_model or '(model_name required)':<20} "
              f"{spec.api_key_env:<20} {(spec.package or 'langchain-openai'):<24} {status}")
        if spec.aliases:
            print(f"{'':<18} also written as: {', '.join(spec.aliases)}")
        if spec.note:
            print(f"{'':<18} {spec.note}")
