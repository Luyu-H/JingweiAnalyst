"""
Retrieval of knowledge fragments (docs/knowledge) for the agent's planning step.

The agent runs this before it plans. `KnowledgeRetriever.retrieve` goes through the
stages: R1 prefilter, R2a selection, R2b condensing, R3 assembly, R4 plan check.
"""

import hashlib
import re
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple

from src.docs.knowledge_index import (KNOWLEDGE_ROOT, SEVERITY_ORDER, check_fragments,
                                      load_fragments, registry_datasets)
from src.models.prompt import PromptTemplates

_SECTION_KEYS = (
    ("key points", "key_points"),
    ("implications for this project", "implications"),
    ("checked in the local files", "local_checks"),
    ("common mistakes", "common_mistakes"),
    ("source excerpts", "source_excerpts"),
)

_SECTION_BUNDLES = {"implications": ("local_checks", "implications")}

_CITATION = re.compile(r"\[([a-z0-9]+(?:-[a-z0-9]+)+)\]")
_SELECTION = re.compile(r"<selected_knowledge>(.*?)</selected_knowledge>", re.S)
_NOTES = re.compile(r"<knowledge_notes>(.*?)</knowledge_notes>", re.S)
_NOTE_HEADER = re.compile(r"^\[([a-z0-9]+(?:-[a-z0-9]+)+)\]\s*(.*)$")
_NOT_APPLICABLE = re.compile(r"^(?:\((?:must|should|info)\)\s*)?(?:not applicable|不适用)\s*[:：]?\s*(.*)$", re.I)
_NUMBER = re.compile(r"(?<![A-Za-z0-9_.])\d+(?:[.,]\d+)*(?:[eE][-+]?\d+)?(?![A-Za-z_\d])")
_BACKTICKED = re.compile(r"`([^`\n]+)`")
_IDENTIFIER = re.compile(r"(?<![0-9.])[A-Za-z_][A-Za-z0-9_]*")
_NOTE_LINE = re.compile(r"^\s*[-*•]\s*(?:Note|注意)?\s*\[([a-z0-9]+(?:-[a-z0-9]+)+)\]", re.M | re.I)
_CJK = re.compile(r"[\u3400-\u9fff]")

DEFAULT_SEVERITY_WEIGHTS = {"must": 1.0, "should": 0.6, "info": 0.3}
DEFAULT_SECTIONS = {
    "must": ["key_points", "implications", "common_mistakes"],
    "should": ["implications", "common_mistakes"],
    "info": ["summary", "implications"],
}


def _split_sections(body: str) -> Dict[str, str]:
    """A fragment body as {section key: text}; text before the first heading is 'preamble'."""
    sections: Dict[str, str] = {}
    current, lines = "preamble", []

    def flush():
        text = "\n".join(lines).strip()
        if text:
            sections[current] = text

    for line in body.splitlines():
        if line.startswith("## "):
            flush()
            heading = line[3:].strip().lower()
            current = next((key for prefix, key in _SECTION_KEYS if heading.startswith(prefix)),
                           heading.replace(" ", "_"))
            lines = []
        else:
            lines.append(line)
    flush()
    return sections


def _term_pattern(term: str, case_sensitive: bool = False) -> re.Pattern:
    """
    A pattern that finds `term` as a whole word in ASCII text, or as a plain
    substring when it contains CJK characters, which carry no word boundaries.
    """
    flags = 0 if case_sensitive else re.I
    if _CJK.search(term):
        return re.compile(re.escape(term), flags)
    return re.compile(rf"(?<![A-Za-z0-9_]){re.escape(term)}(?![A-Za-z0-9_])", flags)


class KnowledgeBase:
    """The fragments under one knowledge root, parsed once, with their R1-R4 operations."""

    def __init__(self, root: Path = KNOWLEDGE_ROOT, severity_weights: Optional[Dict] = None,
                 backtick_only_variables: Iterable[str] = (),
                 dataset_aliases: Optional[Dict[str, List[str]]] = None,
                 always_allowed_numbers: Iterable[float] = ()):
        self.root = Path(root)
        
        self.always_allowed_numbers = {float(n) for n in always_allowed_numbers or ()}
        self.fragments = load_fragments(self.root)
        for fragment in self.fragments:
            fragment["sections"] = _split_sections(fragment["body"])
        self.by_id = {f["id"]: f for f in self.fragments}

        self.severity_weights = dict(DEFAULT_SEVERITY_WEIGHTS)
        self.severity_weights.update(severity_weights or {})
        self.backtick_only = set(backtick_only_variables or ())

        aliases = dataset_aliases or {}
        self.dataset_patterns = {
            name: [_term_pattern(term) for term in [name, *(aliases.get(name) or [])]]
            for name in registry_datasets()
        }
        self._keyword_patterns = {
            f["id"]: [(str(k), _term_pattern(str(k))) for k in (f.get("keywords") or [])]
            for f in self.fragments
        }

    def check(self) -> List[str]:
        """The catalogue checks of src.docs.knowledge_index, run on this root."""
        return check_fragments(self.fragments, self.root)

    def digest(self) -> Dict[str, Any]:
        """A fingerprint of every fragment and of the index, for the run's manifest."""
        files = {f["path"].relative_to(self.root).as_posix():
                 hashlib.sha256(f["path"].read_bytes()).hexdigest() for f in self.fragments}
        index = self.root / "INDEX.md"
        combined = hashlib.sha256("".join(f"{k}:{v}\n" for k, v in sorted(files.items()))
                                  .encode()).hexdigest()
        return {
            "root": str(self.root),
            "fragment_count": len(files),
            "fragments_sha256": combined,
            "index_sha256": hashlib.sha256(index.read_bytes()).hexdigest() if index.exists() else None,
        }

    def fragment_record(self, fragment_id: str) -> Dict[str, Any]:
        fragment = self.by_id[fragment_id]
        return {
            "id": fragment_id,
            "severity": fragment.get("severity"),
            "path": fragment["path"].relative_to(self.root).as_posix(),
            "sha256": hashlib.sha256(fragment["path"].read_bytes()).hexdigest(),
        }

    def _signals(self, fragment: Dict[str, Any], text: str, quoted: set) -> Dict[str, List[str]]:
        """What in `text` points at this fragment."""
        variables = []
        for name in fragment.get("variables") or []:
            name = str(name)
            if name in quoted:
                variables.append(name)
            # A short name that is also an English word ('so', 'ph') counts only
            # when the question marks it as a name.
            elif name not in self.backtick_only and _term_pattern(name, case_sensitive=True).search(text):
                variables.append(name)

        keywords = [k for k, pattern in self._keyword_patterns[fragment["id"]] if pattern.search(text)]
        models = [str(m) for m in fragment.get("models") or []
                  if _term_pattern(str(m)).search(text)]
        datasets = [d for d in fragment.get("datasets") or []
                    if any(p.search(text) for p in self.dataset_patterns.get(d, []))]
        return {"variables": variables, "keywords": keywords, "models": models, "datasets": datasets}

    def prefilter(self, query: str, context: str = "",
                  max_candidates: Optional[int] = None) -> List[Dict[str, Any]]:
        """
        Every fragment scored against the question, best first, cut to max_candidates.
        """
        text = f"{context}\n{query}" if context else query
        quoted = set(re.findall(r"`([^`\s]+)`", text))
        candidates = []
        for fragment in self.fragments:
            signals = self._signals(fragment, text, quoted)
            raw = (2 * len(signals["variables"]) + len(signals["keywords"])
                   + bool(signals["datasets"]) + bool(signals["models"]))
            weight = self.severity_weights.get(fragment.get("severity"), 0.0)
            candidates.append({
                "id": fragment["id"],
                "type": fragment.get("type"),
                "severity": fragment.get("severity"),
                "summary": " ".join(str(fragment.get("summary", "")).split()),
                "score": round(raw * weight, 2),
                "signals": {k: v for k, v in signals.items() if v},
            })
        candidates.sort(key=lambda c: (-c["score"], SEVERITY_ORDER.get(c["severity"], 9), c["id"]))
        return candidates[:max_candidates] if max_candidates else candidates

    @staticmethod
    def lexical_top(candidates: List[Dict[str, Any]], k: int) -> List[Dict[str, str]]:
        """R1's own choice: the k best candidates that matched anything at all."""
        return [{"id": c["id"], "reason": "lexical prefilter score %s" % c["score"]}
                for c in candidates if c["score"] > 0][:k]

    def parse_selection(self, text: str, max_selected: int) -> Optional[List[Dict[str, str]]]:
        """
        The fragments named inside <selected_knowledge>, as [{"id", "reason"}].
        """
        match = _SELECTION.search(text or "")
        if not match:
            return None
        chosen, seen = [], set()
        for line in match.group(1).splitlines():
            line = line.strip().lstrip("-*•").strip()
            if not line:
                continue
            # A Chinese answer separates the reason with a full-width colon.
            fragment_id, _, reason = line.replace("：", ":").partition(":")
            fragment_id = fragment_id.strip().strip("`[]")
            if fragment_id in self.by_id and fragment_id not in seen:
                seen.add(fragment_id)
                chosen.append({"id": fragment_id, "reason": reason.strip()})
        return chosen[:max_selected]

    def render_fragment(self, fragment_id: str, sections: List[str]) -> str:
        """One fragment as a block headed by its id, holding only the given sections."""
        fragment = self.by_id[fragment_id]
        parts = [f"### [{fragment_id}] {fragment.get('title', '')}".rstrip(),
                 f"(type: {fragment.get('type')}; severity: {fragment.get('severity')})"]
        body_sections = fragment["sections"]
        wants_body = any(name != "summary" for name in sections)
        
        if wants_body and body_sections.get("preamble"):
            parts.append(body_sections["preamble"])
        for name in sections:
            if name == "summary":
                parts.append("Summary: " + " ".join(str(fragment.get("summary", "")).split()))
                continue
            for key in _SECTION_BUNDLES.get(name, (name,)):
                if body_sections.get(key):
                    title = next((t for t, k in _SECTION_KEYS if k == key), key).capitalize()
                    parts.append(f"#### {title}\n{body_sections[key]}")
        return "\n".join(parts)

    def assemble(self, selected: List[Dict[str, str]], sections_by_severity: Dict[str, List[str]],
                 max_chars: int) -> Tuple[str, Dict[str, str]]:
        """
        The selected fragments as one text block within max_chars, must first.
        """
        sections = dict(DEFAULT_SECTIONS)
        sections.update(sections_by_severity or {})
        ordered = sorted(selected, key=lambda s: SEVERITY_ORDER.get(self.by_id[s["id"]].get("severity"), 9))

        blocks, used, placement = [], 0, {}
        for item in ordered:
            severity = self.by_id[item["id"]].get("severity")
            for form, names in (("full", sections.get(severity, ["summary"])), ("summary_only", ["summary"])):
                block = self.render_fragment(item["id"], names)
                if not max_chars or used + len(block) + 2 <= max_chars:
                    blocks.append(block)
                    used += len(block) + 2
                    placement[item["id"]] = form
                    break
            else:
                placement[item["id"]] = "dropped"
        return "\n\n".join(blocks), placement

    def parse_notes(self, text: str, expected_ids: List[str]) -> Optional[Dict[str, Dict[str, Any]]]:
        """
        The condensed notes inside <knowledge_notes>, as {id: {"status", "notes", "reason"}}.
        """
        match = _NOTES.search(text or "")
        if not match:
            return None
        expected = set(expected_ids)
        parsed: Dict[str, Dict[str, Any]] = {}
        current = None
        for raw in match.group(1).splitlines():
            line = raw.strip()
            if not line:
                continue
            header = _NOTE_HEADER.match(line)
            if header:
                fragment_id, rest = header.group(1), header.group(2).strip()
                current = fragment_id if fragment_id in expected else None
                if current is None:
                    continue
                entry = parsed.setdefault(current, {"status": "notes", "notes": [], "reason": ""})
                not_applicable = _NOT_APPLICABLE.match(rest)
                if not_applicable:
                    entry.update(status="not_applicable", reason=not_applicable.group(1).strip())
                continue
            if current is None:
                continue
            entry = parsed[current]
            if line[:1] in "-*•":
                entry["notes"].append(line.lstrip("-*•").strip())
            elif entry["notes"]:
                entry["notes"][-1] += " " + line
        return parsed

    def source_text(self, fragment_id: str) -> str:
        """The whole fragment file, front matter included: what a note may quote from."""
        return self.by_id[fragment_id]["path"].read_text(encoding="utf-8")

    @staticmethod
    def _numbers(text: str) -> set:
        values = set()
        for token in _NUMBER.findall(text):
            readings = [token.replace(",", "")] + (token.split(",") if "," in token else [])
            for reading in readings:
                try:
                    value = float(reading)
                except ValueError:
                    continue
                if value != int(value) or value >= 10:
                    values.add(value)
        return values

    def verify_notes(self, fragment_id: str, notes: List[str], question: str,
                     sibling_ids: Iterable[str] = ()) -> List[str]:
        source = self.source_text(fragment_id) + "\n" + (question or "")
        name_source = "\n".join([source] + [self.source_text(i) for i in sibling_ids
                                            if i != fragment_id and i in self.by_id])
        allowed_numbers = self._numbers(source) | self.always_allowed_numbers
        problems = []
        for note in notes:
            for value in sorted(self._numbers(note) - allowed_numbers):
                problems.append(f"number {value:g} is not in the fragment or the question")
            # A backticked span may be an expression built from quoted names
            # (`grid_label=gr`, `Σ tos·areacello / Σ areacello`), so each name in it
            # is checked, not the span as written.
            for span in _BACKTICKED.findall(note):
                for name in _IDENTIFIER.findall(span):
                    if not _term_pattern(name, case_sensitive=True).search(name_source):
                        problems.append(f"`{name}` is not in the selected fragments or the question")
        return list(dict.fromkeys(problems))

    def render_notes(self, fragment_id: str, notes: List[str]) -> str:
        severity = self.by_id[fragment_id].get("severity")
        return "\n".join([f"[{fragment_id}] ({severity})", *(f"- {note}" for note in notes)])

    def render_summary(self, fragment_id: str) -> str:
        fragment = self.by_id[fragment_id]
        summary = " ".join(str(fragment.get("summary", "")).split())
        return f"[{fragment_id}] ({fragment.get('severity')})\n- {summary}"

    def assemble_blocks(self, blocks: List[Tuple[str, str]], max_chars: int) -> Tuple[str, List[str]]:
        """
        Pack (fragment id, block) pairs into max_chars, must first; returns the text
        and the ids that did not fit.
        """
        ordered = sorted(blocks, key=lambda b: SEVERITY_ORDER.get(self.by_id[b[0]].get("severity"), 9))
        kept, used, dropped = [], 0, []
        for fragment_id, block in ordered:
            if not max_chars or used + len(block) + 2 <= max_chars:
                kept.append(block)
                used += len(block) + 2
            else:
                dropped.append(fragment_id)
        return "\n\n".join(kept), dropped

    def check_plan(self, plan_text: str, selected_ids: List[str]) -> Dict[str, Any]:
        """
        Which selected fragments the plan cites as [fragment-id], which it passed
        over, and which it wrote a note for more than once.
        """
        plan_text = plan_text or ""
        cited = list(dict.fromkeys(_CITATION.findall(plan_text)))
        selected = list(selected_ids)
        uncited = [i for i in selected if i not in cited]
        note_lines: Dict[str, int] = {}
        for fragment_id in _NOTE_LINE.findall(plan_text):
            note_lines[fragment_id] = note_lines.get(fragment_id, 0) + 1
        return {
            "cited": [i for i in cited if i in selected],
            "uncited": uncited,
            "uncited_must": [i for i in uncited if self.by_id[i].get("severity") == "must"],
            "cited_not_selected": [i for i in cited if i not in selected and i in self.by_id],
            "note_lines": note_lines,
            "repeated_notes": {i: n for i, n in note_lines.items() if n > 1},
        }


# Statuses of a selected fragment under which it reached the planner.
INJECTED_STATUSES = ("condensed", "summary_fallback", "full", "summary_only")

# ask(name, system, human) -> (answer text, seconds): one recorded model exchange.
Ask = Callable[[str, str, str], Tuple[str, float]]


class KnowledgeRetriever:
    """
    The retrieval run for one question: R1 prefilter, R2a selection, R2b condensing,
    R3 assembly. `settings` is the resolved 'knowledge' section of the agent's
    configuration (see AgentWithTools._setup_knowledge).
    """

    def __init__(self, knowledge_base: KnowledgeBase, settings: Dict[str, Any]):
        self.kb = knowledge_base
        self.settings = settings

    def retrieve(self, query: str, context: str, language: str, ask: Ask) -> Tuple[str, Dict[str, Any]]:
        """
        The text for the planning prompt ("" when nothing applies) and the record
        written to the turn's knowledge.json.
        """
        kb, cfg = self.kb, self.settings
        candidates = kb.prefilter(query, context, cfg["max_candidates"])
        record: Dict[str, Any] = {
            "enabled": True,
            "mode": cfg["selector_mode"],
            "condensing": cfg["selector_mode"] == "llm" and cfg["condense_enabled"],
            "knowledge_base": kb.digest(),
            "candidates": [{k: c[k] for k in ("id", "severity", "score", "signals")} for c in candidates],
        }

        selected = self._select(query, context, language, candidates, ask, record)
        if not selected:
            record.update(selected=[], injected_ids=[], injected_chars=0)
            return "", record

        if record["condensing"]:
            text, entries = self._condense(query, context, language, selected, ask, record)
        else:
            text, placement = kb.assemble(selected, cfg["sections"], cfg["max_chars"])
            entries = [{**kb.fragment_record(item["id"]), "reason": item["reason"],
                        "status": placement[item["id"]]} for item in selected]

        record["selected"] = entries
        record["injected_ids"] = [e["id"] for e in entries if e["status"] in INJECTED_STATUSES]
        record["injected_chars"] = len(text)
        return text, record

    def _ask_until_readable(self, ask: Ask, name: str, system: str, human: str, tag: str,
                            parse: Callable[[str], Any]) -> Tuple[Any, Dict[str, Any]]:
        """
        One exchange, asked again up to `retries` times while the answer cannot be
        read: the call fails, the answer is empty (a reasoning model sometimes
        returns nothing at all), or the tag is missing.

        Returns what `parse` made of the first readable answer, or None, and
        {"attempts": [...], "seconds": total, "fallback_reason": last problem}.
        """
        attempts: List[Dict[str, Any]] = []
        for attempt in range(1 + max(0, int(self.settings.get("retries", 0)))):
            label = name if attempt == 0 else f"{name}_retry{attempt}"
            try:
                answer, seconds = ask(label, system, human)
            except Exception as e:
                attempts.append({"problem": f"the call failed: {type(e).__name__}: {e}"})
                continue
            parsed = parse(answer)
            if parsed is not None:
                attempts.append({"seconds": round(seconds, 1)})
                break
            attempts.append({"seconds": round(seconds, 1),
                             "problem": "the answer is empty" if not (answer or "").strip()
                             else f"the answer has no <{tag}> tag"})
        else:
            parsed = None
        info: Dict[str, Any] = {"attempts": attempts,
                                "seconds": round(sum(a.get("seconds") or 0 for a in attempts), 1)}
        if parsed is None:
            info["fallback_reason"] = f"{name}: {attempts[-1]['problem']} (after {len(attempts)} attempts)"
        return parsed, info

    def _select(self, query: str, context: str, language: str, candidates: List[Dict[str, Any]],
                ask: Ask, record: Dict[str, Any]) -> List[Dict[str, str]]:
        kb, cfg = self.kb, self.settings
        info: Dict[str, Any] = {}
        selected = None
        if cfg["selector_mode"] == "llm":
            system = PromptTemplates.get_knowledge_selection_prompt(
                language, candidates, context, cfg["max_selected"])
            human = PromptTemplates.get_knowledge_selection_input(language, query)
            selected, info = self._ask_until_readable(
                ask, "knowledge_selection", system, human, "selected_knowledge",
                lambda answer: kb.parse_selection(answer, cfg["max_selected"]))
        if selected is None:
            selected = kb.lexical_top(candidates, cfg["fallback_top_k"])
        info["ids"] = [item["id"] for item in selected]
        record["selection"] = info
        return selected

    def _condense(self, query: str, context: str, language: str, selected: List[Dict[str, str]],
                  ask: Ask, record: Dict[str, Any]) -> Tuple[str, List[Dict[str, Any]]]:
        kb, cfg = self.kb, self.settings
        ids = [item["id"] for item in selected]
        sources = "\n\n".join(
            kb.render_fragment(i, cfg["source_sections"].get(kb.by_id[i].get("severity"), ["summary"]))
            for i in ids)
        system = PromptTemplates.get_knowledge_condense_prompt(
            language, sources, context, cfg["max_notes_per_fragment"], cfg["condense_max_chars"])
        human = PromptTemplates.get_knowledge_condense_input(language, query)

        parsed, info = self._ask_until_readable(
            ask, "knowledge_condensing", system, human, "knowledge_notes",
            lambda answer: kb.parse_notes(answer, ids))
        info["source_chars"] = len(sources)
        record["condensing_run"] = info

        question = f"{context}\n{query}" if context else query
        entries, blocks = [], []
        for item in selected:
            fragment_id = item["id"]
            entry = {**kb.fragment_record(fragment_id), "reason": item["reason"]}
            result = parsed.get(fragment_id) if parsed is not None else None
            if parsed is None:
                entry.update(status="summary_fallback", fallback=info["fallback_reason"])
            elif result is None:
                entry.update(status="summary_fallback", fallback="the condenser did not mention it")
            elif result["status"] == "not_applicable":
                entry.update(status="not_applicable", not_applicable_reason=result["reason"])
            else:
                notes = result["notes"][:cfg["max_notes_per_fragment"]]
                if len(result["notes"]) > len(notes):
                    entry["notes_cut"] = len(result["notes"]) - len(notes)
                problems = (kb.verify_notes(fragment_id, notes, question, sibling_ids=ids)
                            if notes else ["no notes"])
                if problems:
                    entry.update(status="summary_fallback", fallback="the notes failed the check",
                                 violations=problems, rejected_notes=notes)
                else:
                    entry.update(status="condensed", notes=notes)

            entries.append(entry)
            if entry["status"] == "condensed":
                blocks.append((fragment_id, kb.render_notes(fragment_id, entry["notes"])))
            elif entry["status"] == "summary_fallback":
                blocks.append((fragment_id, kb.render_summary(fragment_id)))

        text, dropped = kb.assemble_blocks(blocks, cfg["condense_max_chars"])
        for entry in entries:
            if entry["id"] in dropped:
                entry["status"] = "dropped"
        return text, entries
