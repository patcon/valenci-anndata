import asyncio
import warnings
import re
import textwrap
from anndata import AnnData
import pandas as pd
from dataclasses import dataclass
from googletrans import Translator
from io import StringIO
from pathlib import Path
from polis_client import PolisClient
from typing import List, Literal, Optional
from urllib.parse import urlparse

from ..preprocessing import rebuild_vote_matrix
from ..utils import run_async


DEFAULT_BASE = "https://pol.is"

REPORT_RE = re.compile(r"^r[a-z0-9]{15,}$")     # e.g. r4zdxrdscmukmkakmbz3k
CONVO_RE  = re.compile(r"^[0-9][a-z0-9]{8,}$")  # e.g. 4asymkcrjf (starts with digit)

@dataclass
class PolisSource:
    kind: Literal["api", "report", "local"]
    base_url: str | None = None
    conversation_id: str | None = None
    report_id: str | None = None
    path: Path | None = None

def _parse_polis_source(source: str):
    """
    Returns a PolisSource with:
        base_url
        report_id
        conversation_id
    """
    # ───────────────────────────────────────────
    # 0. Local directory case
    # ───────────────────────────────────────────
    path = Path(source)
    if path.exists() and path.is_dir():
        try:
            _find_single_csv(path, "votes.csv")
            _find_single_csv(path, "comments.csv")
        except (FileNotFoundError, ValueError) as e:
            raise ValueError(str(e)) from None

        return PolisSource(
            kind="local",
            path=path,
        )

    source = source.strip()

    # ───────────────────────────────────────────
    # 1. URL case
    # ───────────────────────────────────────────
    if source.startswith("http://") or source.startswith("https://"):
        url = urlparse(source)
        base_url = f"{url.scheme}://{url.netloc}"

        # normalize path segments
        parts = [p for p in url.path.split("/") if p]

        report_id = None
        conversation_id = None

        if len(parts) == 2 and parts[0] == "report":
            # /report/<report_id>
            report_id = parts[1]
        elif len(parts) == 1:
            # /<conversation_id>
            conversation_id = parts[0]

        return PolisSource(
            kind="report" if report_id else "api",
            base_url=base_url,
            report_id=report_id,
            conversation_id=conversation_id,
        )

    # ───────────────────────────────────────────
    # 2. Bare IDs (conversation or report)
    # ───────────────────────────────────────────
    # Starts with digit → conversation_id
    if CONVO_RE.match(source):
        return PolisSource(
            kind="api",
            base_url=DEFAULT_BASE,
            report_id=None,
            conversation_id=source,
        )

    # Starts with "r" → report_id
    if REPORT_RE.match(source):
        return PolisSource(
            kind="report",
            base_url=DEFAULT_BASE,
            report_id=source,
            conversation_id=None,
        )

    raise ValueError(f"Unrecognized Polis source format: {source}")

def _fill_missing_fields_from_api(statements: pd.DataFrame, conversation_id: str, client: PolisClient, fields=("is-seed", "is-meta")) -> pd.DataFrame:
    """
    Fill missing columns in `statements` from the API, only for specified fields.
    Warns the user if any fields are filled.
    """
    missing_fields = [col for col in fields if col not in statements.columns or statements[col].isnull().all()]
    if not missing_fields:
        return statements

    api_statements = client.get_comments(conversation_id=conversation_id)
    if not api_statements:
        warnings.warn(f"Missing fields {', '.join(missing_fields)} but API fetch returned no statements.")
        return statements

    api_df = pd.DataFrame([s.to_dict() for s in api_statements]).rename(columns={
        "tid": "comment-id",
        "pid": "author-id",
        "txt": "comment-body",
        "created": "timestamp",
        "mod": "moderated",
        "is_seed": "is-seed",
        "is_meta": "is-meta",
    }).set_index("comment-id", drop=False)

    # Fill only the missing fields
    for col in missing_fields:
        if col in api_df.columns:
            statements[col] = statements[col].combine_first(api_df[col])

    warnings.warn(
        f"The following fields were missing and have been filled from the API: {', '.join(missing_fields)}"
    )

    return statements

def _find_single_csv(path: Path, suffix: str) -> Path:
    """
    Resolve a CSV file from a directory, preferring an exact filename match.

    If `path / suffix` exists, it is returned immediately. Otherwise, the
    directory is searched for files matching `*{suffix}`. An error is raised
    if zero or multiple matches are found.
    """
    # 1. Prefer exact filename within directory
    exact = path / suffix
    if exact.is_file():
        return exact

    # 2. Fall back to suffix-based search
    matches = sorted(path.glob(f"*{suffix}"))

    if not matches:
        raise FileNotFoundError(f"No *{suffix} file found in {path}")

    if len(matches) > 1:
        raise ValueError(
            f"Multiple *{suffix} files found in {path}: "
            + ", ".join(p.name for p in matches)
        )

    return matches[0]


def load(source: str, *, translate_to: Optional[str] = None, build_X: bool = True) -> AnnData:
    """
    Load a Polis conversation or report into an AnnData object.

    This function accepts either a URL or an ID for a Polis conversation or report,
    fetches raw vote events and statements via the Polis API or CSV export, and
    optionally constructs a participant × statement vote matrix in `adata.X`.

    Parameters
    ----------
    source : str
        The Polis source to load. Supported formats include:

        - Full report URL: `https://pol.is/report/<report_id>`
        - Conversation URL: `https://pol.is/<conversation_id>`
        - Custom host URLs: `https://<host>/report/<report_id>` or `https://<host>/<conversation_id>`
        - Bare IDs:
            - Conversation ID (starts with a digit), e.g., `4asymkcrjf`
            - Report ID (starts with 'r'), e.g., `r4zdxrdscmukmkakmbz3k`
        - Local directory containing CSV exports:
            - *votes.csv
            - *comments.csv

        The function will automatically parse the source to determine whether
        it refers to a conversation or report and fetch the appropriate data.


    translate_to : str or None, optional
        Target language code (e.g., "en", "fr", "es") for translating statement text.
        If provided, the original statement text in `adata.uns["statements"]["comment-body"]`
        is translated and stored in `adata.var["content"]`. The `adata.var["language_current"]`
        field is updated to the target language, and `adata.var["is_translated"]` is set to True.
        Defaults to None (no translation).

    build_X : bool, default True
        If True, constructs a participant × statement vote matrix from the raw
        votes using `rebuild_vote_matrix()`. This populates `adata.obs`,
        `adata.var`, and `adata.X` (with a copy in
        `adata.layers['raw_sparse']`). After the first build, a snapshot of this
        initial matrix is stored in `adata.raw`.

    Returns
    -------
    adata : anndata.AnnData
        An AnnData object containing the loaded Polis data.

        
    pd.DataFrame
        `adata.uns["votes"]`  
        Raw vote events fetched from the API or CSV export.
    dict
        `adata.uns["votes_meta"]`  
        Metadata about the sources of votes, e.g., API vs CSV.
    pd.DataFrame
        `adata.uns["statements"]`  
        Raw statements/comments for the conversation.
    dict
        `adata.uns["statements_meta"]`  
        Metadata about the statements source.
    dict
        `adata.uns["source"]`  
        Basic information about the Polis source (base URL, conversation ID, report ID).
    dict
        `adata.uns["schema"]`  
        High-level description of `X` and `votes`.
    np.ndarray
        `adata.X` (if `build_X=True`)  
        Participant × statement vote matrix (rows = participants, columns = statements).
    pd.DataFrame 
        `adata.obs` (if `build_X=True`)  
        Participant metadata (index = voter IDs).
    pd.DataFrame 
        `adata.var` (if `build_X=True`)  
        Statement metadata (index = statement IDs).
    anndata.AnnData 
        `adata.raw` (if `build_X=True`)  
        Snapshot of the first vote matrix and associated metadata. This allows
        downstream filtering or processing without losing the original vote matrix.

    Notes
    -----
    - If `build_X=False`, only `adata.uns` will be populated, containing the raw
      votes and statements, and `.X`, `.obs`, `.var`, and `.raw` will remain empty.
    - `adata.raw` is assigned only after the first vote matrix build and is intended
      to be immutable.
    - If `translate_to` is provided, `adata.var["content"]` is updated with translated
    text and `adata.var["language_current"]` is set to the target language.
    - The vote matrix is derived from the most recent votes per participant per statement,
      sorted by timestamp.

    Examples
    --------

    Load data from a report or conversation ID or URL.

    ```py
    adata = val.datasets.polis.load("https://pol.is/report/r2dfw8eambusb8buvecjt")
    adata = val.datasets.polis.load("6rphtwwfn4")
    ```

    Load data from an alternative Polis instance via URL.

    ```py
    adata = val.datasets.polis.load("https://polis.tw/6rphtwwfn4")
    ```

    Load data from a path containing Polis CSV export files.

    ```sh
    $ ls exports/my_conversation_2024_11_03
    comments.csv votes.csv summary.csv ...
    ```

    ```py
    adata = val.datasets.polis.load("./exports/my_conversation_2024_11_03")
    ```
    """
    adata = _load_raw_polis_data(source)

    if build_X:
        rebuild_vote_matrix(adata, trim_rule=1.0, inplace=True)
        adata.raw = adata.copy()
        # Store a copy in case we bring something else into X workspace later.
        adata.layers["raw_sparse"] = adata.X # type: ignore[arg-type]

    _populate_var_statements(adata, translate_to=translate_to)

    # if convo_meta.conversation_id:
    #     xids = client.get_xids(conversation_id=convo_meta.conversation_id)
    #     adata.uns["xids"] = pd.DataFrame(xids)

    return adata

def _load_raw_polis_data(source):
    convo_src = _parse_polis_source(source)
    if convo_src.kind == "local":
        return _load_from_local_path(convo_src)

    if convo_src.kind in {"api", "report"}:
        return _load_from_polis(convo_src)

    raise AssertionError("Unreachable")

def _load_from_local_path(convo_src: PolisSource) -> AnnData:
    path = convo_src.path
    assert path is not None

    votes_path = _find_single_csv(path, "votes.csv")
    comments_path = _find_single_csv(path, "comments.csv")

    votes = pd.read_csv(votes_path)
    votes["source"] = "local_csv"
    votes["source_id"] = str(path)

    statements = pd.read_csv(comments_path)
    # TODO: detect if is-seed and is-meta are missing, and augment if not

    adata = AnnData()
    return _load_votes_and_statements(adata, votes, statements, convo_src)

def _load_from_polis(convo_src: PolisSource) -> AnnData:
    assert convo_src.base_url is not None
    client = PolisClient(base_url=convo_src.base_url)

    # ───────────────────────────────────────────
    # Load votes
    # ───────────────────────────────────────────
    if convo_src.report_id:
        votes_csv_text = client.get_export_file(filename="votes.csv", report_id=convo_src.report_id)
        votes = pd.read_csv(StringIO(votes_csv_text))
        votes["source"] = "csv"
        votes["source_id"] = convo_src.report_id
        votes.sort_values("timestamp", inplace=True)

        report = client.get_report(report_id=convo_src.report_id)
        if report:
            convo_src.conversation_id = report.conversation_id or None

    elif convo_src.conversation_id:
        votes_list = client.get_all_votes_slow(conversation_id=convo_src.conversation_id)
        votes = pd.DataFrame([v.to_dict() for v in votes_list])
        votes_rename_map = {
            "modified": "timestamp",
            "pid": "voter-id",
            "tid": "comment-id",
        }
        votes.rename(columns=votes_rename_map, inplace=True)
        votes["source"] = "api"
        votes["source_id"] = convo_src.conversation_id
        votes.sort_values("timestamp", inplace=True)

    else:
        raise ValueError("No votes could be loaded")

    # ───────────────────────────────────────────
    # Load statements
    # ───────────────────────────────────────────
    statements = client.get_comments(conversation_id=convo_src.conversation_id) or []
    statements_df = pd.DataFrame([s.to_dict() for s in statements])

    adata = AnnData()
    adata = _load_votes_and_statements(adata, votes, statements_df, convo_src)

    _maybe_print_attribution(convo_src)
    return adata

def _load_votes_and_statements(
    adata: AnnData,
    votes: pd.DataFrame,
    statements: pd.DataFrame,
    convo_src: PolisSource,
) -> AnnData:
    """
    Common path to populate adata.uns with votes and statements, normalizing columns
    and storing metadata.
    """
    # ───────────────────────────────────────────
    # Normalize votes
    # ───────────────────────────────────────────
    votes = votes.copy()
    votes.sort_values("timestamp", inplace=True)
    adata.uns["votes"] = votes

    adata.uns["votes_meta"] = {
        "sources": {
            convo_src.kind: {
                "via": "filesystem" if convo_src.kind == "local" else "live_api" if convo_src.kind == "api" else "live_csv",
                "conversation_id": convo_src.conversation_id,
                "report_id": convo_src.report_id,
                "path": str(convo_src.path) if convo_src.path else None,
                "base_url": convo_src.base_url,
                "retrieved_at": pd.Timestamp.utcnow().isoformat(),
            }
        },
        "sorted_by": "timestamp",
    }

    # ───────────────────────────────────────────
    # Normalize statements
    # ───────────────────────────────────────────
    statements = statements.copy()

    statements_rename_map = {
        "tid": "comment-id",
        "pid": "author-id",
        "txt": "comment-body",
        "created": "timestamp",
        "mod": "moderated",
        "is_seed": "is-seed",
        "is_meta": "is-meta",
    }
    statements.rename(columns={k: v for k, v in statements_rename_map.items() if k in statements.columns}, inplace=True)

    adata.uns["statements"] = statements.set_index("comment-id", drop=False)

    adata.uns["statements_meta"] = {
        "source": {
            "via": "filesystem" if convo_src.kind == "local" else "live_api" if convo_src.kind == "api" else "live_csv",
            "conversation_id": convo_src.conversation_id,
            "report_id": convo_src.report_id,
            "path": str(convo_src.path) if convo_src.path else None,
            "base_url": convo_src.base_url,
            "retrieved_at": pd.Timestamp.utcnow().isoformat(),
        }
    }

    adata.uns["source"] = {
        "kind": convo_src.kind,
        "base_url": convo_src.base_url,
        "conversation_id": convo_src.conversation_id,
        "report_id": convo_src.report_id,
        "path": str(convo_src.path) if convo_src.path else None,
    }

    adata.uns["schema"] = {
        "X": "participant × statement vote matrix (derived)",
        "votes": "raw vote events",
    }

    return adata

def _maybe_print_attribution(convo_src: PolisSource):
    """
    Print attribution text to satisfy Creative Commons license.
    """
    if not (convo_src.report_id or convo_src.conversation_id):
        return

    base = (
        "Data was gathered using the Polis software "
        "(see: https://compdemocracy.org/polis and "
        "https://github.com/compdemocracy/polis) "
        "and is sub-licensed under CC BY 4.0 with Attribution to "
        "The Computational Democracy Project."
    )

    if convo_src.report_id:
        tail = (
            "The data and more information about how the data was collected "
            f"can be found at the following link: https://pol.is/report/{convo_src.report_id}"
        )
    else:
        tail = (
            f"The data was retrieved from https://pol.is/{convo_src.conversation_id} "
            "and more information can be found at "
            "https://compdemocracy.org/Polis-Conversation-Data/"
        )

    print(format_attribution(base + "\n\n" + tail))

def format_attribution(text: str, *, width: int = 80) -> str:
    return "\n".join(
        textwrap.fill(
            paragraph,
            width=width,
            break_long_words=False,
            break_on_hyphens=False,
        )
        for paragraph in text.split("\n\n")
    )


def _populate_var_statements(adata, translate_to: Optional[str] = None):
    # Statements in adata.uns["statements"] use the Polis CSV export schema
    # (comment-id, author-id, comment-body, is-seed, is-meta, etc.)
    statements_aligned = adata.uns["statements"].copy()
    statements_aligned.index = statements_aligned.index.astype(str)
    statements_aligned = statements_aligned.reindex(adata.var_names)

    # Canonical Polis CSV-style fields
    adata.var["content"] = statements_aligned["comment-body"]
    adata.var["participant_id_authored"] = statements_aligned["author-id"]
    adata.var["created_date"] = statements_aligned["timestamp"]
    adata.var["moderation_state"] = statements_aligned["moderated"]

    # Optional fields (may or may not be present)
    adata.var["is_seed"] = (
        statements_aligned["is-seed"]
        if "is-seed" in statements_aligned.columns
        else None
    )

    adata.var["is_meta"] = (
        statements_aligned["is-meta"]
        if "is-meta" in statements_aligned.columns
        else pd.NA
    )

    adata.var["language_original"] = (
        statements_aligned["lang"]
        if "lang" in statements_aligned.columns
        else pd.NA
    )

    adata.var["language_current"] = adata.var["language_original"]
    adata.var["is_translated"] = False

    if translate_to is not None:
        translate_statements(adata, translate_to=translate_to, inplace=True)

# private async function
async def _translate_texts_async(texts: List[str], dest_lang: str) -> List[str]:
    async with Translator() as translator:
        results = await asyncio.gather(
            *(translator.translate(t, dest=dest_lang) for t in texts)
        )
    return [r.text for r in results]

# public synchronous wrapper
def translate_texts(texts: List[str], dest_lang: str) -> List[str]:
    return run_async(_translate_texts_async(texts, dest_lang))

def translate_statements(
    adata: AnnData,
    translate_to: Optional[str],
    inplace: bool = True
) -> Optional[list[str]]:
    """
    Translate statements in `adata.uns['statements']['comment-body']` into another language,
    or copy originals if translate_to is None.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing `uns['statements']` and `var_names`.
    translate_to : Optional[str]
        Target language code (e.g., "en", "fr", "es").
    inplace : bool, default True
        If True, updates `adata.var['content']` and `adata.var['language_current']`.
        If False, returns a list of translated strings without modifying `adata`.

    Returns
    -------
    translated_texts : list[str] | None
        List of translated texts if `inplace=False`, else None.
    """
    statements_aligned = adata.uns["statements"].copy()
    statements_aligned.index = statements_aligned.index.astype(str)
    statements_aligned = statements_aligned.reindex(adata.var_names)

    original_texts = statements_aligned["comment-body"].tolist()

    # ───────────────────────────────────────────
    # NO-TRANSLATION PATH (explicit)
    # ───────────────────────────────────────────
    if translate_to is None:
        if inplace:
            adata.var["content"] = original_texts
            adata.var["language_current"] = adata.var["language_original"]
            adata.var["is_translated"] = False
            return None
        else:
            return original_texts


    # ───────────────────────────────────────────
    # TRANSLATION PATH
    # ───────────────────────────────────────────
    translated_texts = run_async(
        _translate_texts_async(original_texts, translate_to)
    )

    if inplace:
        adata.var["content"] = translated_texts
        adata.var["language_current"] = translate_to
        adata.var["is_translated"] = True
        return None
    else:
        return translated_texts

__all__ = [
    "load",
    "translate_statements",
]