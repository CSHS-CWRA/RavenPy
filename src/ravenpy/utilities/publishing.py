"""Publishing utilities for RavenPy."""

from __future__ import annotations

import os
import re
from io import StringIO
from pathlib import Path
from typing import TextIO


def publish_release_notes(
    style: str = "md",
    file: os.PathLike[str] | StringIO | TextIO | None = None,
    changes: str | os.PathLike[str] | None = None,
) -> str | None:
    """
    Format release notes in Markdown or ReStructuredText.

    Parameters
    ----------
    style : {"rst", "md"}
        Use ReStructuredText formatting or Markdown. Default: Markdown.
    file : {os.PathLike, StringIO, TextIO}, optional
        If provided, prints to the given file-like object. Otherwise, returns a string.
    changes : str or os.PathLike[str], optional
        If provided, manually points to the file where the changelog can be found.
        Assumes a relative path otherwise.

    Returns
    -------
    str, optional
        If `file` not provided, the formatted release notes.

    Notes
    -----
    This function is used solely for development and packaging purposes.
    """
    if isinstance(changes, str | Path):
        changes_file = Path(changes).absolute()
    else:
        changes_file = Path(__file__).absolute().parents[3].joinpath("CHANGELOG.rst")

    if not changes_file.exists():
        raise FileNotFoundError("Changelog file not found in RavenPy folder tree.")

    with Path(changes_file).open(encoding="utf-8") as hf:
        changes = hf.read()

    if style == "rst":
        hyperlink_replacements = {
            r":issue:`([0-9]+)`": r"`GH/\1 <https://github.com/CSHS-CWRA/RavenPy/issues/\1>`_",
            r":pull:`([0-9]+)`": r"`PR/\1 <https://github.com/CSHS-CWRA/RavenPy/pull/\>`_",
            r":user:`([a-zA-Z0-9_.-]+)`": r"`@\1 <https://github.com/\1>`_",
        }
    elif style == "md":
        hyperlink_replacements = {
            r":issue:`([0-9]+)`": r"[GH/\1](https://github.com/CSHS-CWRA/RavenPy/issues/\1)",
            r":pull:`([0-9]+)`": r"[PR/\1](https://github.com/CSHS-CWRA/RavenPy/pull/\1)",
            r":user:`([a-zA-Z0-9_.-]+)`": r"[@\1](https://github.com/\1)",
        }
    else:
        msg = f"Formatting style not supported: {style}"
        raise NotImplementedError(msg)

    for search, replacement in hyperlink_replacements.items():
        changes = re.sub(search, replacement, changes)

    if style == "md":
        changes = changes.replace("=========\nChangelog\n=========", "# Changelog")

        titles = {r"\n(.*?)\n([\-]{1,})": "-", r"\n(.*?)\n([\^]{1,})": "^"}
        for title_expression, level in titles.items():
            found = re.findall(title_expression, changes)
            for grouping in found:
                fixed_grouping = (
                    str(grouping[0]).replace("(", r"\(").replace(")", r"\)")
                )
                search = rf"({fixed_grouping})\n([\{level}]{'{' + str(len(grouping[1])) + '}'})"
                replacement = f"{'##' if level == '-' else '###'} {grouping[0]}"
                changes = re.sub(search, replacement, changes)

        link_expressions = r"[\`]{1}([\w\s]+)\s<(.+)>`\_"
        found = re.findall(link_expressions, changes)
        for grouping in found:
            search = rf"`{grouping[0]} <.+>`\_"
            replacement = f"[{str(grouping[0]).strip()}]({grouping[1]})"
            changes = re.sub(search, replacement, changes)

    if not file:
        return changes
    if isinstance(file, Path | os.PathLike):
        with Path(file).open("w", encoding="utf-8") as f:
            print(changes, file=f)
    else:
        print(changes, file=file)
    return None
