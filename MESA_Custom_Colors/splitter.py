#!/usr/bin/env python3

import re
import sys

STRUCTURAL_STARTS = (
    "\\begin", "\\end", "\\label", "\\includegraphics", "\\centering",
    "\\hline", "\\toprule", "\\midrule", "\\bottomrule",
)

PROTECTED_COMMANDS = (
    "texttt", "url", "href", "label", "ref", "autoref", "cref", "Cref",
    "cite", "citep", "citet", "citealt", "includegraphics",
    "section", "subsection", "subsubsection", "paragraph",
)

ABBREVIATIONS = (
    "e.g.", "i.e.", "Fig.", "Figs.", "Eq.", "Eqs.", "Sec.", "Secs.",
    "Dr.", "Prof.", "Mr.", "Ms.", "et al.",
)

FILE_EXTENSIONS = (
    ".bin", ".tex", ".pdf", ".png", ".jpg", ".jpeg", ".dat", ".data",
    ".txt", ".fits", ".csv", ".h5", ".hdf5", ".yaml", ".yml",
)

def is_structural_block(block: str) -> bool:
    stripped = block.strip()
    if not stripped:
        return True

    lines = stripped.splitlines()

    if any(line.strip().startswith(STRUCTURAL_STARTS) for line in lines):
        return True

    if "\\begin{" in stripped or "\\end{" in stripped:
        return True

    if "&" in stripped and "\\\\" in stripped:
        return True

    if stripped.startswith("\\[") or stripped.endswith("\\]"):
        return True

    if stripped.startswith("$$") or stripped.endswith("$$"):
        return True

    return False


def protect_latex_fragments(text: str):
    placeholders = {}

    def repl(match):
        key = f"@@PROTECTED_{len(placeholders)}@@"
        placeholders[key] = match.group(0)
        return key

    # Protect common one-brace commands: \texttt{flux_cube.bin}, \ref{...}, etc.
    command_pattern = r"\\(" + "|".join(PROTECTED_COMMANDS) + r")(\[[^\]]*\])?\{[^{}]*\}"
    text = re.sub(command_pattern, repl, text)

    # Protect simple filenames/extensions that appear outside commands.
    for ext in FILE_EXTENSIONS:
        pattern = rf"\b[\w./-]+{re.escape(ext)}\b"
        text = re.sub(pattern, repl, text)

    # Protect common abbreviations.
    for abbr in ABBREVIATIONS:
        text = text.replace(abbr, repl(re.match(r".+", abbr)))

    return text, placeholders


def unprotect(text: str, placeholders: dict):
    for key, value in placeholders.items():
        text = text.replace(key, value)
    return text


def split_sentences(paragraph: str):
    paragraph = " ".join(line.strip() for line in paragraph.splitlines())
    paragraph = re.sub(r"\s+", " ", paragraph).strip()

    protected, placeholders = protect_latex_fragments(paragraph)

    # Split after sentence punctuation only when followed by likely new sentence.
    parts = re.split(r"(?<=[.!?])\s+(?=[A-Z\\])", protected)

    parts = [unprotect(part.strip(), placeholders) for part in parts if part.strip()]
    return parts


def format_block(block: str):
    if is_structural_block(block):
        return block.strip()

    sentences = split_sentences(block)

    if len(sentences) <= 1:
        return sentences[0] if sentences else ""

    return "\n".join(sentences)


def main():
    text = sys.stdin.read() if len(sys.argv) == 1 else open(sys.argv[1]).read()

    blocks = re.split(r"\n\s*\n", text)
    formatted = [format_block(block) for block in blocks]

    print("\n\n".join(block for block in formatted if block.strip()))


if __name__ == "__main__":
    main()
