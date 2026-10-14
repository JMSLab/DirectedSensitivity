"""The markdown checker, and every markdown file in the repository against it."""
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT))
from DirectedSensitivity.MarkdownStyle import CheckMarkdown, FixSentences

MARKDOWN_FILES = sorted(ROOT / name for name in subprocess.run(
    ['git', 'ls-files', '*.md'], cwd=ROOT, capture_output=True, text=True).stdout.splitlines())


def messages(text):
    return [message for _, message in CheckMarkdown(text)]


@pytest.mark.parametrize('path', MARKDOWN_FILES, ids=lambda path: str(path.relative_to(ROOT)))
def test_every_markdown_file_follows_the_conventions(path):
    assert CheckMarkdown(path.read_text()) == []


@pytest.mark.parametrize('text', [
    'One sentence. Two sentences.\n',
    '- An item. With two sentences.\n',
    'The value is $x$. Then it grows.\n',
])
def test_flags_two_sentences_on_one_line(text):
    assert any('one sentence per line' in message for message in messages(text))


@pytest.mark.parametrize('text', [
    'One sentence.\nTwo sentences.\n',
    'Andrews, Essig Aberg, and Jesse M. Shapiro wrote it.\n',
    'As in DellaVigna et al. (2012), e.g. Table 2, cf. Figure 1.\n',
    '**Example**. The nesting model is $\\mathcal{N}$.\n',
    '## 1. Results\n',
    '| A. B | C. D |\n|---|---|\n',
    'See `a.py`. `b.py` too.\n',
    'Text $a$\n\n$$b = c. D$$\n',
])
def test_leaves_single_sentences_alone(text):
    assert messages(text) == []


def test_fix_breaks_sentences_and_keeps_list_indentation():
    fixed = FixSentences('- First. Second.\n\nThird. Fourth.\n')
    assert fixed == '- First.\n  Second.\n\nThird.\nFourth.\n'
    assert messages(fixed) == []


@pytest.mark.parametrize('text, expected', [
    ('Let $g^*_k$ and $h^*_k$ hold.\n', '\\ast'),
    ('Let $\\{a\\}$ hold.\n', '\\lbrace'),
    ('Let $a\\,b$ hold.\n', 'plain space'),
    ('Let $c = \\$5$ hold.\n', 'amount in words'),
    ('Let $\\operatorname{f}$ hold.\n', '\\operatorname'),
    ('| A |\n|---|\n| $a \\mid b$ $c | d$ |\n', '\\mid'),
    ('Text\n$$a = b$$\nmore text.\n', 'own paragraph'),
    ('Pairs $P : \\mathbb{R}_+ \\times \\mathbb{R}_+ \\to [0,1]$ hold.\n', 'emphasis'),
    ('Let $\\bar{d}_{k}$ and $m_{j}$ hold.\n', 'emphasis'),
    ('A width-$\\delta$ interval.\n', 'opens this math'),
    ('Over x/$y$ here.\n', 'opens this math'),
])
def test_flags_math_github_misreads(text, expected):
    assert any(expected in message for message in messages(text))


@pytest.mark.parametrize('text', [
    'Let $g^{\\ast}_k$ and $h^{\\ast}_k$ hold.\n',
    'Pairs $`P : \\mathbb{R}_+ \\times \\mathbb{R}_+ \\to [0,1]`$ hold.\n',
    'Let $a^{--}_k$ and $(x)_k$ hold.\n',
    'Text.\n\n```math\n\\{a\\} \\\\ g^*_k\n```\n',
    'Text.\n\n$$\\bar{d}_{k} + m_{j}$$\n',
    'Prices of $5 and $6.\n',
    'A ($\\delta$) and **$x$** and the $\\delta$-wide interval.\n',
])
def test_leaves_math_github_reads_correctly(text):
    assert messages(text) == []
