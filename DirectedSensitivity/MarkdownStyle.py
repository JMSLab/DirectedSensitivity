"""Checks markdown for the conventions in AGENTS.md that a script can see."""
import re
import sys
import unicodedata
from pathlib import Path

ABBREVIATIONS = {'al', 'cf', 'etc', 'vs', 'Dr', 'Eq', 'Fig', 'Jr', 'Mr', 'Ms', 'No', 'St'}
MATH_ESCAPES = {r'\{': r'\lbrace', r'\}': r'\rbrace', r'\,': 'a plain space',
                r'\;': 'a plain space', r'\!': 'nothing', '\\\\': 'a ```math block',
                r'\$': 'the amount in words'}
UNSUPPORTED_MACROS = {r'\operatorname'}
FENCE = re.compile(r'^\s*(`{3,}|~{3,})')
CODE = re.compile(r'(?P<math>\$(`+).+?\2\$)|(`+).+?\3', re.S)
DISPLAY_MATH = re.compile(r'\$\$(.+?)\$\$', re.S)
INLINE_MATH = re.compile(r'(?<![\\$])\$(?![\s$])((?:\\.|[^$\\])+?)(?<!\s)\$(?!\d)', re.S)
LINK_TARGET = re.compile(r'\]\([^)\s]*\)')
LINE_PREFIX = re.compile(r'^\s*(?:>\s*)*(?:#{1,6}\s+)?(?:(?:\d+[.)]|[-*+])\s+)?(?:\[[ xX]\]\s+)?'
                         r'(?:\*\*[^*]+\*\*[.:]?\s*)?')
LIST_PREFIX = re.compile(r'^\s*(?:>\s*)*(?:(?:\d+[.)]|[-*+])\s+)?')
SENTENCE_END = re.compile(r'([.?!])["\')\]*_]*([ \t]+)(?=[A-Z§\[*_])')


def Main(arguments):
    fix = '--fix' in arguments
    paths = [argument for argument in arguments if argument != '--fix']
    if fix:
        for path in paths:
            Path(path).write_text(FixSentences(Path(path).read_text()))
    problems = [f'{path}:{line}: {message}' for path in paths
                for line, message in CheckMarkdown(Path(path).read_text())]
    print('\n'.join(problems) if problems else f'{len(paths)} file(s) follow AGENTS.md.')
    return 1 if problems else 0


def CheckMarkdown(text):
    """Return (line number, problem) pairs for one markdown document."""
    problems = []
    for first, lines in _Paragraphs(text.splitlines()):
        paragraph = _Paragraph('\n'.join(lines))
        problems += [(first + 1 + line, message) for line, message in paragraph.Problems()]
    return sorted(problems)


def FixSentences(text):
    """Break every line that holds more than one sentence, one sentence per line."""
    lines = text.splitlines()
    for first, block in reversed(list(_Paragraphs(lines))):
        lines[first:first + len(block)] = _Paragraph('\n'.join(block)).Broken().split('\n')
    return '\n'.join(lines) + ('\n' if text.endswith('\n') else '')


class _Paragraph:
    """One run of non-blank lines, with its code and math spans located."""

    def __init__(self, text):
        self.text = text
        self.spans = {'code': [], 'math': []}
        self.openings = []
        masked = list(text)
        for kind, pattern in (('code', CODE), ('math', DISPLAY_MATH), ('math', INLINE_MATH)):
            for match in pattern.finditer(''.join(masked)):
                self.spans[kind].append(match.span())
                if pattern is INLINE_MATH or (pattern is CODE and match.group('math')):
                    self.openings.append(match.start())
                for index in range(*match.span()):
                    if text[index] != '\n':
                        masked[index] = 'c' if kind == 'code' else '§'
        self.masked = ''.join(masked)
        displays = [span for span in self.spans['math'] if text.startswith('$$', span[0])]
        self.whole_display = len(displays) == 1 and text[slice(*displays[0])] == text.strip()
        self.wedged_display = displays[0][0] if displays and not self.whole_display else None
        self.table = text.lstrip().startswith('|')

    def Problems(self):
        """Return (line within the paragraph, problem) pairs."""
        problems = []
        if self.wedged_display is not None:
            problems.append((self._LineOf(self.wedged_display), 'put a $$ block in its own '
                             'paragraph, with a blank line above and below'))
        for begin, end in self.spans['math']:
            problems += [(self._LineOf(begin), message)
                         for message in _MathProblems(self.text[begin:end], self.table)]
        problems += [(self._LineOf(begin), 'put a space or ( before the $ that opens this math')
                     for begin in self.openings if not _OpensMath(self.text, begin)]
        problems += [(self._LineOf(index), 'github.com reads this _ in math as emphasis; '
                      'write the span as $`...`$') for index in self._EmphasisInMath()]
        problems += [(self._LineOf(begin), 'one sentence per line: break after '
                      f'"...{self.text[max(begin - 30, 0):begin].lstrip()}"')
                     for begin, _ in self.SentenceBreaks()]
        return problems

    def SentenceBreaks(self):
        """Return the (begin, end) whitespace runs that separate two sentences on one line."""
        if self.whole_display or self.table:
            return []
        breaks, offset = [], 0
        for line in self.masked.split('\n'):
            breaks += [(offset + begin, offset + end) for begin, end in _SentenceGaps(line)]
            offset += len(line) + 1
        return breaks

    def Broken(self):
        text = self.text
        for begin, end in reversed(self.SentenceBreaks()):
            line = text[text.rfind('\n', 0, begin) + 1:]
            indent = re.sub(r'[^>\s]', ' ', LIST_PREFIX.match(line).group())
            text = text[:begin] + '\n' + indent + text[end:]
        return text

    def _EmphasisInMath(self):
        """Openers of underscore pairs that github.com would read as emphasis, touching math."""
        if self.whole_display:
            return []
        inside_code = _Membership(len(self.text), self.spans['code'])
        inside_math = _Membership(len(self.text), self.spans['math'])
        openers, flagged = [], []
        for index, char in enumerate(self.text):
            if char != '_' or inside_code[index] or (index and self.text[index - 1] == '\\'):
                continue
            left, right, before, after = _Flanking(self.text, index)
            if right and (not left or _IsPunctuation(after)) and openers:
                opener = openers.pop()
                if inside_math[opener] or inside_math[index]:
                    flagged.append(opener)
            elif left and (not right or _IsPunctuation(before)):
                openers.append(index)
        return flagged

    def _LineOf(self, index):
        return self.text.count('\n', 0, index)


def _Paragraphs(lines):
    """Yield (first line index, lines) for each run of non-blank lines outside fences."""
    fence, first, block = None, None, []
    for index, line in enumerate(lines):
        match = FENCE.match(line)
        if fence:
            if match and match.group(1)[0] == fence[0] and len(match.group(1)) >= len(fence):
                fence = None
            continue
        if match or not line.strip():
            if block:
                yield first, block
            fence = match.group(1) if match else None
            first, block = None, []
            continue
        first = index if first is None else first
        block.append(line)
    if block:
        yield first, block


def _MathProblems(math, table):
    messages = [r'write \ast, not *, inside math'] if '*' in math else []
    commands = set(re.findall(r'\\(?:[A-Za-z]+|.)', math))
    messages += [f'write {MATH_ESCAPES[command]}, not {command}, inside $ math'
                 for command in sorted(commands & MATH_ESCAPES.keys())]
    messages += [f'do not use {command}' for command in sorted(commands & UNSUPPORTED_MACROS)]
    if table and '|' in math:
        messages.append(r'write \mid or \vert, not |, inside math in a table')
    return messages


def _SentenceGaps(line):
    body = LINK_TARGET.sub(lambda match: ']' + 'u' * (len(match.group()) - 1), line)
    if body.lstrip().startswith(('|', '<')):
        return []
    skip = len(LINE_PREFIX.match(body).group())
    gaps = []
    for match in SENTENCE_END.finditer(body, skip):
        word = re.search(r'(\S*)$', body[skip:match.start(1)]).group(1)
        initial = re.fullmatch(r'\(?[A-Z]', word) or re.fullmatch(r'(?:[A-Za-z]\.)+[A-Za-z]', word)
        if not initial and word.strip('(') not in ABBREVIATIONS:
            gaps.append(match.span(2))
    return gaps


def _OpensMath(text, index):
    """Whether github.com starts inline math at the $ at ``index``."""
    before = text[:index]
    return not before or before[-1].isspace() or before[-1] == '(' or before.endswith('**')


def _Flanking(text, index):
    """CommonMark's left- and right-flanking tests for the delimiter at ``index``."""
    before = text[index - 1] if index else '\n'
    after = text[index + 1] if index + 1 < len(text) else '\n'
    left = not after.isspace() and (not _IsPunctuation(after)
                                    or before.isspace() or _IsPunctuation(before))
    right = not before.isspace() and (not _IsPunctuation(before)
                                      or after.isspace() or _IsPunctuation(after))
    return left, right, before, after


def _Membership(length, spans):
    member = [False] * length
    for begin, end in spans:
        member[begin:end] = [True] * (end - begin)
    return member


def _IsPunctuation(char):
    return unicodedata.category(char)[0] in 'PS'


if __name__ == '__main__':
    sys.exit(Main(sys.argv[1:]))
