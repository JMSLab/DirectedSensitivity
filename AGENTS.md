# Conventions for coding agents

The first convention below applies to every file you write; the rest apply to every markdown file in this repository and to every markdown file that [`WORKFLOW.md`](workflow/WORKFLOW.md) asks you to write.
They keep the source easy to diff and the page correct on github.com.

## Ending a file

End every text file you write with a newline, including generated output such as JSON.
Some tools and checks reject a file whose last byte is not a newline.

## Checking a file

Run the checker on each markdown file you write or edit, and fix everything it reports:

```
python -m DirectedSensitivity.MarkdownStyle path/to/file.md
```

With `--fix` it first breaks every line that holds more than one sentence; the other problems need a manual fix.
The checker cannot see every rule below, so follow all of them regardless.

## Lines and paragraphs

- Write one sentence per line.
  GitHub joins consecutive lines into one paragraph, so the page is unchanged and diffs stay small.
- Separate paragraphs with a blank line.
- Do not force a line break with trailing spaces or a trailing backslash.
- Use relative links for files in the repository, so that they resolve on github.com.

## Math

GitHub parses markdown before it renders math, so some characters inside `$` math are consumed or misread.

- Write inline math as `$...$`, and keep it short.
- Put a space or `(` before the `$` that opens inline math; after any other character, such as `-` or `/`, github.com does not render it.
- Put each displayed equation in its own `$$...$$` block on its own lines, with a blank line above and below it.
  Use one relation per block, and do not chain equations with `\qquad`.
- Use a fenced block opened with ```` ```math ```` for display math that needs `\\` line breaks or alignment.
- Write `\ast`, not `*`.
- Write `\lbrace` and `\rbrace`, not `\{` and `\}`.
- Write a plain space, not `\,` or `\;`, and omit `\!`.
- Do not put a dollar sign inside math; write an amount as `$5{,}000$ dollars`.
- A subscript after `}`, `)` or `]`, as in `\bar{d}_{k}` or `\mathbb{R}_+`, can pair with a later `_` as emphasis and break the math.
  Write such a span as `` $`...`$ `` instead of `$...$`; the checker finds these.
- Inside math in a table cell, write `\mid` or `\vert`, not `|`.
- Do not use `\operatorname`.
