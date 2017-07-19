# goose-article

A customized article class for LaTeX.

## Usage

`goose-article` is a customized class designed for scientific articles. The usage is similar to the default `article`-class while the class takes care of formatting.

By default most of the standard LaTeX-packages are loaded. Any of these packages can be re-loaded, with other defaults, without problems. In addition, the title, the authors and their affiliations, contact information, and optionally a header should be specified.

This results in the following structure:

```latex
%!TEX program = xelatex
\documentclass[...]{goose-article}

% The title (also used as PDF-title).
\title{...}

% The different authors. The optional number(s) correspond to the affiliations.
% N.B. the authblk-package is used, see that package for more information.
\author[1,2]{...}
\author[2]{...}

% The different affiliations. Use "\nl" to enforce a line-break.
\affil[1]{...}
\affil[2]{...}

% The contact information.
\contact{...} % E.g. "\href{mailto:tom@geus.me}{tom@geus.me}"

% The name to put in the PDF-information.
\hypersetup{pdfauthor={...}}

% Text to put in the header. The page number is always used.
\header{...}

% %%%%%%%%%%%%%%
\begin{document}
% %%%%%%%%%%%%%%

\maketitle

\begin{abstract}
...
\end{abstract}

\keywords{...}

...

\bibliography{...}

% %%%%%%%%%%%%
\end{document}
% %%%%%%%%%%%%
```

>   Note that the first line `%!TEX program = xelatex` is only needed if a non-LaTeX-standard font is selected. In fact, only when an editor is used which supports compiler selection this way. For those unfamiliar, XeLaTeX is similar to `pdflatex` but is allows for the usage of TrueType-fonts.

## Options

*   `garamond`, `times`, `verdana`

    Choose a font. The default computer-modern font is used if no font is specified.

*   `narrow`

    Widen the margins of the page, useful during the review process.

*   `doublespacing`

    Set the line-spacing to double, useful during the review process.

*   `namecite`

    Use names instead of numbers to cite to references.

## Citations

Citations and references are handled using [natbib](http://ctan.org/pkg/natbib). In this class, the `unsrtnat` layout is used. Thereby, the extended `unsrtnat.bst` is available that includes output for the `eprint` field. The `goose-article` class creates commands to convert the `doi` and `eprint` fields to links (to `doi.org` and `arxiv.org` respectively).

Following standard natbib, one can use `\cite{...}` or `\citep{...}` for normal citations and `\citet{...}` to include the name. [See also this cheat-sheet](http://merkel.texture.rocks/Latex/natbib.php).

Note that the outputted reference-list depends largely on the content of the included `bib`-file. A simple command-line tool, [bibparse](https://github.com/tdegeus/bibparse), is available to clean-up arbitrary `bib`-files.
