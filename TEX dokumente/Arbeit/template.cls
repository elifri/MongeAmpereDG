\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{template}

\LoadClass[a4paper,bibliography=totoc,BCOR=5mm]{scrbook}

%% auxiliary packages
\usepackage[english]{babel}
\usepackage{bold-extra}
\usepackage{cite}
%\usepackage{amsthm}
%\usepackage{lscape}

\usepackage{hyperref}

\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}

\usepackage{framed}
\usepackage[framed,thmmarks]{ntheorem}

\usepackage{amsmath}
\usepackage{amssymb}
%\usepackage{a4wide}
\usepackage{relsize}
\usepackage{tikz}
\usepackage{todonotes}
\usepackage{multirow}
\usepackage{bigdelim}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{ltxtable}
\usepackage{float}

\usepackage{color}
\usepackage{stmaryrd} % for llbracket


%% page layout
\setlength{\parindent}{0pt}
% let latex do its best...
\setlength{\parskip}{1ex plus 0.25ex minus 0.25ex}
\setlength{\footskip}{20pt}

%\addtolength{\textheight}{0.75in}

% but please dont do this...
\widowpenalty=2000
\clubpenalty=2000
% -----

%% tikz settings

\definecolor{rwthblue}{HTML}{0066CC}
\usetikzlibrary{calc}

% drawing hypergraphs
\tikzstyle{node} = [draw, circle]
\tikzstyle{edge} = [draw, rectangle]
% -----


%% theorems
\theoremstyle{break}
\theorembodyfont{\rmfamily}
\newcounter{thmc}[section] %theorem counter
\renewcommand{\thethmc}{\thechapter.\arabic{thmc}}


\newtheorem{theorem}[thmc]{Theorem}
\newtheorem{example}[thmc]{Example}
\newtheorem{lemma}[thmc]{Lemma}
\newtheorem{definition}[thmc]{Definition}
\newtheorem{claim}[thmc]{claim}
\newtheorem{corollary}[thmc]{Corollary}

\theoremstyle{nonumberplain}
\theoremstyle{nonumberbreak}
\theoremsymbol{\ensuremath{\square}}

\newtheorem{proof}{Proof}

% -----
