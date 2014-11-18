\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{template}

\LoadClass[a4paper,bibliography=totoc,BCOR=5mm]{scrbook}

%% auxiliary packages
\usepackage[english, ngerman]{babel}
\usepackage{bold-extra}
\usepackage{cite}
%\usepackage{amsthm}
%\usepackage{lscape}

\usepackage[toc, page]{appendix}

\usepackage{hyperref} %insertes hyperref to mentioned sections

\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}

\usepackage{framed}
\usepackage[framed,thmmarks]{ntheorem}

\usepackage{amsmath}
\usepackage{amssymb}
%\usepackage{a4wide}
\usepackage{relsize}
\usepackage{tikz} %for generating images

\usepackage{todonotes}

\usepackage{bigdelim}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{ltxtable}
\usepackage{float}
\usepackage{pgfplotstable, booktabs} %to generate tables from data files

\usepackage{color}
\usepackage{stmaryrd} % for llbracket

\usepackage[chapter]{algorithm}

\usepackage{enumerate} %to allow different kind of enumerations
 \usepackage{algpseudocode}
 
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
\newcounter{thmc}[chapter] %theorem counter
\renewcommand{\thethmc}{\thechapter.\arabic{thmc}}


\newtheorem{theorem}[thmc]{Theorem}
\newtheorem{proposition}[thmc]{Proposition}
\newtheorem{example}[thmc]{Example}
\newtheorem{lemma}[thmc]{Lemma}
\newtheorem{definition}[thmc]{Definition}
\newtheorem{claim}[thmc]{claim}
\newtheorem{corollary}[thmc]{Corollary}
\newtheorem{test}[thmc]{Test}

\theoremstyle{nonumberplain}
\theoremstyle{nonumberbreak}
\theoremsymbol{\ensuremath{\square}}

\newtheorem{proof}{Proof}

% -----
%---define pythonplotstyle

\usepackage[procnames]{listings} %to include python code
\definecolor{keywords}{RGB}{255,0,90}
\definecolor{comments}{RGB}{0,0,113}
\definecolor{red}{RGB}{160,0,0}
\definecolor{green}{RGB}{0,150,0}
	
\lstset{language=Python, 
        basicstyle=\ttfamily\small, 
        keywordstyle=\color{keywords},
        commentstyle=\color{comments},
        stringstyle=\color{red},
        showstringspaces=false,
        breaklines     = true,
%        identifierstyle=\color{green},
        procnamekeys={def,class}}

% -----
%---define pgfplottablestyle

\pgfkeys{
	/pgfplots/table/string type in dec sep align/.style={
		string type,
		postproc cell content/.code={%
			\ifnum\pgfplotstablepartno=0%
			\pgfkeys{/pgfplots/table/@cell content/.add={}{&}}
			\fi
		}%
	}
}

\pgfplotsset{compat=1.8}
\pgfplotstableset{
    columns/iterations/.style={column name={$h$},% Set the name to be used for the first column
    	/pgf/number format/precision=0,
    	preproc/expr={2^##1}, %preposs (calculate 2^i) 
		postproc cell content/.append style={
			/pgfplots/table/@cell content/.add={1/{}}{}, %add 1/ before every number
		},
    },
    %make sure, dec sep align does not affect the iterations column
    columns/iterations/in dec sep align/.style={
		postproc cell content/.code={%
			\ifnum\pgfplotstablepartno=0%
			\pgfkeys{/pgfplots/table/@cell content/.add={}{&}}
			\fi
		}%
	},
	columns/l2error/.style={
        column name=$L^2$ error,  % ... and the second
		sci,sci zerofill,sci sep align,
%		dec sep align,      % align on the decimal marker
        /pgf/number format/sci e, 
        /pgf/number format/sci precision=4 ,    % print 4 digits
        /pgf/number format/precision=4     % print 4 digits
            },
    columns/h1error/.style={
		column name=$H^1$ error,  % ... and the second
		sci,sci zerofill,sci sep align,
        /pgf/number format/sci e, 
        /pgf/number format/fixed zerofill,  % print trailing zeros
        /pgf/number format/sci precision=4,     % print 4 digits
        /pgf/number format/precision=4     % print 4 digits
            },
    every row 0 column N/.style={set content={-}},
    every head row/.style={
        before row=\toprule,    % booktabs rules
        after row=\midrule
    },
    every last row/.style={
        after row=\bottomrule
    }} %{../../FEniCS/data/MA1_Neilan_deg22_l2errornorm}
