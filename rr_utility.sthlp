{smcl}
help for {hi:rr_utility} version 0.1.1 (Date: 25.04.2024). {right: (Caspar Kaiser & Anthony Lepinteur)}
{hline}
{title:Utility to check for coefficient reversals in regressions of ordinal dependent variables}

{title:Syntax}

{p 4 4 2}
{cmd:rr_utility}
[{cmd:,} {it:options}]

{synoptset tabbed}{...}
{synopthdr}
{synoptline}
{syntab:General {help rr_utility##opt_general:[+]}}
{synopt:{opt py:thon}}Finds least non-linear transformation for each coefficient (where possible). Requires that Stata can call Python and that all dependencies are installed. {p_end}
{synopt:{opt ga:mma}}Specifies a shift by gamma (see below for more details). Causes additional output to be displayed. {p_end}
{synopt:{opt fast}}Runs much faster but does not compute critical values for p-values. {p_end}
{synopt:{cmd:start(}{it:real}{cmd:)}}Specifies the smallest value of c over which should be searched. Default is -2. {p_end}
{synopt:{cmd:end(}{it:real}{cmd:)}}Specifies the largest value of c over which should be searched. Default is 2. {p_end}
{synopt:{cmd:prec:ision(}{it:real}{cmd:)}}Specifies the 'density' or 'precision' of the grid of c. For example prec(0.1) says that we evaluate values of c in steps of 0.1. Default is 0.1. {p_end}
{synopt:{cmd:range(}{it:{help numlist}}{cmd:)}}Specifies the range we want the scale to be on. Two numbers must be supplied, where the first number is the desired minimum of the scale and the second number is desire maximum of the scale. Reversal cost and c-values are normally unaffected unless the gamma option is used. Default is the range of the original scale {p_end}
{synopt:{cmd:critval(}{it:real}{cmd:)}}Specifies the alpha level where we speak of statistical significance. Default is 0.05.{p_end}

{syntab:More options {help rr_utility##opt_more:[+]}}
{synopt:{cmd:dstub(}{it:str}{cmd:)}}Specifies that the binary dummies should be saved and storted in a stub specified by string. {p_end}
{synopt:{cmd:keep(}{it:{help varlist}}{cmd:)}}Specifies list of variables to be kept in the displayed results table(s). Default is all variables in the original regression. {p_end}

{synoptline}
{p 4 4} Things to say. {p_end}
{p 4 4} More things to say. {p_end}

{marker introduction}{...}
{title:Introduction}

{p 4 4}{cmd:rr_utility} is a command to do things. More things to say. {p_end}
{p 4 4}At minimum, {cmd:rr_utility} requires Stata version 17.0. {p_end}

{marker description}{...}
{title:Description}

{p 4 4}Even more things to say.

{marker description_eq1}{...}
{p 6} Maybe an equation? y=bx+e	s.t. b<0. (1).	

{p 4 4} And more text.  


{title:Options}
{marker opt_general}{...}
{dlgtab:General}

{p 4 4} {opt py:thon} Finds least non-linear transformation for each coefficient (where possible). Requires that Stata can call Python and that all dependencies are installed. I can reference equations. See {help rr_utility##description_eq1:(1)}.

{marker opt_more}{...}
{dlgtab:More options}
{p 4 4} I also can reference other parts of the help file. See {help rr_utility##introduction:introduction}.


