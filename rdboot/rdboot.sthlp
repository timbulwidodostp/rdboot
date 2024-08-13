{smcl}
{* *! version 1.1.1}{...}
{title:Title}

{phang}
{bf:rdboot} {hline 2} Executes estimation and robust inference for treatment effects in the sharp and fuzzy mean regression discontinuity designs (RDD) based on multiplier bootstrap and bias correction.

    {bf:rdboot} {hline 2} Robust Bootstrap Method for Sharp and Fuzzy Regression Discontinuity Designs.

{marker syntax}{...}
{title:Syntax}

{p 4 17 2}
{cmd:rdboot}
{it:depvar}
{it:runvar}
{ifin}
[{cmd:,} {bf:c}({it:cutoff}) {bf:fuzzy}({it:fuzzyvar})]


{marker description}{...}
{title:Description}

{phang}
{cmd:rdboot} excecutes estimation and robust inference for treatment effects in the sharp and fuzzy mean regression discontinuity designs (RDD) based on the method of multiplier bootstrap and bias correction developed in 
{browse "https://arxiv.org/pdf/1702.04430.pdf":Chiang et al. (2017)}. 
A detailed introduction to this command is given in Yang et al. (2017). 


{marker options}{...}
{title:Options}


{phang}
{bf:c({it:cutoff})} sets the cutoff value used for the RD designs. The default value is {bf: c(0)}.

{phang}
{bf:fuzzy({it:fuzzyvar})} takes the treatment variable used for estimation of the fuzzy case. Receiving no input for this option, the command treats the data as a sharp design by default.



{marker examples}{...}
{title:Examples}

{phang}
({bf:y} outcome variable, {bf:x} running variable, {bf:d} treatment variable)

{phang}Setup

{phang}{cmd:. import excel "https://raw.githubusercontent.com/timbulwidodostp/rdboot/main/rdboot/rdboot.xlsx", sheet("Sheet1") firstrow clear}{p_end}

{phang}Estimation for Sharp RD designs

{phang}{cmd:. rdboot Outcome Running}{p_end}

{phang}Estimation for Fuzzy RD designs

{phang}{cmd:. rdboot Outcome Running, fuzzy(Treatment)}{p_end}

{title:References}

{p 4 8}Chiang, H.D., Y.-C. Hsu, Y. Sasaki, and F. Yang. 2017. A Unified Robust Bootstrap
Method for Sharp/Fuzzy Mean/Quantile Regression Discontinuity/Kink Designs.
{browse "https://arxiv.org/abs/1702.04430v3":arXiv.1702.04430v3}.
{p_end}

{p 4 8}Yang, F., H.D. Chiang, Y.-C. Hsu, and Y. Sasaki. 2017. rdboot: Software Using Robust Bootstrap Method for Sharp and Fuzzy Regression Discontinuity Designs. 

{title:Authors}

{p 4 8}Timbul Widodo, Olah Data Semarang, Indonesian, ID.{p_end}

{p 4 8}Fangzhu Yang, Johns Hopkins University, Baltimore, MD.{p_end}

{p 4 8}Harold. D. Chiang, Vanderbilt University, Nashville, TN.{p_end}

{p 4 8}Yu-Chin Hsu, Academia Sinica, Taipei, Taiwan.{p_end}

{p 4 8}Yuya Sasaki, Vanderbilt University, Nashville, TN.{p_end}



