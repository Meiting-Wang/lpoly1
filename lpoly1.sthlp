{smcl}
{right:Created time: Oct 1, 2021}
{* -----------------------------title------------------------------------ *}{...}
{p 0 17 2}
{bf:[W-19] lpoly1} {hline 2} Perform one-dimensional local polynomial regression like {help lpoly}, but more coefficients can be additionally obtained. You can view source code in {browse "https://github.com/Meiting-Wang/lpoly1":github}.


{* -----------------------------Syntax------------------------------------ *}{...}
{title:Syntax}

{p 8 10 2}
{cmd:lpoly1} {it:{help varname:yvar}} {it:{help varname:xvar}}, {opth at:(varname)} [{opth ker:nel(lpoly1##Kernel:kernel)} {opt bw:idth(real>0)} {opt d:egree(integer>=0)} {opth keep:(strings:string)}]


{* -----------------------------Contents------------------------------------ *}{...}
{title:Contents}

{p 4 4 2}
{help lpoly1##Description:Description}{break}
{help lpoly1##Options:Options}{break}
{help lpoly1##Kernel:Kernels}{break}
{help lpoly1##Examples:Examples}{break}
{help lpoly1##Author:Author}


{* -----------------------------Description------------------------------------ *}{...}
{marker Description}{title:Description}

{p 4 4 2}
{cmd:lpoly1} can perform one-dimensional local polynomial regression like {help lpoly}, but derivatives of the fitted value with respect to the independent variable with different orders can be additionally obtained.

{p 4 4 2}
If we set {opt degree} to {bf:p}, {bf:beta0}, {bf:beta1},..., {bf:betap} variables will be generated by default, which represent the {bf:0th} order derivative, {bf:1st} order derivative,..., {bf:pth} order derivative, respectively.

{p 4 4 2}
It is worth noting that this command can be only used in version 16.0 or later.


{* -----------------------------Options------------------------------------ *}{...}
{marker Options}{title:Options}

{synoptset 20}{...}
{synopthdr}
{synoptline}
{synopt :{opth at:(varname)}}Specify the values to calculate derivatives({bf:0th} order to {bf:pth} order).{p_end}
{synopt :{opth ker:nel(lpoly1##Kernel:kernel)}}Specify the kernel function, {bf:gaussian} as the default.{p_end}
{synopt :{opt bw:idth(real>0)}}Specify kernel bandwidth. A positive real number is required, 0.1 as the default.{p_end}
{synopt :{opt d:egree(integer>=0)}}Specify the degree of the polynomial. A nonnegative integer number is required, 3 as the default.{p_end}
{synopt :{opth keep:(strings:string)}}Choose which variables in {bf:beta0}, {bf:beta1},..., {bf:betap} to keep.{p_end}
{synoptline}


{* -----------------------------kernel------------------------------------ *}{...}
{marker Kernel}{title:Kernels}

{synoptset 20}{...}
{synopthdr :kernel}
{synoptline}
{synopt :{opt gaussian}}Gaussian kernel function; the default{p_end}
{synopt :{opt gaussian_m}}Modified gaussian kernel function for test{p_end}
{synopt :{opt epanechnikov}}Epanechnikov kernel function{p_end}
{synopt :{opt epan2}}alternative Epanechnikov kernel function{p_end}
{synopt :{opt biweight}}biweight kernel function{p_end}
{synopt :{opt cosine}}cosine trace kernel function{p_end}
{synopt :{opt parzen}}Parzen kernel function{p_end}
{synopt :{opt rectangle}}rectangle kernel function{p_end}
{synopt :{opt triangle}}triangle kernel function{p_end}
{synoptline}


{* -----------------------------Examples------------------------------------ *}{...}
{marker Examples}{title:Examples}

{p 4 4 2}Do a local polynomial regression of Y on X at x{p_end}
{p 8 10 2}. {bf:lpoly1 Y X, at(x) bwidth(0.1) degree(3) kernel(gaussian)}{p_end}

{p 4 4 2}Use the official command {help lpoly} to do the same to verify the correctness of the result above.{p_end}
{p 8 10 2}. {bf:lpoly Y X, at(x) bwidth(0.1) degree(3) kernel(gaussian) generate(beta0_cheak) nograph}{p_end}


{* -----------------------------Author------------------------------------ *}{...}
{marker Author}{title:Author}

{p 4 4 2}
Meiting Wang{break}
Institute for Economic and Social Research, Jinan University{break}
Guangzhou, China{break}
wangmeiting92@gmail.com
