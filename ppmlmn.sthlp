{smcl}
{* version 1.1.0  <fill in date>}{...}
{vieweralsosee "[R] poisson" "help poisson"}{...}
{viewerjumpto "Syntax" "ppmlmn##syntax"}{...}
{viewerjumpto "Description" "ppmlmn##description"}{...}
{viewerjumpto "Options" "ppmlmn##options"}{...}
{viewerjumpto "Examples" "ppmlmn##examples"}{...}
{viewerjumpto "Stored results" "ppmlmn##results"}{...}
{viewerjumpto "References" "ppmlmn##references"}{...}
{viewerjumpto "Author" "ppmlmn##author"}{...}
{title:Title}

{p 4 8 2}{cmd:ppmlmn} {hline 2} PPML gravity estimation with two-way fixed
effects and a country-level {it:m}-out-of-{it:n} bootstrap for heavy-tailed
inference{p_end}


{marker syntax}{title:Syntax}

{p 8 16 2}
{cmd:ppmlmn} {depvar} {indepvars} {ifin}{cmd:,}
{cmdab:ex:porter(}{varname}{cmd:)}
{cmdab:im:porter(}{varname}{cmd:)}
[{it:options}]{p_end}

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt :{cmdab:ex:porter(}{varname}{cmd:)}}origin (exporter) identifier{p_end}
{synopt :{cmdab:im:porter(}{varname}{cmd:)}}destination (importer) identifier{p_end}

{syntab:Bootstrap}
{synopt :{cmdab:r:eps(}{it:#}{cmd:)}}number of bootstrap subsamples {it:M}; default {cmd:reps(999)}{p_end}
{synopt :{cmd:g(}{it:#}{cmd:)}}countries drawn per subsample; default {cmd:floor(N^0.7)}{p_end}
{synopt :{cmdab:s:eed(}{it:#}{cmd:)}}random-number seed; if omitted the current seed is used{p_end}
{synopt :{cmd:studentizer(}{it:str}{cmd:)}}{cmd:hc0} or {cmd:hc1}; default {cmd:hc1}{p_end}
{synopt :{cmd:nodots}}suppress the bootstrap progress dots{p_end}

{syntab:Reporting}
{synopt :{cmdab:l:evel(}{it:#}{cmd:)}}confidence level; default {cmd:level(95)}{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{title:Description}

{pstd}
{cmd:ppmlmn} estimates the Poisson pseudo-maximum-likelihood (PPML) gravity
model

{p 8 8 2}E[{it:y_od}] = exp({it:a_o} + {it:b_d} + {it:x_od}'{bf:b}){p_end}

{pstd}
with exporter effects {it:a_o} and importer effects {it:b_d}, and reports two
95% confidence intervals for each slope coefficient: the conventional HC0
heteroskedasticity-robust ("sandwich") interval, and a country-level
{it:m}-out-of-{it:n} bootstrap interval.{p_end}

{pstd}
When bilateral shocks have a heavy (Pareto) tail with index {it:alpha} in (1,2),
the PPML score has infinite variance and the sandwich interval under-covers: its
coverage falls well below the nominal level. The bootstrap interval remains valid
in this regime. It draws {it:G} countries without replacement, keeps the induced
subnetwork of dyads, refits the two-way fixed-effect PPML on that subnetwork, and
self-normalizes each subsample slope by its own sandwich standard error,
{it:T_j} = ({it:b_j} - {it:b})/{it:se_j}. The interval is symmetric about the
coefficient, {it:b} -/+ {it:q} {it:se}, where {it:q} is the {cmd:level()}
quantile of the absolute statistics |{it:T_j}|.{p_end}

{pstd}
Fixed effects are profiled out by Sinkhorn (matrix balancing); slopes are
updated by Newton-Raphson. Countries that are perfectly separated (a zero
exporter or importer marginal) are dropped, as in standard PPML. No additional
packages are required.{p_end}


{marker options}{title:Options}

{phang}{cmd:exporter(}{varname}{cmd:)} and {cmd:importer(}{varname}{cmd:)}
specify the origin and destination identifiers (string or numeric). They define
the two-way fixed effects and the country (vertex) resampling unit. Both are
required.{p_end}

{phang}{cmd:reps(}{it:#}{cmd:)} sets the number of bootstrap subsamples {it:M}
(default 999). Larger {it:M} gives more stable tail quantiles. The accompanying
paper uses {cmd:reps(2000)}.{p_end}

{phang}{cmd:g(}{it:#}{cmd:)} sets the number of countries drawn per subsample.
The default {cmd:floor(N^0.7)} satisfies the {it:m}-out-of-{it:n} rate condition
(the subsample grows slower than {it:N}); for {it:N} = 79 this is 21.{p_end}

{phang}{cmd:seed(}{it:#}{cmd:)} sets the random-number seed for reproducibility.
If omitted, the current RNG state is used.{p_end}

{phang}{cmd:studentizer(hc0|hc1)} selects the standard error used to
self-normalize (default {cmd:hc1}, which applies the small-sample {it:n}/({it:n}-{it:k})
factor and matches the "HC0" column reported in many gravity packages).{p_end}

{phang}{cmd:nodots} suppresses the progress dots printed during the bootstrap.{p_end}

{phang}{cmd:level(}{it:#}{cmd:)} sets the confidence level (default 95).{p_end}


{marker examples}{title:Examples}

{pstd}Aggregate bilateral trade, two-way FE, default bootstrap:{p_end}
{p 8 12 2}{cmd:. ppmlmn trade ln_dist contig comlang_off colony, exporter(iso_o) importer(iso_d)}{p_end}

{pstd}A reproducible run with 2,000 draws and a fixed seed:{p_end}
{p 8 12 2}{cmd:. ppmlmn trade ln_dist contig comlang_off colony, ex(iso_o) im(iso_d) reps(2000) seed(42)}{p_end}

{pstd}Migration flows with a user-set subsample size, quietly:{p_end}
{p 8 12 2}{cmd:. ppmlmn mig ln_dist contig comlang_off colony, ex(orig) im(dest) g(40) nodots}{p_end}

{pstd}After estimation, the HC0 variance is in {cmd:e(V)}, so standard
postestimation works:{p_end}
{p 8 12 2}{cmd:. test contig comlang_off}{p_end}
{p 8 12 2}{cmd:. matrix list e(ci_boot)}{p_end}


{marker results}{title:Stored results}

{pstd}{cmd:ppmlmn} stores the following in {cmd:e()}:{p_end}

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt :{cmd:e(N)}}number of observations (dyads used){p_end}
{synopt :{cmd:e(N_countries)}}number of countries {it:N}{p_end}
{synopt :{cmd:e(G)}}subsample size (countries per draw){p_end}
{synopt :{cmd:e(reps)}}bootstrap replications {it:M}{p_end}
{synopt :{cmd:e(level)}}confidence level{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt :{cmd:e(cmd)}}{cmd:ppmlmn}{p_end}
{synopt :{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt :{cmd:e(exporter)}}{cmd:e(importer)} fixed-effect identifiers{p_end}
{synopt :{cmd:e(studentizer)}}{cmd:hc0} or {cmd:hc1}{p_end}

{p2col 5 22 26 2: Matrices}{p_end}
{synopt :{cmd:e(b)}}coefficient vector{p_end}
{synopt :{cmd:e(V)}}HC0 sandwich variance{p_end}
{synopt :{cmd:e(ci_hc0)}}HC0 confidence interval ({it:k} x 2){p_end}
{synopt :{cmd:e(ci_boot)}}symmetric bootstrap confidence interval ({it:k} x 2){p_end}
{synopt :{cmd:e(se_boot)}}studentizer used by the bootstrap{p_end}
{synopt :{cmd:e(nrep)}}valid bootstrap draws per coefficient{p_end}


{marker references}{title:References}

{phang}Egger, P., T. Ji, and Y. Wang. "PPML and Heavy-Tailed Trade Flows: Why
Standard Inference Fails and How to Fix It." Working paper.{p_end}


{marker author}{title:Author}

{pstd}Yulong Wang{p_end}
{pstd}{browse "https://sites.google.com/site/yulongwanghome/"}{p_end}
