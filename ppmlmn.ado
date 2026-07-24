*! ppmlmn 1.0.0  <fill in date>   Yulong Wang
*! PPML gravity estimation with two-way fixed effects and a country-level
*! m-out-of-n bootstrap for heavy-tailed inference.
*! https://sites.google.com/site/yulongwanghome/
*!
*! Reference: Egger, Ji, and Wang, "PPML and Heavy-Tailed Trade Flows:
*! Why Standard Inference Fails and How to Fix It."
*-----------------------------------------------------------------------------

program define ppmlmn, eclass sortpreserve
    version 14.0

    syntax varlist(numeric min=2) [if] [in], ///
        EXporter(varname) IMporter(varname)  ///
        [ Reps(integer 999) Seed(integer -1) G(integer -1) ///
          Level(cilevel) Studentizer(string) noDOTS ]

    marksample touse
    markout `touse' `exporter' `importer'
    gettoken dep indep : varlist

    quietly count if `touse'
    if r(N) == 0 {
        di as err "no observations"
        exit 2000
    }

    if "`studentizer'" == "" local studentizer hc1
    if !inlist("`studentizer'", "hc0", "hc1") {
        di as err "studentizer() must be {bf:hc0} or {bf:hc1}"
        exit 198
    }
    local tau = (100 - `level')/100
    if `reps' < 50 {
        di as err "reps() must be at least 50 for a usable bootstrap"
        exit 198
    }
    local dots = cond("`dots'" == "nodots", 0, 1)
    if `seed' >= 0  set seed `seed'

    * string copies of the FE identifiers for the country (vertex) union
    tempvar estr istr
    capture confirm numeric variable `exporter'
    if !_rc  quietly gen strL `estr' = strofreal(`exporter') if `touse'
    else     quietly gen strL `estr' = `exporter'             if `touse'
    capture confirm numeric variable `importer'
    if !_rc  quietly gen strL `istr' = strofreal(`importer') if `touse'
    else     quietly gen strL `istr' = `importer'             if `touse'

    if `dots' di as txt _n "Bootstrapping (M = `reps'):" _c

    mata: ppmlmn_driver("`dep'", "`indep'", "`estr'", "`istr'", "`touse'", ///
                        `reps', `g', `tau', "`studentizer'", `dots')

    if `dots' di as txt " done."

    * ---- collect results from r() ----
    tempname b V ci0 cib seb nrep
    matrix `b'    = r(b)
    matrix `V'    = r(V)
    matrix `ci0'  = r(ci_hc0)
    matrix `cib'  = r(ci_boot)
    matrix `seb'  = r(se_boot)
    matrix `nrep' = r(nrep)
    local nobs = r(N)
    local Nc   = r(Nc)
    local Gd   = r(G)

    matrix colnames `b'   = `indep'
    matrix colnames `V'   = `indep'
    matrix rownames `V'   = `indep'

    ereturn post `b' `V', esample(`touse') depname(`dep') obs(`nobs')
    ereturn matrix ci_hc0  = `ci0'
    ereturn matrix ci_boot = `cib'
    ereturn matrix se_boot = `seb'
    ereturn matrix nrep    = `nrep'
    ereturn scalar N_countries = `Nc'
    ereturn scalar G     = `Gd'
    ereturn scalar reps  = `reps'
    ereturn scalar level = `level'
    ereturn local  exporter "`exporter'"
    ereturn local  importer "`importer'"
    ereturn local  studentizer "`studentizer'"
    ereturn local  depvar "`dep'"
    ereturn local  title "PPML gravity, country-level m-out-of-n bootstrap"
    ereturn local  cmd "ppmlmn"

    ppmlmn_display, level(`level')
end


*-----------------------------------------------------------------------------
* Display: side-by-side HC0 and bootstrap intervals
*-----------------------------------------------------------------------------
program define ppmlmn_display
    syntax , Level(cilevel)
    local dep   "`e(depvar)'"
    local names : colnames e(b)
    tempname b ci0 cib seb
    matrix `b'   = e(b)
    matrix `ci0' = e(ci_hc0)
    matrix `cib' = e(ci_boot)
    matrix `seb' = e(se_boot)
    local p = colsof(`b')

    di ""
    di as txt "PPML gravity with two-way FE and country-level " _c
    di as txt "{it:m}-out-of-{it:n} bootstrap"
    di as txt "{hline 79}"
    di as txt "Dependent variable: " as res %-12s "`dep'" ///
        _col(48) as txt "Number of obs   = " as res %10.0fc e(N)
    di as txt "Exporter FE: " as res %-8s "`e(exporter)'" as txt " Importer FE: " ///
        as res %-8s "`e(importer)'" _col(48) as txt "Countries (N)   = " as res %10.0f e(N_countries)
    di as txt "Bootstrap reps: " as res %-6.0f e(reps) as txt " Studentizer: " ///
        as res %-4s "`e(studentizer)'" _col(48) as txt "Subsample (G)   = " as res %10.0f e(G)
    di as txt "{hline 79}"
    di as txt _col(28) "HC0 sandwich" _col(56) "{it:m}-out-of-{it:n} bootstrap"
    di as txt %12s abbrev("`dep'",12) " {c |}" ///
        _col(17) "Coef." _col(26) "Std. Err." ///
        _col(38) "[`level'% C.I.]" _col(60) "[`level'% C.I.]"
    di as txt "{hline 13}{c +}{hline 65}"
    forvalues k = 1/`p' {
        local nm : word `k' of `names'
        local cf  = `b'[1,`k']
        local se  = `seb'[1,`k']
        local l0  = `ci0'[`k',1]
        local h0  = `ci0'[`k',2]
        local lb  = `cib'[`k',1]
        local hb  = `cib'[`k',2]
        di as txt %12s abbrev("`nm'",12) " {c |}" ///
            as res _col(15) %8.4f `cf' _col(26) %8.4f `se' ///
            _col(35) %8.4f `l0' _col(45) %8.4f `h0' ///
            _col(56) %9.4f `lb' _col(67) %9.4f `hb'
    }
    di as txt "{hline 13}{c +}{hline 65}"
    di as txt "Note: bootstrap interval is robust to Pareto-tailed shocks " ///
        "(tail index in (1,2));"
    di as txt "      the HC0 sandwich interval is not. See {help ppmlmn}."
end


* ===========================================================================
*  Mata engine
* ===========================================================================
version 14.0
mata:
mata set matastrict on

struct grp {
    real colvector code
    real scalar    K
    real colvector perm
    real matrix    info
}
struct ppmlfit {
    real colvector beta, mu, res
    real matrix    Xt, Hi
    real scalar    n, p, Ko, Kd, ok
    real colvector se_hc0, se_hc1
}

real colvector _codes_str(string colvector v, string colvector keys)
{
    transmorphic   A
    real scalar    i, K
    real colvector c
    K = rows(keys)
    A = asarray_create("string")
    for (i=1; i<=K; i++) asarray(A, keys[i], i)
    c = J(rows(v), 1, .)
    for (i=1; i<=rows(v); i++) c[i] = asarray(A, v[i])
    return(c)
}
struct grp _group(real colvector id)
{
    struct grp     g
    transmorphic   A
    real colvector u, c
    real scalar    i
    u   = uniqrows(id); g.K = rows(u)
    A   = asarray_create("real")
    for (i=1; i<=g.K; i++) asarray(A, u[i], i)
    c = J(rows(id), 1, .)
    for (i=1; i<=rows(id); i++) c[i] = asarray(A, id[i])
    g.code = c; g.perm = order(c, 1); g.info = panelsetup(c[g.perm], 1)
    return(g)
}
real colvector _accum(struct grp g, real colvector v) return(panelsum(v[g.perm], g.info))

real colvector _safesolve(real matrix H, real colvector b, real scalar p)
{
    real matrix R
    real scalar ridge
    R = cholesky(H)
    if (hasmissing(R)) {
        ridge = 1e-9*(trace(H)/p + 1e-12)
        if (hasmissing(luinv(H + ridge*I(p)))) return(J(0,1,.))
        return(luinv(H + ridge*I(p)) * b)
    }
    return(cholsolve(H, b))
}

struct ppmlfit ppml2way(real colvector y, real matrix X, struct grp go, struct grp gd)
{
    struct ppmlfit f
    real scalar    n, p, Ko, Kd, it, s, j, k, mxs, t, d, K, dof
    real colvector Ro, Cd, logRo, logCd, la, lb, beta, w, mu, smo, smd, lb_new, xb, g, step
    real matrix    A, B, Bn, Xt, H, meat, V
    real colvector oc, dc

    f.ok = 1; n = rows(y); p = cols(X)
    Ko = go.K; Kd = gd.K; oc = go.code; dc = gd.code
    Ro = _accum(go, y); Cd = _accum(gd, y)
    if (colmin(Ro) <= 0 | colmin(Cd) <= 0) { f.ok = 0; return(f) }
    logRo = ln(Ro); logCd = ln(Cd)
    la = J(Ko,1,0); lb = J(Kd,1,0); beta = J(p,1,0)
    A = J(Ko,p,0); B = J(Kd,p,0); Xt = X; H = I(p); mu = y

    for (it=1; it<=250; it++) {
        xb = X*beta
        if (hasmissing(xb) | colmax(abs(xb)) > 300) { f.ok = 0; return(f) }
        w = exp(xb)
        for (s=1; s<=60; s++) {
            la = logRo :- ln(_accum(go, exp(lb[dc]):*w))
            lb_new = logCd :- ln(_accum(gd, exp(la[oc]):*w))
            if (colmax(abs(lb_new - lb)) < 1e-12) { lb = lb_new; break }
            lb = lb_new
        }
        mu = exp(la[oc] + lb[dc]) :* w
        smo = _accum(go, mu); smd = _accum(gd, mu)
        for (j=1; j<=100; j++) {
            for (k=1; k<=p; k++) A[,k] = _accum(go, mu:*(X[,k] - B[dc,k])) :/ smo
            Bn = B
            for (k=1; k<=p; k++) Bn[,k] = _accum(gd, mu:*(X[,k] - A[oc,k])) :/ smd
            if (colmax(abs(vec(Bn - B))) < 1e-13) { B = Bn; break }
            B = Bn
        }
        Xt = X - A[oc,] - B[dc,]
        g  = Xt'(y - mu); H = (Xt:*mu)'Xt
        step = _safesolve(H, g, p)
        if (rows(step) == 0) { f.ok = 0; return(f) }
        mxs = colmax(abs(X*step)); if (mxs > 10) step = step*(10/mxs)
        t = 1
        for (d=1; d<=40; d++) { if (colmax(abs(X*(beta+t*step))) <= 250) break; t = t/2 }
        beta = beta + t*step
        if (colmax(abs(t*step)) < 1e-10) break
    }
    f.beta = beta; f.Xt = Xt; f.mu = mu; f.res = y - mu; f.Hi = luinv(H)
    f.n = n; f.p = p; f.Ko = Ko; f.Kd = Kd
    meat = (Xt:*(f.res:^2))'Xt; V = f.Hi*meat*f.Hi
    K = p + (Ko-1) + (Kd-1) + 1; dof = n/(n-K)
    f.se_hc0 = sqrt(diagonal(V)); f.se_hc1 = sqrt(diagonal(V)*dof)
    return(f)
}

real scalar _quantile(real colvector x, real scalar q)
{
    real colvector s
    real scalar    h, lo
    s = sort(select(x, x :< .), 1)
    if (rows(s) == 0) return(.)
    h = (rows(s)-1)*q + 1; lo = floor(h)
    if (lo >= rows(s)) return(s[rows(s)])
    return(s[lo] + (h - lo)*(s[lo+1] - s[lo]))
}
real scalar _fragility(real colvector xk, real scalar n)
{
    real colvector u
    real scalar    i, mn, c
    u = uniqrows(xk); if (rows(u) > 3) return(n)
    mn = n
    for (i=1; i<=rows(u); i++) { c = sum(xk :== u[i]); if (c < mn) mn = c }
    return(mn)
}
real scalar _whichmin(real colvector v)
{
    real scalar i, mi, mv
    mi = 1; mv = v[1]
    for (i=2; i<=rows(v); i++) if (v[i] < mv) { mv = v[i]; mi = i }
    return(mi)
}
real colvector _drop(real colvector v, real scalar i)
{
    if (rows(v) == 1) return(J(0,1,.))
    if (i == 1) return(v[2::rows(v)])
    if (i == rows(v)) return(v[1::rows(v)-1])
    return(v[1::i-1] \ v[i+1::rows(v)])
}

real matrix country_boot(real colvector y, real matrix X, real colvector uo,
    real colvector ud, real scalar Nc, real colvector beta, real colvector se,
    real scalar M, real scalar G, real scalar tau, string scalar stud,
    real colvector nrep, real scalar dots)
{
    real scalar    p, j, ka, k, wi, drop, every, did
    real matrix    T, ci, Xs
    real colvector sel, msk, rows_in, ys, active, frag, bj, sj, bad, lo, hi, pick, perm
    struct grp     go, gd
    struct ppmlfit fj

    p = cols(X); T = J(M, p, .)
    every = max((1, floor(M/50))); did = 0
    for (j=1; j<=M; j++) {
        if (dots & mod(j, every) == 0) { printf("."); displayflush() }
        perm = order(runiform(Nc,1), 1); pick = perm[1::G]
        sel = J(Nc,1,0); sel[pick] = J(G,1,1)
        msk = sel[uo] :& sel[ud]; rows_in = selectindex(msk)
        if (rows(rows_in) < p + 4) continue
        ys = y[rows_in]; Xs = X[rows_in,]
        active = J(0,1,.)
        for (k=1; k<=p; k++) if (colmax(Xs[,k]) > colmin(Xs[,k])) active = active \ k
        if (rows(active) == 0) continue
        frag = J(rows(active),1,.)
        for (ka=1; ka<=rows(active); ka++) frag[ka] = _fragility(Xs[,active[ka]], rows(ys))
        go = _group(uo[rows_in]); gd = _group(ud[rows_in])
        drop = 1
        while (rows(active) > 0 & drop) {
            fj = ppml2way(ys, Xs[,active], go, gd)
            if (fj.ok == 0) {
                wi = _whichmin(frag); active = _drop(active, wi); frag = _drop(frag, wi); continue
            }
            if (stud == "hc1") sj = fj.se_hc1; else sj = fj.se_hc0
            bj = fj.beta; bad = J(0,1,.)
            for (ka=1; ka<=rows(active); ka++)
                if (!(bj[ka] < . & sj[ka] > 0 & sj[ka] < .)) bad = bad \ ka
            if (rows(bad) > 0) {
                wi = bad[_whichmin(frag[bad])]; active = _drop(active, wi); frag = _drop(frag, wi); continue
            }
            for (ka=1; ka<=rows(active); ka++) { k = active[ka]; T[j,k] = (bj[ka]-beta[k])/sj[ka] }
            drop = 0; did = did + 1
        }
    }
    lo = J(p,1,.); hi = J(p,1,.)
    for (k=1; k<=p; k++) {
        nrep[k] = sum(T[,k] :< .)
        if (nrep[k] >= 20) {
            lo[k] = beta[k] - _quantile(T[,k], 1-tau/2)*se[k]
            hi[k] = beta[k] - _quantile(T[,k], tau/2)*se[k]
        }
    }
    ci = (lo, hi); return(ci)
}

void ppmlmn_driver(string scalar dep, string scalar indeps, string scalar evar,
    string scalar ivar, string scalar touse, real scalar M, real scalar Gin,
    real scalar tau, string scalar stud, real scalar dots)
{
    real colvector  y, uo, ud, se, nrep
    real matrix     X, V, ci_hc0, ci_boot
    string colvector e, i, keys, names
    real scalar     Nc, G, p, z, K, dof
    struct grp      go, gd
    struct ppmlfit  f

    y = st_data(., dep, touse); names = tokens(indeps)
    X = st_data(., names, touse)
    e = st_sdata(., evar, touse); i = st_sdata(., ivar, touse)

    keys = uniqrows(e \ i); Nc = rows(keys)
    uo = _codes_str(e, keys); ud = _codes_str(i, keys)

    go = _group(uo); gd = _group(ud)
    f  = ppml2way(y, X, go, gd)
    if (f.ok == 0) {
        errprintf("\nppmlmn: full-sample PPML failed (perfect separation or non-convergence).\n")
        errprintf("        Check for collinear regressors or fully separated countries.\n")
        exit(498)
    }
    p = f.p
    if (stud == "hc1") se = f.se_hc1; else se = f.se_hc0
    V = f.Hi*((f.Xt:*(f.res:^2))'f.Xt)*f.Hi
    if (stud == "hc1") { K = p + (f.Ko-1) + (f.Kd-1) + 1; dof = f.n/(f.n-K); V = V*dof }
    z = invnormal(1 - tau/2)
    ci_hc0 = (f.beta - z*se, f.beta + z*se)

    G = (Gin > 0 ? Gin : floor(Nc^0.7))
    nrep = J(p,1,0)
    ci_boot = country_boot(y, X, uo, ud, Nc, f.beta, se, M, G, tau, stud, nrep, dots)

    st_matrix("r(b)", f.beta'); st_matrix("r(V)", V)
    st_matrix("r(ci_hc0)", ci_hc0); st_matrix("r(ci_boot)", ci_boot)
    st_matrix("r(se_boot)", se'); st_matrix("r(nrep)", nrep')
    st_numscalar("r(N)", f.n); st_numscalar("r(Nc)", Nc); st_numscalar("r(G)", G)
}
end
