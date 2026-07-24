# ppmlmn — installation

`ppmlmn` estimates two-way fixed-effect PPML gravity models and reports a
country-level *m*-out-of-*n* bootstrap confidence interval that stays valid when
bilateral shocks are heavy-tailed (Pareto, tail index in (1,2)), a case in which
the usual HC0/sandwich interval under-covers.

Files in this package:

```
ppmlmn.ado      the command (self-contained; no other packages needed)
ppmlmn.sthlp    help file        ( . help ppmlmn )
ppmlmn.pkg      net-install manifest
stata.toc       net-install table of contents
```

---

## For users — installing

### Option A — automatic install (recommended)

Once the four files are hosted in a single web folder (see "For the maintainer"
below), users install with one line:

```stata
. net install ppmlmn, from("https://raw.githubusercontent.com/yulongwang06/Tail-PPML/main/") replace
```

or interactively:

```stata
. net from "https://raw.githubusercontent.com/yulongwang06/Tail-PPML/main/"
. net install ppmlmn, replace
```

Then `. help ppmlmn` and you are ready.

### Option B — manual install (from the ZIP on the website)

1. Download `ppmlmn.zip` and unzip it to a folder, e.g. `C:\ppmlmn` (or `~/ppmlmn`).
2. In Stata:

```stata
. net install ppmlmn, from("C:/ppmlmn") replace
```

(Forward slashes work on every platform.) Alternatively, copy `ppmlmn.ado` and
`ppmlmn.sthlp` into your `PERSONAL` ado-directory — run `. personal dir` to find
it — and they will be available in every session.

### Quick check

```stata
. which ppmlmn
. help ppmlmn
```

---


## Citation

Egger, P., T. Ji, and Y. Wang, "PPML and Heavy-Tailed Trade Flows: Why Standard
Inference Fails and How to Fix It."

Author / support: Yulong Wang — https://sites.google.com/site/yulongwanghome/
