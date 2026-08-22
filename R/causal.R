# The causal layer: dagitty does the graph theory (adjustment sets,
# d-separation, layout); AlphaSDM does the measurement (embedding-derived
# nodes carrying reliabilities) and the reporting. The DAG is the user's
# stated hypothesis — typically drawn visually at dagitty.net and pasted
# here as model code, coordinates included.

#' Assemble a causal model from a DAG and measured node data
#'
#' @param dag A dagitty object or dagitty model string (e.g. pasted from
#'   the dagitty.net visual editor; coordinates are kept for reporting).
#' @param data Data frame with one column per measured DAG node. Nodes
#'   without columns are treated as latent (fine for graph math,
#'   excluded from testing/estimation).
#' @param measurements Optional named list, one entry per node, each a
#'   list with `reliability` (0-1, e.g. a held-out R^2) and `source`
#'   (free text, e.g. "embedding index vs SRTM"). Nodes from the user's
#'   own data default to reliability 1.
#' @return A causal_model object.
#' @export
causal_model <- function(dag, data, measurements = list()) {
  if (!requireNamespace("dagitty", quietly = TRUE))
    stop("The causal layer needs the dagitty package: install.packages(\"dagitty\")",
         call. = FALSE)
  if (is.character(dag)) dag <- dagitty::dagitty(dag)
  nodes <- names(dagitty::coordinates(dag)$x)
  if (length(nodes) == 0) nodes <- setdiff(unique(unlist(dagitty::edges(dag)[, 1:2])), NA)
  measured <- intersect(nodes, names(data))
  latent <- setdiff(nodes, measured)
  meas <- lapply(stats::setNames(measured, measured), function(n) {
    m <- measurements[[n]]
    list(reliability = if (!is.null(m$reliability)) m$reliability else 1,
         source = if (!is.null(m$source)) m$source else "user data")
  })
  structure(list(dag = dag, data = data, nodes = nodes,
                 measured = measured, latent = latent,
                 measurements = meas),
            class = "causal_model")
}

#' Valid adjustment sets, ranked by measurement reliability
#'
#' All minimal sufficient adjustment sets from the backdoor criterion
#' (dagitty), ranked by the reliability of their members — per Arif &
#' MacNeil (2023): when several sets are valid, use the one measured
#' most accurately. Sets containing latent nodes are listed but marked
#' unusable.
#' @export
causal_adjustment_sets <- function(cm, exposure, outcome,
                                   effect = c("total", "direct")) {
  effect <- match.arg(effect)
  sets <- dagitty::adjustmentSets(cm$dag, exposure = exposure,
                                  outcome = outcome, effect = effect,
                                  type = "minimal")
  if (length(sets) == 0)
    return(data.frame(set = "(none needed)", usable = TRUE,
                      min_reliability = 1, mean_reliability = 1))
  out <- do.call(rbind, lapply(sets, function(s) {
    s <- as.character(s)
    rel <- vapply(s, function(n)
      if (n %in% cm$measured) cm$measurements[[n]]$reliability else NA_real_,
      numeric(1))
    data.frame(set = if (length(s)) paste(s, collapse = " + ") else "{}",
               usable = !anyNA(rel),
               min_reliability = if (length(rel)) min(rel) else 1,
               mean_reliability = if (length(rel)) mean(rel) else 1)
  }))
  out[order(-out$usable, -out$min_reliability, -out$mean_reliability), ]
}

#' Test the DAG against the measured data
#'
#' Every conditional independence the DAG implies among measured nodes,
#' tested with dagitty::localTests. Estimates with confidence intervals
#' are returned as they are; an interval excluding zero contradicts the
#' DAG (or the data/measurements) and is flagged.
#' @export
causal_dag_test <- function(cm, conf = 0.95) {
  d <- cm$data[, cm$measured, drop = FALSE]
  d <- d[stats::complete.cases(d), ]
  lt <- tryCatch(
    dagitty::localTests(cm$dag, data = d, type = "cis", conf.level = conf),
    error = function(e) NULL)
  if (is.null(lt) || nrow(lt) == 0)
    return(data.frame(note = "no testable implications among measured nodes"))
  lt$contradicts_dag <- lt$`2.5%` > 0 | lt$`97.5%` < 0
  lt
}

#' Estimate a causal effect with the best usable adjustment set
#'
#' Fits `outcome ~ exposure + adjusters (+ controls)` with the family
#' given, using the top-ranked usable adjustment set. `controls` is for
#' the detection/effort subgraph only (observer expertise, duration,
#' protocol...) — nodes that affect observation, never the ecology; they
#' are kept out of the DAG's causal claims. Reports the effect with the
#' adjusters' reliabilities alongside, because imperfectly measured
#' adjusters leave residual confounding (attenuation is visible, not
#' hidden).
#' @export
causal_effect <- function(cm, exposure, outcome,
                          effect = c("total", "direct"),
                          family = stats::gaussian(), controls = NULL,
                          engine = c("glm", "gam", "bam")) {
  effect <- match.arg(effect); engine <- match.arg(engine)
  sets <- causal_adjustment_sets(cm, exposure, outcome, effect)
  use <- sets[sets$usable, , drop = FALSE]
  if (nrow(use) == 0)
    stop("No usable adjustment set: every valid set contains an unmeasured ",
         "node (", paste(sets$set, collapse = " | "), "). Consider the ",
         "frontdoor route through a measured mediator.", call. = FALSE)
  adj <- if (use$set[1] %in% c("{}", "(none needed)")) character(0)
         else strsplit(use$set[1], " \\+ ")[[1]]
  rhs <- c(exposure, adj, controls)
  f <- stats::reformulate(rhs, response = outcome)
  d <- cm$data
  if (engine %in% c("gam", "bam") && !requireNamespace("mgcv", quietly = TRUE))
    stop("engine = '", engine, "' needs mgcv", call. = FALSE)
  fit <- switch(engine,
    gam = mgcv::gam(f, family = family, data = d),
    bam = mgcv::bam(f, family = family, data = d, discrete = TRUE),
    glm = stats::glm(f, family = family, data = d))
  co <- if (engine == "glm") summary(fit)$coefficients else summary(fit)$p.table
  ex_rows <- grep(paste0("^", exposure), rownames(co))
  rel <- vapply(adj, function(n) cm$measurements[[n]]$reliability, numeric(1))
  structure(list(exposure = exposure, outcome = outcome, effect = effect,
                 estimate = co[ex_rows, 1][1], se = co[ex_rows, 2][1],
                 p = co[ex_rows, 4][1],
                 adjustment_set = adj, adjuster_reliability = rel,
                 alternatives = sets, fit = fit,
                 attenuation_note = if (length(rel) && min(rel) < 1)
                   sprintf(paste0("Adjuster reliability as low as %.2f: ",
                                  "residual confounding through imperfectly ",
                                  "measured adjusters may bias this estimate."),
                           min(rel)) else NULL),
            class = "causal_effect")
}

#' @export
print.causal_effect <- function(x, ...) {
  cat(sprintf("%s -> %s (%s effect): %+.4f (se %.4f, p %.3g)\n",
              x$exposure, x$outcome, x$effect, x$estimate, x$se, x$p))
  cat(sprintf("adjusted for: %s\n",
              if (length(x$adjustment_set)) paste(
                sprintf("%s [rel %.2f]", x$adjustment_set,
                        x$adjuster_reliability), collapse = ", ")
              else "(nothing — no open backdoor)"))
  if (!is.null(x$attenuation_note)) cat(x$attenuation_note, "\n")
  invisible(x)
}

#' Interactive HTML report for a causal model
#'
#' Renders the DAG as a clickable SVG using the dagitty layout (draw the
#' graph at dagitty.net and the coordinates carry through), with node
#' measurement provenance, ranked adjustment sets, DAG-data tests, and
#' estimated effects. Self-contained file, no external assets.
#'
#' @param cm A causal_model.
#' @param effects Optional named list of causal_effect objects.
#' @param file Output path ("causal_report.html").
#' @param title Report title.
#' @export
causal_report <- function(cm, effects = list(), file = "causal_report.html",
                          title = "Causal model report") {
  g <- cm$dag
  co <- dagitty::coordinates(g)
  if (anyNA(co$x) || length(co$x) == 0) {
    g <- dagitty::graphLayout(g); co <- dagitty::coordinates(g)
  }
  nodes <- names(co$x)
  ex <- dagitty::exposures(g); oc <- dagitty::outcomes(g)
  # scale layout into a 900x520 viewBox (dagitty y grows downward, like SVG)
  rx <- range(co$x); ry <- range(co$y)
  sx <- function(x) 80 + (x - rx[1]) / max(diff(rx), 1e-9) * 740
  sy <- function(y) 70 + (y - ry[1]) / max(diff(ry), 1e-9) * 380
  role <- function(n) {
    if (n %in% ex) "exposure" else if (n %in% oc) "outcome"
    else if (n %in% cm$latent) "latent" else "measured"
  }
  esc <- function(s) gsub("<", "&lt;", gsub("&", "&amp;", s))

  ed <- dagitty::edges(g)
  edge_svg <- paste(vapply(seq_len(nrow(ed)), function(i) {
    x1 <- sx(co$x[ed$v[i]]); y1 <- sy(co$y[ed$v[i]])
    x2 <- sx(co$x[ed$w[i]]); y2 <- sy(co$y[ed$w[i]])
    dx <- x2 - x1; dy <- y2 - y1; L <- sqrt(dx^2 + dy^2)
    # retract ends so arrowheads sit at the node border, not its center
    x1r <- x1 + dx / L * 46; y1r <- y1 + dy / L * 20
    x2r <- x2 - dx / L * 46; y2r <- y2 - dy / L * 20
    if (ed$e[i] == "<->") {
      mx <- (x1r + x2r) / 2 - dy / L * 40; my <- (y1r + y2r) / 2 + dx / L * 40
      sprintf('<path d="M%.0f %.0f Q%.0f %.0f %.0f %.0f" class="edge bidir" marker-start="url(#arr)" marker-end="url(#arr)"/>',
              x1r, y1r, mx, my, x2r, y2r)
    } else
      sprintf('<line x1="%.0f" y1="%.0f" x2="%.0f" y2="%.0f" class="edge" marker-end="url(#arr)"/>',
              x1r, y1r, x2r, y2r)
  }, character(1)), collapse = "\n")

  node_svg <- paste(vapply(nodes, function(n) {
    sprintf(paste0('<g class="node %s" data-node="%s" tabindex="0">',
                   '<ellipse cx="%.0f" cy="%.0f" rx="52" ry="20"/>',
                   '<text x="%.0f" y="%.0f">%s</text></g>'),
            role(n), esc(n), sx(co$x[n]), sy(co$y[n]),
            sx(co$x[n]), sy(co$y[n]) + 4, esc(n))
  }, character(1)), collapse = "\n")

  node_info <- lapply(stats::setNames(nodes, nodes), function(n) {
    m <- cm$measurements[[n]]
    list(role = role(n),
         source = if (is.null(m)) "unmeasured (latent)" else m$source,
         reliability = if (is.null(m)) NA else m$reliability)
  })
  eff_info <- lapply(effects, function(e) {
    lo <- e$estimate - 1.96 * e$se; hi <- e$estimate + 1.96 * e$se
    lnk <- tryCatch(e$fit$family$link, error = function(x) NA)
    list(exposure = e$exposure, outcome = e$outcome, effect = e$effect,
         estimate = e$estimate, lo = lo, hi = hi, p = e$p,
         pct = if (identical(lnk, "log")) 100 * (exp(e$estimate) - 1) else NA,
         adjustment = e$adjustment_set,
         reliability = as.list(e$adjuster_reliability),
         note = e$attenuation_note)
  })
  adj <- causal_adjustment_sets(cm, ex[1], oc[1])
  tests <- causal_dag_test(cm)
  payload <- list(nodes = node_info, effects = eff_info,
                  adjustment_sets = adj,
                  tests = if ("estimate" %in% names(tests))
                    cbind(implication = rownames(tests), tests) else tests)
  json <- jsonlite::toJSON(payload, auto_unbox = TRUE, na = "null", digits = 4)

  html <- paste0('<title>', esc(title), '</title>
<style>
:root{--bg:#FAF9F6;--ink:#1E1B16;--mut:#6B675F;--card:#FFFFFF;--line:#E5E1D8;
--acc:#2F6B4F;--out:#7A4A2B;--warn:#A33B2E}
:root:not([data-theme="light"]){}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#171512;
--ink:#EDEAE3;--mut:#9B968C;--card:#211E1A;--line:#37332C;--acc:#5FA383;--out:#C08A5F;--warn:#D06A5B}}
:root[data-theme="dark"]{--bg:#171512;--ink:#EDEAE3;--mut:#9B968C;--card:#211E1A;
--line:#37332C;--acc:#5FA383;--out:#C08A5F;--warn:#D06A5B}
body{background:var(--bg);color:var(--ink);font:16px/1.55 "Iowan Old Style",Georgia,serif;
margin:0;padding:2rem 1rem}
main{max-width:60rem;margin:0 auto;display:flex;flex-direction:column;gap:1.2rem}
h1{font-size:1.7rem;margin:0}h2{font-size:1.1rem;margin:0 0 .5rem}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:1rem 1.2rem}
svg{width:100%;height:auto;display:block}
.edge{stroke:var(--mut);stroke-width:1.6;fill:none}
.bidir{stroke-dasharray:5 4}
.node ellipse{fill:var(--card);stroke:var(--mut);stroke-width:1.5}
.node.exposure ellipse{stroke:var(--acc);stroke-width:2.5}
.node.outcome ellipse{stroke:var(--out);stroke-width:2.5}
.node.latent ellipse{stroke-dasharray:4 3}
.node text{text-anchor:middle;font:13px "Iowan Old Style",Georgia,serif;fill:var(--ink)}
.node{cursor:pointer}.node:hover ellipse,.node:focus ellipse{fill:var(--line)}
#detail{min-height:3.2rem;color:var(--mut)}
#detail b{color:var(--ink)}
table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}
th,td{text-align:left;padding:.3rem .6rem;border-bottom:1px solid var(--line);font-size:.92rem}
.flag{color:var(--warn);font-weight:bold}
.mono{color:var(--mut);font-size:.85rem}
</style>
<main>
<h1>', esc(title), '</h1>
<div class="card"><h2>Causal diagram</h2>
<svg viewBox="0 0 900 520" role="img" aria-label="causal DAG">
<defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7"
markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="var(--mut)"/></marker></defs>
', edge_svg, '
', node_svg, '
</svg>
<div id="detail">Click a node for its measurement provenance.</div></div>
<div class="card"><h2>Adjustment sets (backdoor, ranked by reliability)</h2><div id="adj"></div></div>
<div class="card"><h2>DAG&ndash;data consistency</h2><div id="tests"></div></div>
<div class="card"><h2>Estimated effects</h2><div id="effects"></div></div>
</main>
<script>
const D = ', json, ';
const fmt = (x,d=3) => x==null ? "&ndash;" : Number(x).toFixed(d);
document.querySelectorAll(".node").forEach(el => el.addEventListener("click", () => {
  const n = el.dataset.node, i = D.nodes[n];
  document.getElementById("detail").innerHTML =
    `<b>${n}</b> &mdash; ${i.role}. Measured by: ${i.source}` +
    (i.reliability!=null ? ` (reliability ${fmt(i.reliability,2)})` : "");
}));
document.getElementById("adj").innerHTML = "<table><tr><th>set</th><th>usable</th>"+
 "<th>min reliability</th><th>mean</th></tr>" + D.adjustment_sets.map(s =>
 `<tr><td>${s.set}</td><td>${s.usable?"yes":"no"}</td><td>${fmt(s.min_reliability,2)}</td>`+
 `<td>${fmt(s.mean_reliability,2)}</td></tr>`).join("") + "</table>";
document.getElementById("tests").innerHTML = (D.tests.length && D.tests[0].implication) ?
 "<table><tr><th>implied independence</th><th>estimate</th><th>95% CI</th><th></th></tr>" +
 D.tests.map(t => `<tr><td>${t.implication}</td><td>${fmt(t.estimate)}</td>`+
 `<td>[${fmt(t["2.5%"])}, ${fmt(t["97.5%"])}]</td>`+
 `<td>${t.contradicts_dag?"<span class=flag>contradicts DAG</span>":"consistent"}</td></tr>`).join("")+
 "</table>" : `<p class="mono">${D.tests.note ?? (D.tests[0]&&D.tests[0].note) ?? "no testable implications among measured nodes"}</p>`;
document.getElementById("effects").innerHTML = Object.keys(D.effects).length ?
 "<table><tr><th>panel</th><th>effect</th><th>estimate</th><th>95% CI</th><th>%</th><th>p</th><th>adjusted for</th></tr>"+
 Object.entries(D.effects).map(([k,e]) =>
 `<tr><td>${k}</td><td>${e.exposure} &rarr; ${e.outcome} (${e.effect})</td>`+
 `<td>${fmt(e.estimate)}</td><td>[${fmt(e.lo)}, ${fmt(e.hi)}]</td>`+
 `<td>${e.pct!=null?fmt(e.pct,1)+"%":"&ndash;"}</td><td>${e.p<0.001?e.p.toExponential(1):fmt(e.p)}</td>`+
 `<td class="mono">${(e.adjustment||[]).map(a=>`${a} [${fmt(e.reliability[a],2)}]`).join(", ")}</td></tr>`).join("")+
 "</table>" + Object.values(D.effects).filter(e=>e.note).map(e=>`<p class="mono">${e.note}</p>`).join("")
 : `<p class="mono">no effects estimated yet</p>`;
</script>
')
  writeLines(html, file)
  invisible(file)
}
