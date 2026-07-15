// run_python.js
// Adds an inline "Run" button to every Python code block on the page.
// Uses Pyodide (WebAssembly Python) to execute in the browser, no server.
//
// MkDocs Material + pymdownx.highlight emits this structure:
//   <div class="language-python highlight">
//     <pre><span></span><code>...code with line anchors...</code></pre>
//   </div>
// So we hook onto the outer .language-python container.
//
// Heavy packages (torch, mace, ase, pymatgen, torch_geometric) are not
// available in Pyodide. For code that imports those, the button shows a
// link pointing the reader to Google Colab instead.
//
// --- Widget mode -----------------------------------------------------------
// Authors can opt a code block into "widget mode" by preceding it with a
// YAML config block whose first line is the marker `# widget-config`.
// Example markdown:
//
//   ```yaml
//   # widget-config
//   sliders:
//     L:     {min: 0.5, max: 5.0, step: 0.1, default: 1.0, label: "Box L (Å)"}
//     n_max: {min: 1,   max: 8,   step: 1,   default: 4,   label: "States"}
//   ```
//
//   ```python
//   # widget — uses L, n_max bound from the sliders above
//   import numpy as np
//   ...
//   ```
//
// The python block must also begin with a `# widget` line so that the
// pairing is unambiguous and existing tutorial code is unaffected.

const HEAVY_PACKAGES = [
  "torch", "mace", "mace_torch", "ase", "pymatgen",
  "torch_geometric", "botorch", "gpytorch", "mp_api",
  // Not packaged for Pyodide either — point readers to Colab rather than
  // letting the run fail with a confusing ModuleNotFoundError.
  "MDAnalysis", "pymbar", "umap", "umap_learn"
];

let pyodideReadyPromise = null;

// Pyodide is fetched at runtime. We prefer a SAME-ORIGIN self-hosted copy
// (docs/assets/pyodide/ — core + the numpy/scipy/matplotlib dependency
// closure), which is guaranteed reachable whenever the page itself loaded, so
// the Run feature no longer depends on any third-party CDN. If that copy is
// somehow unavailable we fall back to several *complete* jsDelivr backends.
// cdn.jsdelivr.net alone is unreliable on some networks — notably mainland
// China, where its DNS is intermittently polluted, surfacing as:
//   "Failed to fetch dynamically imported module: .../pyodide.asm.js"
// Each fallback is a different hostname mirroring the FULL build, so the wheels
// loadPackage() needs live alongside the core (npmmirror/cdnjs ship core only).

// Derive the self-hosted base from this script's own URL, so it works no matter
// what sub-path the site is deployed under.
const LOCAL_PYODIDE_BASE = (function () {
  try {
    const scripts = document.getElementsByTagName("script");
    for (const s of scripts) {
      if (s.src && s.src.indexOf("assets/javascripts/run_python.js") !== -1) {
        return s.src.replace(/assets\/javascripts\/run_python\.js.*$/, "assets/pyodide/");
      }
    }
  } catch (e) { /* fall through to CDNs */ }
  return null;
})();

const PYODIDE_BASES = [
  ...(LOCAL_PYODIDE_BASE ? [LOCAL_PYODIDE_BASE] : []),
  "https://gcore.jsdelivr.net/pyodide/v0.26.0/full/",
  "https://fastly.jsdelivr.net/pyodide/v0.26.0/full/",
  "https://cdn.jsdelivr.net/pyodide/v0.26.0/full/"
];

// Inject a <script src> and resolve once it has executed (reject on error).
function loadScript(src) {
  return new Promise((resolve, reject) => {
    const s = document.createElement("script");
    s.src = src;
    s.async = true;
    s.onload = () => resolve();
    s.onerror = () => { s.remove(); reject(new Error("failed to load " + src)); };
    document.head.appendChild(s);
  });
}

// Probe a mirror with a bounded timeout so a blocked/polluted host is skipped
// quickly instead of hanging on DNS. jsDelivr sends permissive CORS headers,
// so this cross-origin fetch is allowed.
async function mirrorReachable(base, timeoutMs) {
  const ctrl = new AbortController();
  const timer = setTimeout(() => ctrl.abort(), timeoutMs);
  try {
    const r = await fetch(base + "pyodide.js", { signal: ctrl.signal });
    return r.ok;
  } catch (e) {
    return false;
  } finally {
    clearTimeout(timer);
  }
}

function loadPyodideOnce() {
  if (pyodideReadyPromise) return pyodideReadyPromise;
  pyodideReadyPromise = (async () => {
    // Pick the first mirror that actually responds, then load everything from
    // it. Probing first means loadPyodide() is called exactly once (calling it
    // repeatedly across mirrors is not supported).
    let base = null;
    for (const candidate of PYODIDE_BASES) {
      if (await mirrorReachable(candidate, 8000)) { base = candidate; break; }
      console.warn("[run_python] Pyodide mirror unreachable, trying next: " + candidate);
    }
    if (!base) {
      pyodideReadyPromise = null; // allow a later click to retry from scratch
      throw new Error(
        "Could not reach any Pyodide mirror (gcore / fastly / cdn jsDelivr). " +
        "This is usually a temporary network block — check your connection or " +
        "VPN/proxy and click Run again."
      );
    }

    if (typeof loadPyodide === "undefined") {
      await loadScript(base + "pyodide.js");
    }
    const py = await loadPyodide({ indexURL: base });
    // Preload the full set of Pyodide packages used anywhere in the handbook
    // (numpy 92, matplotlib 45, scipy 14, pandas 4 blocks) once, so every "Run"
    // works regardless of which one the reader clicks. loadPackagesFromImports
    // below then covers any snippet that needs something beyond this base set.
    await py.loadPackage(["numpy", "scipy", "matplotlib", "pandas"]);
    py.runPython(`
# Pre-bind the conventional aliases (np, pd, plt) and the bare module names, so a
# snippet that uses them WITHOUT an explicit import still runs — e.g. a fragment
# that just does "u = np.array(...)". Every Run block shares one namespace, so
# binding these once at startup covers all of them.
import numpy as np
import numpy
import scipy
import pandas as pd
import pandas
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import io, base64, sys

def __capture_plot():
    if not plt.get_fignums():
        return ""
    buf = io.BytesIO()
    plt.gcf().savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close("all")
    return base64.b64encode(buf.getvalue()).decode()
`);
    return py;
  })();
  return pyodideReadyPromise;
}

function escapeHtml(s) {
  return String(s)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;");
}

function usesHeavyPackage(code) {
  for (const pkg of HEAVY_PACKAGES) {
    const re = new RegExp("(^|\\W)(import\\s+" + pkg + "|from\\s+" + pkg + ")(\\W|$)", "m");
    if (re.test(code)) return pkg;
  }
  return null;
}

// ---------------------------------------------------------------------------
// Tiny YAML subset parser
// ---------------------------------------------------------------------------
// We only support the exact shape we document for widget-config:
//
//   sliders:
//     NAME: {min: A, max: B, step: C, default: D, label: "..."}
//     ...
//
// This is enough to avoid pulling in js-yaml and keeps the page offline-safe.
function parseWidgetConfig(text) {
  const sliders = [];
  const lines = text.split(/\r?\n/);
  let inSliders = false;
  for (let raw of lines) {
    const line = raw.replace(/\t/g, "  ");
    if (!line.trim() || line.trim().startsWith("#")) continue;
    if (/^sliders\s*:\s*$/.test(line)) {
      inSliders = true;
      continue;
    }
    if (!inSliders) continue;
    // A slider entry: "  NAME: { ... }"
    const m = line.match(/^\s+([A-Za-z_][A-Za-z0-9_]*)\s*:\s*\{(.+)\}\s*$/);
    if (!m) continue;
    const name = m[1];
    const body = m[2];
    const spec = { name: name, min: 0, max: 1, step: 0.01, default: 0, label: name };
    // Pull key: value pairs. Values can be numbers or quoted strings.
    const re = /([A-Za-z_]+)\s*:\s*("([^"]*)"|'([^']*)'|[-+0-9.eE]+)/g;
    let mm;
    while ((mm = re.exec(body)) !== null) {
      const key = mm[1];
      const valStr = mm[2];
      let value;
      if (mm[3] !== undefined) value = mm[3];
      else if (mm[4] !== undefined) value = mm[4];
      else value = parseFloat(valStr);
      spec[key] = value;
    }
    sliders.push(spec);
  }
  return { sliders: sliders };
}

// ---------------------------------------------------------------------------
// Pyodide execution helper (shared by run-button and widget paths)
// ---------------------------------------------------------------------------
async function runPythonInto(py, code, output, globalsObj) {
  py.runPython("import sys, io; sys.stdout = io.StringIO(); sys.stderr = sys.stdout");

  // Inject slider values as module-level globals via Pyodide's globals proxy.
  if (globalsObj) {
    for (const [k, v] of Object.entries(globalsObj)) {
      py.globals.set(k, v);
    }
  }

  try {
    // Auto-load any Pyodide packages this snippet imports, so a run never fails
    // just because the code imports a module without explicitly loading it.
    // Packages already preloaded are a no-op; anything else in the self-hosted
    // set is fetched on demand and cached for later runs.
    await py.loadPackagesFromImports(code);
    await py.runPythonAsync(code);
  } catch (e) {
    const stderr = py.runPython("sys.stdout.getvalue()");
    output.innerHTML =
      (stderr ? "<pre class='pyodide-stdout'>" + escapeHtml(stderr) + "</pre>" : "") +
      "<pre class='pyodide-error'>" + escapeHtml(e.toString()) + "</pre>";
    return false;
  }

  const stdout = py.runPython("sys.stdout.getvalue()");
  let plotB64 = "";
  try {
    plotB64 = py.runPython("__capture_plot()");
  } catch (_) { /* ignore */ }

  let html = "";
  if (stdout) html += "<pre class='pyodide-stdout'>" + escapeHtml(stdout) + "</pre>";
  if (plotB64) html += "<img class='pyodide-plot' src='data:image/png;base64," + plotB64 + "' alt='matplotlib output'>";
  if (!html) html = "<div class='pyodide-note'>Ran successfully — no output.</div>";
  output.innerHTML = html;
  return true;
}

// ---------------------------------------------------------------------------
// Plain Run button
// ---------------------------------------------------------------------------
function attachRunButton(blockEl) {
  // blockEl is the <div class="language-python ..."> wrapper.
  if (blockEl.dataset.runAttached === "1") return;
  blockEl.dataset.runAttached = "1";

  const codeEl = blockEl.querySelector("pre > code");
  if (!codeEl) return;

  // Ensure positioning context for absolute button
  blockEl.style.position = "relative";
  blockEl.classList.add("pyodide-host");

  // Run button — anchored to top-right of the code block
  const btn = document.createElement("button");
  btn.className = "pyodide-run-btn";
  btn.type = "button";
  btn.innerHTML = "<span class='pyodide-icon'>▶</span> Run";
  blockEl.appendChild(btn);

  // Output panel inserted AFTER the code block, sibling-wise
  const output = document.createElement("div");
  output.className = "pyodide-output";
  output.style.display = "none";
  blockEl.parentNode.insertBefore(output, blockEl.nextSibling);

  btn.addEventListener("click", async () => {
    const code = codeEl.textContent;
    const heavy = usesHeavyPackage(code);
    if (heavy) {
      output.style.display = "block";
      output.innerHTML =
        "<div class='pyodide-note'>This snippet imports <code>" +
        heavy +
        "</code>, which is too heavy for in-browser Python.<br>" +
        "Open the corresponding chapter in <a href='https://colab.research.google.com/' target='_blank' rel='noopener'>Google Colab</a> to run it with the full scientific stack (torch / ase / pymatgen / mace).</div>";
      return;
    }

    btn.disabled = true;
    btn.innerHTML = "<span class='pyodide-icon'>⏳</span> Loading…";
    output.style.display = "block";
    output.innerHTML =
      "<div class='pyodide-note'>Loading Python runtime in your browser. First run takes ~10 s while NumPy / SciPy / matplotlib are fetched and cached. Subsequent runs are instant.</div>";

    try {
      const py = await loadPyodideOnce();
      btn.innerHTML = "<span class='pyodide-icon'>⚙</span> Running…";
      const ok = await runPythonInto(py, code, output, null);
      btn.disabled = false;
      btn.innerHTML = ok
        ? "<span class='pyodide-icon'>▶</span> Run again"
        : "<span class='pyodide-icon'>▶</span> Run";
    } catch (err) {
      output.innerHTML = "<pre class='pyodide-error'>" + escapeHtml(err.toString()) + "</pre>";
      btn.disabled = false;
      btn.innerHTML = "<span class='pyodide-icon'>▶</span> Run";
    }
  });
}

// ---------------------------------------------------------------------------
// Widget mode: sliders that drive a python re-run
// ---------------------------------------------------------------------------
function attachWidget(yamlBlockEl, pythonBlockEl, config) {
  if (pythonBlockEl.dataset.widgetAttached === "1") return;
  pythonBlockEl.dataset.widgetAttached = "1";
  yamlBlockEl.dataset.widgetAttached = "1";

  const codeEl = pythonBlockEl.querySelector("pre > code");
  if (!codeEl) return;

  // Hide the yaml config block — readers do not need to see the slider spec.
  yamlBlockEl.style.display = "none";

  // Build the slider panel and place it BEFORE the python block.
  const panel = document.createElement("div");
  panel.className = "pyodide-widget-panel";

  const header = document.createElement("div");
  header.className = "pyodide-widget-header";
  header.innerHTML =
    "<span class='pyodide-widget-title'>Interactive — drag to update</span>" +
    "<span class='pyodide-widget-status'>idle</span>";
  panel.appendChild(header);

  const status = header.querySelector(".pyodide-widget-status");
  const inputs = {};

  for (const s of config.sliders) {
    const row = document.createElement("label");
    row.className = "pyodide-slider-row";

    const labelText = document.createElement("span");
    labelText.className = "pyodide-slider-label";
    labelText.textContent = s.label || s.name;
    row.appendChild(labelText);

    const input = document.createElement("input");
    input.type = "range";
    input.min = String(s.min);
    input.max = String(s.max);
    input.step = String(s.step);
    input.value = String(s.default);
    input.className = "pyodide-slider";
    row.appendChild(input);

    const valOut = document.createElement("span");
    valOut.className = "pyodide-slider-value";
    valOut.textContent = formatVal(parseFloat(input.value), s.step);
    row.appendChild(valOut);

    input.addEventListener("input", () => {
      valOut.textContent = formatVal(parseFloat(input.value), s.step);
      scheduleRun();
    });

    inputs[s.name] = { input: input, spec: s };
    panel.appendChild(row);
  }

  pythonBlockEl.parentNode.insertBefore(panel, pythonBlockEl);

  // Output goes after the python block.
  const output = document.createElement("div");
  output.className = "pyodide-output pyodide-widget-output";
  pythonBlockEl.parentNode.insertBefore(output, pythonBlockEl.nextSibling);

  let pending = null;
  let running = false;
  let queued = false;

  function gatherGlobals() {
    const g = {};
    for (const [name, entry] of Object.entries(inputs)) {
      const v = parseFloat(entry.input.value);
      // Snap to int if the step is integer.
      if (Number.isInteger(entry.spec.step) && Number.isInteger(entry.spec.min)) {
        g[name] = Math.round(v);
      } else {
        g[name] = v;
      }
    }
    return g;
  }

  async function doRun() {
    if (running) { queued = true; return; }
    running = true;
    status.textContent = "running…";
    try {
      const py = await loadPyodideOnce();
      const code = codeEl.textContent;
      const heavy = usesHeavyPackage(code);
      if (heavy) {
        output.innerHTML =
          "<div class='pyodide-note'>This widget uses <code>" + heavy +
          "</code>, which is not available in browser Python.</div>";
        status.textContent = "skipped";
        running = false;
        return;
      }
      const ok = await runPythonInto(py, code, output, gatherGlobals());
      status.textContent = ok ? "ready" : "error";
    } catch (err) {
      output.innerHTML = "<pre class='pyodide-error'>" + escapeHtml(err.toString()) + "</pre>";
      status.textContent = "error";
    } finally {
      running = false;
      if (queued) {
        queued = false;
        scheduleRun();
      }
    }
  }

  function scheduleRun() {
    if (pending) clearTimeout(pending);
    pending = setTimeout(doRun, 150);
  }

  // Kick off an initial run (lazy — only when reader scrolls near it would
  // be nicer, but Pyodide loadOnce already debounces the heavy work and
  // most readers will engage with the widget immediately).
  status.textContent = "loading runtime…";
  doRun();
}

function formatVal(v, step) {
  // Very large or very small magnitudes — use scientific notation.
  const absv = Math.abs(v);
  if (absv !== 0 && (absv >= 1e4 || absv < 1e-3)) {
    return v.toExponential(2);
  }
  if (Number.isInteger(step)) return String(Math.round(v));
  const decimals = Math.max(0, Math.min(6, -Math.floor(Math.log10(step))));
  return v.toFixed(decimals);
}

// ---------------------------------------------------------------------------
// First-line scan: returns the first non-empty text line of a code block.
// ---------------------------------------------------------------------------
function firstLine(blockEl) {
  const codeEl = blockEl.querySelector("pre > code");
  if (!codeEl) return "";
  const text = codeEl.textContent || "";
  for (const line of text.split(/\r?\n/)) {
    if (line.trim()) return line.trim();
  }
  return "";
}

function isWidgetYaml(blockEl) {
  return /language-yaml/.test(blockEl.className) &&
         firstLine(blockEl) === "# widget-config";
}

function isWidgetPython(blockEl) {
  return /language-python/.test(blockEl.className) &&
         /^#\s*widget\b/.test(firstLine(blockEl));
}

function nextElementBlock(el) {
  let n = el.nextElementSibling;
  // Skip empty text-ish wrappers if any.
  while (n && n.nodeType === 1 && n.children.length === 0 && !n.textContent.trim()) {
    n = n.nextElementSibling;
  }
  return n;
}

// ---------------------------------------------------------------------------
// Init
// ---------------------------------------------------------------------------
function initWidgets() {
  // Find every yaml block tagged as widget-config and look at its next sibling.
  const yamls = document.querySelectorAll("div[class*='language-yaml']");
  yamls.forEach((y) => {
    if (!isWidgetYaml(y)) return;
    const nxt = nextElementBlock(y);
    if (!nxt || !isWidgetPython(nxt)) return;
    const codeText = (y.querySelector("pre > code") || {}).textContent || "";
    let config;
    try {
      config = parseWidgetConfig(codeText);
    } catch (e) {
      console.error("widget-config parse error:", e);
      return;
    }
    if (!config.sliders || config.sliders.length === 0) return;
    attachWidget(y, nxt, config);
  });
}

function initRunButtons() {
  // pymdownx.highlight emits <div class="language-python highlight"><pre><code>...
  // Fallback for plain pandoc-style <pre><code class="language-python">.
  const containers = document.querySelectorAll(
    "div[class*='language-python'], pre > code.language-python, pre > code.python"
  );
  containers.forEach((el) => {
    if (el.tagName === "DIV") {
      // Skip blocks that are already wired up as widgets.
      if (el.dataset.widgetAttached === "1") return;
      attachRunButton(el);
    } else {
      // <code> element — promote to its <pre>'s parent (or the pre itself if no wrapper)
      const pre = el.closest("pre");
      if (pre) {
        const parent = pre.parentElement;
        if (parent && parent.classList && parent.classList.contains("language-python")) {
          if (parent.dataset.widgetAttached === "1") return;
          attachRunButton(parent);
        } else {
          attachRunButton(pre);
        }
      }
    }
  });
}

function initAll() {
  // Widgets first, so the run-button pass can skip widget-bound python blocks.
  initWidgets();
  initRunButtons();
}

if (document.readyState !== "loading") {
  initAll();
} else {
  document.addEventListener("DOMContentLoaded", initAll);
}

// Material for MkDocs uses navigation.instant — re-attach on page change.
if (typeof document$ !== "undefined") {
  document$.subscribe(() => initAll());
}
