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

const HEAVY_PACKAGES = [
  "torch", "mace", "mace_torch", "ase", "pymatgen",
  "torch_geometric", "botorch", "gpytorch", "mp_api"
];

let pyodideReadyPromise = null;

function loadPyodideOnce() {
  if (pyodideReadyPromise) return pyodideReadyPromise;
  pyodideReadyPromise = (async () => {
    const py = await loadPyodide({
      indexURL: "https://cdn.jsdelivr.net/pyodide/v0.26.0/full/"
    });
    await py.loadPackage(["numpy", "scipy", "matplotlib"]);
    py.runPython(`
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

      py.runPython("import sys, io; sys.stdout = io.StringIO(); sys.stderr = sys.stdout");
      try {
        await py.runPythonAsync(code);
      } catch (e) {
        const stderr = py.runPython("sys.stdout.getvalue()");
        output.innerHTML =
          (stderr ? "<pre class='pyodide-stdout'>" + escapeHtml(stderr) + "</pre>" : "") +
          "<pre class='pyodide-error'>" + escapeHtml(e.toString()) + "</pre>";
        btn.disabled = false;
        btn.innerHTML = "<span class='pyodide-icon'>▶</span> Run";
        return;
      }

      const stdout = py.runPython("sys.stdout.getvalue()");
      let plotB64 = "";
      try {
        plotB64 = py.runPython("__capture_plot()");
      } catch (_) {}

      let html = "";
      if (stdout) html += "<pre class='pyodide-stdout'>" + escapeHtml(stdout) + "</pre>";
      if (plotB64) html += "<img class='pyodide-plot' src='data:image/png;base64," + plotB64 + "' alt='matplotlib output'>";
      if (!html) html = "<div class='pyodide-note'>Ran successfully — no output.</div>";

      output.innerHTML = html;
      btn.disabled = false;
      btn.innerHTML = "<span class='pyodide-icon'>▶</span> Run again";
    } catch (err) {
      output.innerHTML = "<pre class='pyodide-error'>" + escapeHtml(err.toString()) + "</pre>";
      btn.disabled = false;
      btn.innerHTML = "<span class='pyodide-icon'>▶</span> Run";
    }
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
      attachRunButton(el);
    } else {
      // <code> element — promote to its <pre>'s parent (or the pre itself if no wrapper)
      const pre = el.closest("pre");
      if (pre) {
        const parent = pre.parentElement;
        if (parent && parent.classList && parent.classList.contains("language-python")) {
          attachRunButton(parent);
        } else {
          attachRunButton(pre);
        }
      }
    }
  });
}

if (document.readyState !== "loading") {
  initRunButtons();
} else {
  document.addEventListener("DOMContentLoaded", initRunButtons);
}

// Material for MkDocs uses navigation.instant — re-attach on page change.
if (typeof document$ !== "undefined") {
  document$.subscribe(() => initRunButtons());
}
