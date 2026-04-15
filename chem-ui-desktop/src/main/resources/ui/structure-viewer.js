let viewer = null;
let currentPayload = null;
let currentStylePreset = "ball-stick";
let viewerAvailable = false;
const bridgePrefix = "__CHEMCALC__:";

function sendBridge(type, payload) {
  try {
    const message = JSON.stringify({ type, ...(payload || {}) });
    console.info(`${bridgePrefix}${message}`);
  } catch (error) {
    console.info(`${bridgePrefix}{"type":"bridge-error"}`);
  }
}

function ensureViewer() {
  if (viewer) {
    return viewer;
  }
  try {
    viewer = $3Dmol.createViewer("viewer", {
      backgroundColor: "#ffffff"
    });
    viewerAvailable = true;
  } catch (error) {
    viewer = null;
    viewerAvailable = false;
  }
  return viewer;
}

function styleSpec(preset) {
  switch (preset) {
    case "stick":
      return { stick: { radius: 0.18, colorscheme: "Jmol" } };
    case "sphere":
      return { sphere: { scale: 0.32, colorscheme: "Jmol" } };
    case "line":
      return { line: { linewidth: 2.0, colorscheme: "Jmol" } };
    default:
      return {
        stick: { radius: 0.16, colorscheme: "Jmol" },
        sphere: { scale: 0.24, colorscheme: "Jmol" }
      };
  }
}

function updateMessages(messages) {
  const root = document.getElementById("messages");
  root.innerHTML = "";
  (messages || []).forEach(message => {
    const item = document.createElement("div");
    item.className = `message ${message.severity || "INFO"}`;
    item.textContent = message.message;
    root.appendChild(item);
  });
}

function showOverlay(title, detail) {
  const overlay = document.getElementById("overlay");
  document.getElementById("overlay-title").textContent = title;
  document.getElementById("overlay-detail").textContent = detail;
  overlay.classList.remove("hidden");
}

function hideOverlay() {
  document.getElementById("overlay").classList.add("hidden");
}

function renderSvgFallback(payload) {
  const host = document.getElementById("viewer");
  host.classList.add("fallback");
  host.innerHTML = payload.svgData || "";
}

function clearViewerHost() {
  const host = document.getElementById("viewer");
  if (host.classList.contains("fallback")) {
    host.classList.remove("fallback");
    host.innerHTML = "";
  }
}

function installAtomClickBridge() {
  if (!viewer || !viewerAvailable) {
    return;
  }
  viewer.setClickable({}, true, atom => {
    const label = `${atom.elem || "Atom"} #${atom.serial != null ? atom.serial : atom.index}`;
    sendBridge("atom-selected", {
      atomIndex: String(atom.index),
      label
    });
  });
}

window.renderStructure = function(payload) {
  currentPayload = payload;
  const activeViewer = ensureViewer();
  if (!activeViewer || !viewerAvailable) {
    renderSvgFallback(payload);
    updateMessages(payload.messages || []);
    hideOverlay();
    return;
  }
  clearViewerHost();
  activeViewer.removeAllModels();
  activeViewer.addModel(payload.modelData, payload.modelFormat || "mol");
  activeViewer.setStyle({}, styleSpec(currentStylePreset));
  installAtomClickBridge();
  activeViewer.zoomTo();
  activeViewer.render();
  updateMessages(payload.messages || []);
  hideOverlay();
};

window.showViewerMessage = function(title, detail) {
  currentPayload = null;
  const activeViewer = ensureViewer();
  clearViewerHost();
  if (activeViewer && viewerAvailable) {
    activeViewer.removeAllModels();
    activeViewer.render();
  }
  updateMessages([]);
  showOverlay(title, detail);
};

window.setStylePreset = function(preset) {
  currentStylePreset = preset || "ball-stick";
  if (!viewer || !currentPayload || !viewerAvailable) {
    return;
  }
  viewer.setStyle({}, styleSpec(currentStylePreset));
  viewer.render();
};

window.addEventListener("DOMContentLoaded", () => {
  ensureViewer();
  showOverlay("No structure loaded", "Load a formula, SMILES, Molfile/SDF, or InChI to populate the viewer.");
  sendBridge("viewer-ready", {});
});
