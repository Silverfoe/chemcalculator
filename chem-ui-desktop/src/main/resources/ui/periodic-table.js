const CATEGORY_COLORS = {
  "Alkali Metal": "#f8c7d4",
  "Alkaline Earth Metal": "#f7ddb8",
  "Transition Metal": "#f7ee9d",
  "Post-Transition Metal": "#b8dceb",
  "Metalloid": "#d8d8d8",
  "Nonmetal": "#fbf5c4",
  "Halogen": "#bfe9bf",
  "Noble Gas": "#ccf1f3",
  "Lanthanide": "#d5b0e4",
  "Actinide": "#d9c3e4",
  "Metal": "#d9e7f7"
};

const legendElement = document.getElementById("legend");
const statsElement = document.getElementById("stats");
const mainGrid = document.getElementById("main-grid");
const lanthanides = document.getElementById("lanthanides");
const actinides = document.getElementById("actinides");
const emptyState = document.getElementById("empty-state");

function createLegend() {
  legendElement.innerHTML = "";
  Object.entries(CATEGORY_COLORS).forEach(([label, color]) => {
    const chip = document.createElement("div");
    chip.className = "legend-chip";
    chip.style.background = color;
    chip.textContent = label;
    legendElement.appendChild(chip);
  });
}

function tileColor(category) {
  return CATEGORY_COLORS[category] || "#edf2f7";
}

function createTile(element, selectedSymbol) {
  const button = document.createElement("button");
  button.className = "tile";
  if (selectedSymbol && selectedSymbol.toLowerCase() === element.symbol.toLowerCase()) {
    button.classList.add("selected");
  }
  button.style.background = tileColor(element.category);
  button.title = `${element.name} (${element.symbol})\n${element.category}\n${element.standardState}\n${element.atomicWeight}`;
  button.innerHTML = `
    <div class="tile-number">${element.atomicNumber}</div>
    <div class="tile-symbol">${element.symbol}</div>
    <div class="tile-name">${element.name}</div>
    <div class="tile-meta">${element.atomicWeight}</div>
  `;
  button.addEventListener("click", () => {
    if (window.periodicBridge && typeof window.periodicBridge.selectElement === "function") {
      window.periodicBridge.selectElement(element.symbol);
    }
  });
  return button;
}

function renderSeries(container, elements, selectedSymbol) {
  container.innerHTML = "";
  elements.forEach(element => {
    container.appendChild(createTile(element, selectedSymbol));
  });
}

window.renderTable = function(elements, selectedSymbol) {
  mainGrid.innerHTML = "";

  const mainElements = [];
  const lanthanideElements = [];
  const actinideElements = [];

  elements.forEach(element => {
    if (element.category === "Lanthanide") {
      lanthanideElements.push(element);
    } else if (element.category === "Actinide") {
      actinideElements.push(element);
    } else {
      mainElements.push(element);
    }
  });

  mainElements.forEach(element => {
    const tile = createTile(element, selectedSymbol);
    if (element.group !== null) {
      tile.style.gridColumn = String(element.group);
      tile.style.gridRow = String(element.period);
    }
    mainGrid.appendChild(tile);
  });

  renderSeries(lanthanides, lanthanideElements, selectedSymbol);
  renderSeries(actinides, actinideElements, selectedSymbol);

  if (!elements.length) {
    emptyState.classList.remove("hidden");
  } else {
    emptyState.classList.add("hidden");
  }

  const selectedLabel = selectedSymbol ? ` Selected: ${selectedSymbol}.` : "";
  statsElement.textContent = `${elements.length} element${elements.length === 1 ? "" : "s"} visible.${selectedLabel}`;
};

createLegend();
