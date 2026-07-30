'use strict';

import { $, on, setStatus, state, notifyDataUpdated, onDataUpdated } from './main.js';

/* ── Helpers ─────────────────────────────────────────────── */
const vadd   = (a, b) => a.map((v, i) => v + b[i]);
const vscale = (a, s) => a.map(v => v * s);
const vsum   = a => a.reduce((x, y) => x + y, 0);

const gauss = () => {
  let u, v;
  do { u = Math.random(); } while (!u);
  do { v = Math.random(); } while (!v);
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
};

/* ── Generic RK4 integrator ─────────────────────────────── */
const simRK4 = (grid, deriv, y0, Cfun) => {
  const out = []; let y = [...y0]; out.push(vsum(y));
  for (let i = 1; i < grid.length; i++) {
    const [t0, h] = [grid[i - 1], grid[i] - grid[i - 1]];
    const [C0, Cm, C1] = [Cfun(t0), Cfun(t0 + h / 2), Cfun(t0 + h)];
    const k1 = deriv(y, C0);
    const k2 = deriv(vadd(y, vscale(k1, h / 2)), Cm);
    const k3 = deriv(vadd(y, vscale(k2, h / 2)), Cm);
    const k4 = deriv(vadd(y, vscale(k3, h)),     C1);
    y = vadd(y, vscale(vadd(vadd(k1, vscale(k2, 2)), vadd(vscale(k3, 2), k4)), h / 6));
    out.push(vsum(y));
  }
  return out;
};

/* ── Binding models as ODE derivatives ───────────────────────
   makeDeriv(model, Rmax, gv) where gv(id) returns a numeric value
   for the given input id. gv defaults to reading the live DOM. */
function makeDeriv(model, Rmax, gv) {
  gv = gv || (id => +$(id).value);

  if (model === "langmuir") {
    const ka = gv("ka"), kd = gv("kd");
    return {
      size: 1,
      deriv: (y, C) => { const R = y[0]; return [ka * C * (Rmax - R) - kd * R]; },
      fluxCoef: (y) => ({ a: ka * (Rmax - y[0]), b: -kd * y[0] })
    };
  }
  if (model === "hetLigand") {
    const ka1 = gv("hetka1"), kd1 = gv("hetkd1"),
          ka2 = gv("hetka2"), kd2 = gv("hetkd2"), Rmax2 = gv("Rmax2");
    return {
      size: 2,
      deriv: (y, C) => { const R1 = y[0], R2 = y[1];
        return [ka1 * C * (Rmax - R1) - kd1 * R1, ka2 * C * (Rmax2 - R2) - kd2 * R2]; },
      fluxCoef: (y) => ({ a: ka1 * (Rmax - y[0]) + ka2 * (Rmax2 - y[1]),
                          b: -(kd1 * y[0] + kd2 * y[1]) })
    };
  }
  if (model === "bivAnalyte") {
    const ka1 = gv("bivka1"), kd1 = gv("bivkd1"),
          ka2 = gv("bivka2"), kd2 = gv("bivkd2");
    return {
      size: 2,
      deriv: (y, C) => { const R1 = y[0], R2 = y[1], free = Rmax - R1 - 2 * R2;
        return [2 * ka1 * C * free - kd1 * R1 - ka2 * R1 * free + 2 * kd2 * R2,
                2 * ka2 * R1 * free - 2 * kd2 * R2]; },
      fluxCoef: (y) => { const free = Rmax - y[0] - 2 * y[1];
        return { a: 2 * ka1 * free, b: -kd1 * y[0] }; }
    };
  }
  // two-state conformational change:  A + B <-> AB <-> AB*
  const ka1 = gv("ka1"), kd1 = gv("kd1"), ka2 = gv("ka2"), kd2 = gv("kd2");
  return {
    size: 2,
    deriv: (y, C) => { const AB = y[0], ABs = y[1], free = Rmax - AB - ABs;
      return [ka1 * C * free - kd1 * AB - ka2 * AB + kd2 * ABs, ka2 * AB - kd2 * ABs]; },
    fluxCoef: (y) => { const free = Rmax - y[0] - y[1];
      return { a: ka1 * free, b: -kd1 * y[0] }; }
  };
}

/* ── Mass-transport modifier. Wraps ANY model's derivative. ── */
function withTransport(base, kt) {
  return {
    size: base.size,
    deriv: (y, C) => {
      const { a, b } = base.fluxCoef(y);
      let Cs = (kt * C - b) / (kt + a);
      if (!isFinite(Cs) || Cs < 0) Cs = 0;
      return base.deriv(y, Cs);
    }
  };
}

/* ── Formatting ──────────────────────────────────────────── */
const parseConcs = str =>
  (str || "").split(/[\s,;]+/).map(Number).filter(v => Number.isFinite(v) && v > 0);

export const fmtConc = Cnm => (Cnm >= 1 ? Cnm : Cnm.toPrecision(3)) + " nM";

/* ── Model metadata ──────────────────────────────────────── */
const MODEL_HINTS = {
  langmuir:   "Simplest case: one analyte binding one immobilised ligand.",
  hetLigand:  "Two available binding sites with two completely independent dynamics.",
  twostate:   "Binding followed by a conformational change that locks the complex — note the slow dissociation.",
  bivAnalyte: "The analyte may, at sufficient density, bind two membrane receptors simultaneously."
};

/* ── Inputs that drive the (single) simulated region ────────
   Any change to one of these re-runs the simulation directly
   from the live DOM values — no per-region storage/switching. */
const MODEL_INPUT_IDS = [
  "model", "ka", "kd", "ka1", "kd1", "ka2", "kd2",
  "hetka1", "hetkd1", "hetka2", "hetkd2",
  "bivka1", "bivkd1", "bivka2", "bivkd2",
  "kt", "Rmax", "Rmax2"
];

function simulate() {
  const tA   = +$("tBase").value;
  const tD   = tA + +$("tAssoc").value;
  const tEnd = tD + +$("tDissoc").value;
  const grid = Array.from({ length: Math.round(tEnd) + 1 }, (_, i) => i);

  const noiseOn = $("noiseOn").checked;
  const noiseSd = +$("noiseSd").value || 0;
  const drift   = +$("drift").value   || 0;

  // Analyte concentrations, ascending -> the stack runs lowest to highest.
  const concs = parseConcs($("concSeries").value).sort((a, b) => a - b);

  const model = $("model").value;
  const Rmax  = +$("Rmax").value || 0;

  const base   = makeDeriv(model, Rmax);
  const engine = ($("mtlOn") && $("mtlOn").checked) ? withTransport(base, +$("kt").value || 0) : base;

  const traces = concs.map(Cnm => {
    const C    = Cnm * 1e-9;
    const Cfun = t => (t >= tA && t < tD) ? C : 0;
    let y = simRK4(grid, engine.deriv, new Array(engine.size).fill(0), Cfun);
    if (noiseOn) y = y.map((v, k) => v + noiseSd * gauss() + drift * (grid[k] / tEnd));
    return y;                       // one time-series per concentration
  });

  // Single fixed region — placement used to be user-configurable
  // (region select + coordinate/radius UI); now a constant so the
  // stack-image compositing and export format stay unchanged.
  const region = { idx: 1, x: REGION_X, y: REGION_Y, r: REGION_R, traces };

  let gMin = Infinity, gMax = -Infinity;
  traces.forEach(tr => tr.forEach(v => {
    if (v < gMin) gMin = v;
    if (v > gMax) gMax = v;
  }));
  if (!isFinite(gMin)) { gMin = 0; gMax = 0; }

  const nSpots  = concs.length;
  const nFrames = grid.length;

  state.parsed   = { times: new Float64Array(grid), grid, nFrames, nSpots, concs,
                      globalMin: gMin, globalMax: gMax, regions: [region] };
  state.lastData = { grid };

  setStatus(`${nFrames} pts × ${nSpots} conc · signal ${gMax > 0 ? "on" : "black (Rmax = 0)"}`);

  // Was: updateStackImage() called directly. Firing an event instead means
  // this file doesn't need to import/know about the preview-rendering code
  // (or any future real-data loader that also calls notifyDataUpdated()).
  notifyDataUpdated();
}

function setModelVisibility() {
  const m = $("model").value;
  const groupMap = {
    simple:     m === "langmuir",
    twostate:   m === "twostate",
    bivAnalyte: m === "bivAnalyte",
    hetLigand:  m === "hetLigand",
  };
  document.querySelectorAll("[data-group]").forEach(el => {
    el.style.display = groupMap[el.dataset.group] ? "" : "none";
  });
  const hint = $("modelHint");
  if (hint) hint.textContent = MODEL_HINTS[m] ?? "";
}

function genDilution() {
  const [top, f, n] = ["dilTop","dilFactor","dilN"].map(id => +$(id).value);
  const pts = Array.from({ length: Math.max(1, Math.round(n)) }, (_, i) =>
    +(top / f ** i).toPrecision(4)
  );
  $("concSeries").value = pts.join(", ");
  simulate();
}

/* ══════════════════════════════════════════════════════════
   STACK IMAGE — composites the region's disk into each frame
   (Slated to move into its own preview.js — left here for now
   while main.js is the first split.)
   ══════════════════════════════════════════════════════════ */
   //to be moved over to the render.js file

export const IMG_W = 480, IMG_H = 640, MAX16 = 65535;

/* Fixed placement for the single simulated region — previously set via
   the region-select/place UI (#coordX, #coordY, #radius), now a constant. */
export const REGION_X = IMG_W / 2, REGION_Y = IMG_H / 2, REGION_R = 140;

function updateStackImage() {
  const results = $("results");
  if (!state.parsed || !state.parsed.regions || state.parsed.nSpots === 0 || state.parsed.nFrames === 0) {
    if (results) results.style.display = "none";
    return;
  }

  const totalFrames = state.parsed.nFrames * state.parsed.nSpots;
  const slider = $("frame-slider");
  if (slider) {
    slider.max = totalFrames - 1;
    if (+slider.value > totalFrames - 1) slider.value = 0;
  }

  const totalEl = $("total-frames");
  if (totalEl) totalEl.textContent =
    `${totalFrames}  (${state.parsed.nFrames} time points × ${state.parsed.nSpots} concentrations)`;
  if (results) results.style.display = "block";
  renderPreview(slider ? +slider.value : 0);
}

export function decodeFrame(globalFrame) {
  return {
    concIdx: Math.floor(globalFrame / state.parsed.nFrames),
    timeIdx: globalFrame % state.parsed.nFrames
  };
}

/* Brightness (0..MAX16) for one region at one (concIdx, timeIdx). */
function regionBrightness16(rg, concIdx, timeIdx) {
  if (concIdx >= rg.traces.length) return 0;
  const { globalMin, globalMax } = state.parsed;
  const denom = (globalMax - globalMin) || 1;
  const ru    = rg.traces[concIdx][timeIdx];
  const norm  = Math.max(0, (ru - globalMin) / denom);
  return Math.round(norm * MAX16);
}

/* Stamp a filled disk of value `val` centred at (cx, cy) with radius r. */
function stampDisk(mat, cx, cy, r, val) {
  if (!isFinite(cx) || !isFinite(cy) || !isFinite(r) || r <= 0) return;
  const x0 = Math.max(0, Math.floor(cx - r)), x1 = Math.min(IMG_W - 1, Math.ceil(cx + r));
  const y0 = Math.max(0, Math.floor(cy - r)), y1 = Math.min(IMG_H - 1, Math.ceil(cy + r));
  const r2 = r * r;
  for (let py = y0; py <= y1; py++) {
    for (let px = x0; px <= x1; px++) {
      const dx = px - cx, dy = py - cy;
      if (dx * dx + dy * dy <= r2) mat[py * IMG_W + px] = val;
    }
  }
}

/* Full IMG_H × IMG_W Uint16 frame: black background + the region's disk. */
export function getMatrix16(globalFrame) {
  const mat = new Uint16Array(IMG_H * IMG_W);   // zero = black
  if (!state.parsed || !state.parsed.regions) return mat;
  const { concIdx, timeIdx } = decodeFrame(globalFrame);
  for (const rg of state.parsed.regions) {
    const b16 = regionBrightness16(rg, concIdx, timeIdx);
    stampDisk(mat, rg.x, rg.y, rg.r, b16);       // disk sits on top of the noise
  }
  return mat;
}

/* Preview is built from the SAME matrix as the exported .stk frame. */
function renderPreview(globalFrame) {
  if (!state.parsed || !state.parsed.regions) return;
  const canvas = $("img-canvas");
  if (!canvas) return;
  canvas.width = IMG_W; canvas.height = IMG_H;

  const mat = getMatrix16(globalFrame);
  const ctx = canvas.getContext("2d");
  const img = ctx.createImageData(IMG_W, IMG_H);
  for (let i = 0; i < mat.length; i++) {
    const g = mat[i] >> 8;            // 16-bit -> 8-bit grey
    const j = i * 4;
    img.data[j] = g; img.data[j + 1] = g; img.data[j + 2] = g; img.data[j + 3] = 255;
  }
  ctx.putImageData(img, 0, 0);

  const { concIdx, timeIdx } = decodeFrame(globalFrame);
  const frameVal = $("frame-val");
  if (frameVal) frameVal.textContent = globalFrame;

  const concTag = $("conc-tag");
  if (concTag) {
    const cLabel = concIdx < state.parsed.concs.length ? fmtConc(state.parsed.concs[concIdx])
                                                 : `spot ${concIdx + 1}`;
    concTag.textContent = `${cLabel}  —  time point ${timeIdx} / ${state.parsed.nFrames - 1}`;
  }
}

export function findPeakInjectionFrame() {
  const tA = +$("tBase").value;
  const tD = tA + +$("tAssoc").value;
  const { regions, nFrames, nSpots, times } = state.parsed;

  let bestFrame = 0, bestRU = -Infinity;
  for (let f = 0; f < nFrames * nSpots; f++) {
    const { concIdx, timeIdx } = decodeFrame(f);
    const t = times[timeIdx];
    if (t < tA || t >= tD) continue;
    let ru = 0;
    for (const rg of regions)
      if (concIdx < rg.traces.length) ru = Math.max(ru, rg.traces[concIdx][timeIdx]);
    if (ru > bestRU) { bestRU = ru; bestFrame = f; }
  }
  return { frame: bestFrame, ru: bestRU };
}

/* ══════════════════════════════════════════════════════════
   EVENT LISTENERS
   ══════════════════════════════════════════════════════════ */

// Model/kinetics inputs: re-simulate directly from live DOM values.
MODEL_INPUT_IDS.forEach(id => {
  const el = $(id);
  if (!el) { console.warn("Missing element:", id); return; }
  el.addEventListener("input", () => {
    if (id === "model") setModelVisibility();
    simulate();
  });
});

on("mtlOn", "change", () => {
  const kt = $("ktField");
  if (kt) kt.style.display = $("mtlOn").checked ? "" : "none";
  simulate();
});

// Global inputs (shared by the whole image, incl. the concentration series).
["concSeries", "tBase", "tAssoc", "tDissoc", "noiseSd", "drift"]
  .forEach(id => on(id, "input", simulate));

on("noiseOn", "change", () => {
  const o = $("noiseOn").checked;
  const fields = $("noiseFields");
  if (fields) {
    fields.style.opacity       = o ? "1" : ".45";
    fields.style.pointerEvents = o ? "auto" : "none";
  }
  simulate();
});

on("genDil", "click", genDilution);
on("frame-slider", "input", function () { renderPreview(+this.value); });

// Whenever state.parsed/state.lastData is (re)populated — by simulate()
// here, or eventually by a real-data loader elsewhere — redraw the stack
// image. No direct import needed between the producer and this listener.
onDataUpdated(updateStackImage);

/* ── Init ────────────────────────────────────────────────── */
setModelVisibility();
simulate();