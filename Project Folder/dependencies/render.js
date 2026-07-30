'use strict';

import {
  $, on, state, onDataUpdated, fmtConc, IMG_W, IMG_H, MAX16
} from './main.js';
import { simulate } from './kinetics.js';

export function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export function generateCapacityField({
  m, n,
  capacityBins,
  targetConfluence,
  rMin, rMax,
  maxCircles = 5000,
  edgeFloor = 0.15,
  seed
}) {
  const rng = mulberry32(seed);

  const covered = new Uint8Array(m * n);
  let coveredCount = 0;

  const s = new Float32Array(m * n);
  const binField = new Int8Array(m * n).fill(-1);
  const circles = [];

  function placeOneCircle() {
    const ci = rng() * m;
    const cj = rng() * n;
    const r  = rMin + rng() * (rMax - rMin);
    const cap = capacityBins[(rng() * capacityBins.length) | 0];
    const binIdx = capacityBins.indexOf(cap);
    const r2 = r * r;
    const iLo = Math.max(0, Math.floor(ci - r)), iHi = Math.min(m - 1, Math.ceil(ci + r));
    const jLo = Math.max(0, Math.floor(cj - r)), jHi = Math.min(n - 1, Math.ceil(cj + r));
    for (let i = iLo; i <= iHi; i++) {
      for (let j = jLo; j <= jHi; j++) {
        const di = i - ci, dj = j - cj;
        const d2 = di * di + dj * dj;
        if (d2 <= r2) {
          const k = i * n + j;
          if (covered[k] === 0) { covered[k] = 1; coveredCount++; }
          // edge profile: floor at center (d=0) rising to 1.0 at rim (d=r)
          const profile = edgeFloor + (1 - edgeFloor) * (d2 / r2);
          const weighted = cap * profile;      // this cell's contribution here
          if (weighted > s[k]) {                // MAX composition on the WEIGHTED value
            s[k] = weighted;
            binField[k] = binIdx;               // winning cell sets the bin too
          }
        }
      }
    }
    circles.push({ ci, cj, r, capacity: cap, binIdx });
  }

  let iter = 0;
  while (coveredCount / (m * n) < targetConfluence && iter < maxCircles) {
    placeOneCircle();
    iter++;
  }

  return { s, binField, circles, achievedConfluence: coveredCount / (m * n) };
}

/* ══════════════════════════════════════════════════════════
   STACK IMAGE — composites the region's disk into each frame.
   Moved from kinetics.js. Reads state.parsed; never writes it.
   ══════════════════════════════════════════════════════════ */

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
  // Per-concentration local axis starts at each injection's own
  // association onset (t=0), not a shared baseline — see simulate()
  // in kinetics.js.
  const tA = 0;
  const tD = +$("tAssoc").value;
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
   EVENT LISTENERS — frame scrubbing + play/pause
   ══════════════════════════════════════════════════════════ */

on("frame-slider", "input", function () { renderPreview(+this.value); });

let playing = false, raf = null;

function stepFrame() {
  if (!playing || !state.parsed) { playing = false; return; }
  const slider = $("frame-slider");
  const totalFrames = state.parsed.nFrames * state.parsed.nSpots;
  if (!slider || totalFrames <= 0) { playing = false; return; }

  const t = (+slider.value + 4) % totalFrames;
  slider.value = t;
  renderPreview(t);
  raf = requestAnimationFrame(stepFrame);
}

on("play", "click", () => {
  const btn = $("Play");
  playing = !playing;
  if (btn) btn.textContent = playing ? "Pause" : "Play";
  if (playing) stepFrame(); else cancelAnimationFrame(raf);
});

onDataUpdated(updateStackImage);

simulate();