'use strict';

import {
  $, on, state, onDataUpdated, fmtConc, IMG_W, IMG_H, MAX16
} from './main.js';
import { simulate } from './kinetics.js';

function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function generateCapacityField({
  m, n,
  capacityBins,
  targetConfluence,
  rMin, rMax,
  maxCircles = 5000,
  seed
}) {
  const rng = mulberry32(seed);
  const noOverlap = !allowOverlap();   // read once per call — see allowOverlap() below

  const covered = new Uint8Array(m * n);
  let coveredCount = 0;

  const s = new Float32Array(m * n);
  const circles = [];

  const MAX_OVERLAP_RETRIES = 30;

  function placeOneCircle() {
    let ci, cj, r, cap, tries = 0;
    do {
      ci = rng() * m;
      cj = rng() * n;
      r  = rMin + rng() * (rMax - rMin);
      cap = capacityBins[(rng() * capacityBins.length) | 0];
      tries++;
    } while (noOverlap && circleOverlapsAny(ci, cj, r, circles) && tries < MAX_OVERLAP_RETRIES);

    if (noOverlap && circleOverlapsAny(ci, cj, r, circles)) return;   // couldn't fit one — skip

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
          if (cap > s[k]) s[k] = cap;   // MAX composition, flat value within the disk
        }
      }
    }
    circles.push({ ci, cj, r, capacity: cap });
  }

  let iter = 0;
  while (coveredCount / (m * n) < targetConfluence && iter < maxCircles) {
    placeOneCircle();
    iter++;
  }

  return { s, circles, achievedConfluence: coveredCount / (m * n) };
}

function circleOverlapsAny(ci, cj, r, circles) {
  for (const c of circles) {
    const dx = ci - c.ci, dy = cj - c.cj;
    const minDist = r + c.r;
    if (dx * dx + dy * dy < minDist * minDist) return true;
  }
  return false;
}

function allowOverlap() {
  const el = $("overlap");
  return el ? el.checked : true;
}

on("overlap", "change", () => refreshCapacityFieldIfNeeded());

function getSeed() {
  const el = $("inputSeed");
  const v = el ? +el.value : NaN;
  return Number.isFinite(v) ? v : 0;
}

on("genSeed", "click", () => {
  const el = $("inputSeed");
  if (el) el.value = Math.floor(Math.random() * 2 ** 32);
  refreshCapacityFieldIfNeeded();
});

on("inputSeed", "input", refreshCapacityFieldIfNeeded);

const CAPACITY_BINS = [0.4, 0.6, 0.8, 1.0];
const CIRCLE_R_MIN = 15, CIRCLE_R_MAX = 25;

let lastFieldSeed, lastFieldConfluence, lastFieldOverlap;

function getTargetConfluence() {
  const el = $("Confluency");
  const pct = el ? +el.value : NaN;
  return (Number.isFinite(pct) ? pct : 0) / 100;
}

function refreshCapacityFieldIfNeeded() {
  const seed = getSeed();
  const targetConfluence = getTargetConfluence();
  const overlap = allowOverlap();

  const unchanged = state.capacityField
    && seed === lastFieldSeed
    && targetConfluence === lastFieldConfluence
    && overlap === lastFieldOverlap;
  if (unchanged) return;

  state.capacityField = generateCapacityField({
    m: IMG_H, n: IMG_W,
    capacityBins: CAPACITY_BINS,
    targetConfluence,
    rMin: CIRCLE_R_MIN, rMax: CIRCLE_R_MAX,
    seed
  });

  lastFieldSeed = seed;
  lastFieldConfluence = targetConfluence;
  lastFieldOverlap = overlap;
}

on("Confluency", "input", refreshCapacityFieldIfNeeded);

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

function regionBrightness16(rg, concIdx, timeIdx) {
  if (concIdx >= rg.traces.length) return 0;
  const { globalMin, globalMax } = state.parsed;
  const denom = (globalMax - globalMin) || 1;
  const ru    = rg.traces[concIdx][timeIdx];
  const norm  = Math.max(0, (ru - globalMin) / denom);
  return Math.round(norm * MAX16);
}

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

function stampCapacityField(mat, rg, b16) {
  const field = state.capacityField;
  if (!field) return;
  const { s } = field;
  const cx = rg.x, cy = rg.y, r = rg.r;
  if (!isFinite(cx) || !isFinite(cy) || !isFinite(r) || r <= 0) return;
  const x0 = Math.max(0, Math.floor(cx - r)), x1 = Math.min(IMG_W - 1, Math.ceil(cx + r));
  const y0 = Math.max(0, Math.floor(cy - r)), y1 = Math.min(IMG_H - 1, Math.ceil(cy + r));
  const r2 = r * r;
  for (let py = y0; py <= y1; py++) {
    for (let px = x0; px <= x1; px++) {
      const dx = px - cx, dy = py - cy;
      if (dx * dx + dy * dy <= r2) {
        const k = py * IMG_W + px;
        let v = Math.round(s[k] * b16);
        if (v < 0) v = 0; else if (v > MAX16) v = MAX16;
        mat[k] = v;
      }
    }
  }
}

export function getMatrix16(globalFrame) {
  const mat = new Uint16Array(IMG_H * IMG_W);   // zero = black
  if (!state.parsed || !state.parsed.regions) return mat;
  const { concIdx, timeIdx } = decodeFrame(globalFrame);
  for (const rg of state.parsed.regions) {
    const b16 = regionBrightness16(rg, concIdx, timeIdx);
    stampCapacityField(mat, rg, b16);
  }
  return mat;
}

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
  const btn = $("play");
  playing = !playing;
  if (btn) btn.textContent = playing ? "Pause" : "Play";
  if (playing) stepFrame(); else cancelAnimationFrame(raf);
});

onDataUpdated(updateStackImage);
simulate();
refreshCapacityFieldIfNeeded();