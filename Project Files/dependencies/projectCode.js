'use strict';

/* ── Helpers ─────────────────────────────────────────────── */
const $      = id => document.getElementById(id);
const on     = (id, evt, fn) => { const el = $(id); if (el) el.addEventListener(evt, fn); };
const setStatus = msg => { const el = $("status"); if (el) el.textContent = msg; };

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

const fmtConc = Cnm => (Cnm >= 1 ? Cnm : Cnm.toPrecision(3)) + " nM";

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

/* ── State ───────────────────────────────────────────────── */
let lastData = null;   // { grid }        — kept for timestamp helpers
let parsed   = null;   // { region, ... } — everything the image needs

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

  parsed   = { times: new Float64Array(grid), grid, nFrames, nSpots, concs,
               globalMin: gMin, globalMax: gMax, regions: [region] };
  lastData = { grid };

  setStatus(`${nFrames} pts × ${nSpots} conc · signal ${gMax > 0 ? "on" : "black (Rmax = 0)"}`);
  updateStackImage();
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
   ══════════════════════════════════════════════════════════ */

const IMG_W = 480, IMG_H = 640, MAX16 = 65535;

/* Fixed placement for the single simulated region — previously set via
   the region-select/place UI (#coordX, #coordY, #radius), now a constant. */
const REGION_X = IMG_W / 2, REGION_Y = IMG_H / 2, REGION_R = 140;

function updateStackImage() {
  const results = $("results");
  if (!parsed || !parsed.regions || parsed.nSpots === 0 || parsed.nFrames === 0) {
    if (results) results.style.display = "none";
    return;
  }

  const totalFrames = parsed.nFrames * parsed.nSpots;
  const slider = $("frame-slider");
  if (slider) {
    slider.max = totalFrames - 1;
    if (+slider.value > totalFrames - 1) slider.value = 0;
  }

  const totalEl = $("total-frames");
  if (totalEl) totalEl.textContent =
    `${totalFrames}  (${parsed.nFrames} time points × ${parsed.nSpots} concentrations)`;
  if (results) results.style.display = "block";
  renderPreview(slider ? +slider.value : 0);
}

function decodeFrame(globalFrame) {
  return {
    concIdx: Math.floor(globalFrame / parsed.nFrames),
    timeIdx: globalFrame % parsed.nFrames
  };
}

/* Brightness (0..MAX16) for one region at one (concIdx, timeIdx). */
function regionBrightness16(rg, concIdx, timeIdx) {
  if (concIdx >= rg.traces.length) return 0;
  const { globalMin, globalMax } = parsed;
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
function getMatrix16(globalFrame) {
  const mat = new Uint16Array(IMG_H * IMG_W);   // zero = black
  if (!parsed || !parsed.regions) return mat;
  const { concIdx, timeIdx } = decodeFrame(globalFrame);
  for (const rg of parsed.regions) {
    const b16 = regionBrightness16(rg, concIdx, timeIdx);
    stampDisk(mat, rg.x, rg.y, rg.r, b16);       // disk sits on top of the noise
  }
  return mat;
}

/* Preview is built from the SAME matrix as the exported .stk frame. */
function renderPreview(globalFrame) {
  if (!parsed || !parsed.regions) return;
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
    const cLabel = concIdx < parsed.concs.length ? fmtConc(parsed.concs[concIdx])
                                                 : `spot ${concIdx + 1}`;
    concTag.textContent = `${cLabel}  —  time point ${timeIdx} / ${parsed.nFrames - 1}`;
  }
}

function findPeakInjectionFrame() {
  const tA = +$("tBase").value;
  const tD = tA + +$("tAssoc").value;
  const { regions, nFrames, nSpots, times } = parsed;

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

/* ── base64 / raw-deflate helpers (unchanged) ────────────── */
function u8ToBase64(u8) {
  const CHUNK = 0x8000;
  let binary = "";
  for (let i = 0; i < u8.length; i += CHUNK) {
    binary += String.fromCharCode.apply(null, u8.subarray(i, i + CHUNK));
  }
  return btoa(binary);
}

async function deflateRawCompress(u8) {
  const cs = new CompressionStream('deflate-raw');
  const writer = cs.writable.getWriter();
  writer.write(u8);
  writer.close();
  const chunks = [];
  const reader = cs.readable.getReader();
  while (true) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
  }
  const totalLen = chunks.reduce((n, c) => n + c.length, 0);
  const out = new Uint8Array(totalLen);
  let offset = 0;
  for (const c of chunks) { out.set(c, offset); offset += c.length; }
  return out;
}

async function deflateRawDecompress(u8) {
  const ds = new DecompressionStream('deflate-raw');
  const writer = ds.writable.getWriter();
  writer.write(u8);
  writer.close();
  const chunks = [];
  const reader = ds.readable.getReader();
  while (true) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
  }
  const totalLen = chunks.reduce((n, c) => n + c.length, 0);
  const out = new Uint8Array(totalLen);
  let offset = 0;
  for (const c of chunks) { out.set(c, offset); offset += c.length; }
  return out;
}

async function encodeCompressedData(typedArrayOrBuffer) {
  const u8 = typedArrayOrBuffer instanceof ArrayBuffer
    ? new Uint8Array(typedArrayOrBuffer)
    : new Uint8Array(typedArrayOrBuffer.buffer, typedArrayOrBuffer.byteOffset, typedArrayOrBuffer.byteLength);
  const compressed = await deflateRawCompress(u8);
  return u8ToBase64(compressed);
}

async function decodeCompressedData(base64Str) {
  const binary = atob(base64Str);
  const compressed = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i++) compressed[i] = binary.charCodeAt(i);
  return deflateRawDecompress(compressed);
}

async function encodeTimeInput() {
  const total = parsed.nFrames * parsed.nSpots;
  const f32 = new Float32Array(total);
  for (let g = 0; g < total; g++) f32[g] = g + 10;
  const dataB64 = await encodeCompressedData(f32);
  return (
    `  <Input>\n` +
    `    <Name>Time</Name>\n` +
    `    <Data>${dataB64}</Data>\n` +
    `  </Input>`
  );
}

/* One <Input> (Roi{n}) per region — with a single region this is one entry. */
async function encodeResponseInput() {
  const { regions, nFrames, nSpots } = parsed;
  const parts = [];
  for (const rg of regions) {
    const f32 = new Float32Array(nFrames * nSpots);
    for (let concIdx = 0; concIdx < nSpots; concIdx++) {
      for (let timeIdx = 0; timeIdx < nFrames; timeIdx++) {
        const ru = concIdx < rg.traces.length ? rg.traces[concIdx][timeIdx] : 0;
        f32[concIdx * nFrames + timeIdx] = ru / 240;
      }
    }
    const dataB64 = await encodeCompressedData(f32);
    parts.push(
      `  <Input>\n` +
      `    <Name>Roi${rg.idx}</Name>\n` +
      `    <Data>${dataB64}</Data>\n` +
      `  </Input>`
    );
  }
  return parts.join("\n");
}

function formatTimestamp(when = new Date()) {
  const pad = n => String(n).padStart(2, '0');
  return `${pad(when.getMonth()+1)}/${pad(when.getDate())}/${when.getFullYear()} ` +
         `${pad(when.getHours())}:${pad(when.getMinutes())}:${pad(when.getSeconds())}`;
}

function stkTimeOffset(concIdx) {
  const g = lastData.grid, n = parsed.nFrames;
  const step = n > 1 ? (g[1] - g[0]) : 1;
  const span = g[n - 1] - g[0];
  return concIdx * (span + step);
}

function stkFileName(concIdx) {
  const c = parsed.concs[concIdx];
  const tag = (c != null) ? fmtConc(c).replace(/[^0-9A-Za-z.]+/g, '') : `spot${concIdx + 1}`;
  return `spr_stack_${tag}.stk`;
}

/* ── Builders ────────────────────────────────────────────── */

function buildStkBuffer(baseDate = new Date(), concIdx = 0) {
  const FRAME_TYPE_SPR_GRAY16 = 101;
  const nFrames       = parsed.nFrames;
  const bytesPerFrame = IMG_W * IMG_H * 2;
  const timeOffset    = stkTimeOffset(concIdx);

  const startTimeStr  = formatTimestamp(new Date(baseDate.getTime() + timeOffset * 1000));

  const enc = new TextEncoder();
  const c        = parsed.concs[concIdx];
  const concStr  = (c != null) ? fmtConc(c) : `spot ${concIdx + 1}`;
  const labelStr = 'SPR simulation';
  const descStr  = 'single-region simulation';

  const strBytes = s => enc.encode(s);
  const strSize  = s => 1 + strBytes(s).length;

  const headerSize =
    4 + strSize(startTimeStr) + strSize(concStr) + strSize(labelStr) + strSize(descStr) + 4 * 12;
  const frameHeaderSize = 4 + 4 + 4 + 4;
  const totalSize = headerSize + nFrames * (frameHeaderSize + bytesPerFrame);

  const buf  = new ArrayBuffer(totalSize);
  const view = new DataView(buf);
  const u8   = new Uint8Array(buf);
  let   pos  = 0;
  const LE   = true;

  const writeInt32   = v => { view.setInt32(pos, v, LE);   pos += 4; };
  const writeFloat32 = v => { view.setFloat32(pos, v, LE); pos += 4; };
  const writeString  = s => { const b = strBytes(s); u8[pos++] = b.length; u8.set(b, pos); pos += b.length; };

  // Header
  writeInt32(2);
  writeString(startTimeStr);
  writeString(concStr);
  writeString(labelStr);
  writeString(descStr);
  writeFloat32(0.0);   // flowRate
  writeFloat32(0.0);   // leadVolume
  writeFloat32(0.0);   // rinseVolume
  writeFloat32(0.0);   // exposureVolume
  writeFloat32(25.0);  // targetTemperature
  writeFloat32(0.0);   // angle
  writeFloat32(1.0);   // sprExposure
  writeFloat32(1.0);   // sprGain
  writeFloat32(0.0);   // bfExposure
  writeFloat32(1.0);   // bfRedGain
  writeFloat32(1.0);   // bfGreenGain
  writeFloat32(1.0);   // bfBlueGain

  // Frames for this concentration/spot only
  for (let timeIdx = 0; timeIdx < nFrames; timeIdx++) {
    const f         = concIdx * nFrames + timeIdx;
    const timestamp = timeOffset + lastData.grid[timeIdx];
    const mat       = getMatrix16(f);

    writeInt32(FRAME_TYPE_SPR_GRAY16);
    writeFloat32(timestamp);
    writeInt32(IMG_W);
    writeInt32(IMG_H);
    for (let i = 0; i < mat.length; i++) { view.setUint16(pos, mat[i], LE); pos += 2; }
  }
  return buf;
}

function addStkFilesToFolder(folder, startTimeStr = formatTimestamp()) {
  for (let c = 0; c < parsed.nSpots; c++) {
    folder.file(stkFileName(c), buildStkBuffer(startTimeStr, c));
  }
}

async function buildRoiXml(timestamp = formatTimestamp()) {
  const { frame: peakFrame } = findPeakInjectionFrame();

  const mat16   = getMatrix16(peakFrame);
  const grayBuf = new ArrayBuffer(mat16.length * 2);
  const grayDV  = new DataView(grayBuf);
  for (let i = 0; i < mat16.length; i++) grayDV.setUint16(i * 2, mat16[i], true);
  const grayCompressed = await deflateRawCompress(new Uint8Array(grayBuf));
  const grayB64  = u8ToBase64(grayCompressed);
  const sprGrayW = IMG_W, sprGrayH = IMG_H;

  const BF_W = sprGrayW * 2, BF_H = sprGrayH * 2;
  const bfBuf = new Uint8Array(BF_W * BF_H * 3).fill(128);
  const bfCompressed = await deflateRawCompress(bfBuf);
  const bfB64 = u8ToBase64(bfCompressed);

  // One ROI polygon per region (a square bounding box around each disk) —
  // with a single fixed region this always produces exactly one entry.
  const regions = parsed.regions;
  let roiEntries = "";
  regions.forEach(rg => {
    let x0, y0, x1, y1;
    if (isFinite(rg.x) && isFinite(rg.y) && isFinite(rg.r) && rg.r > 0) {
      x0 = Math.max(0, Math.round(rg.x - rg.r));
      y0 = Math.max(0, Math.round(rg.y - rg.r));
      x1 = Math.min(sprGrayW, Math.round(rg.x + rg.r));
      y1 = Math.min(sprGrayH, Math.round(rg.y + rg.r));
    } else {
      const w = Math.round(sprGrayW / 2), h = Math.round(sprGrayH / 2);
      x0 = Math.round((sprGrayW - w) / 2); y0 = Math.round((sprGrayH - h) / 2);
      x1 = x0 + w; y1 = y0 + h;
    }
    const poly = `${x0} ${y0} ${x0} ${y1} ${x1} ${y1} ${x1} ${y0}`;
    roiEntries +=
      `  <Roi>\n` +
      `    <Polygon>${poly}</Polygon>\n` +
      `    <Sensitivity>1</Sensitivity>\n` +
      `  </Roi>\n`;
  });

  const winW = Math.round(sprGrayW / 2), winH = Math.round(sprGrayH / 2);
  const winX0 = Math.round((sprGrayW - winW) / 2), winY0 = Math.round((sprGrayH - winH) / 2);
  const sprWindow = `${winX0}, ${winY0}, ${winX0 + winW}, ${winY0 + winH}`;

  return `<?xml version="1.0" encoding="utf-8"?>
          <RoiGroup>
            <Timestamp>${timestamp}</Timestamp>
            <SprWindow>${sprWindow}</SprWindow>
          ${roiEntries}
            <Snapshot>
              <Width>${sprGrayW}</Width>
              <Height>${sprGrayH}</Height>
              <Type>SprGray16</Type>
              <Data>${grayB64}</Data>
            </Snapshot>
            <Snapshot>
              <Width>${BF_W}</Width>
              <Height>${BF_H}</Height>
              <Type>BrightFieldRgb24</Type>
              <Data>${bfB64}</Data>
            </Snapshot>
          </RoiGroup>`;
}

async function buildBiXml(timestamp = formatTimestamp()) {
  const timeBlock  = await encodeTimeInput();
  const roiEntries = await encodeResponseInput();

  return `<?xml version="1.0" encoding="utf-8"?>
          <SPRm-Realtime>
            <Version>2.8.2</Version>
            <StartTime>${timestamp}</StartTime>
          ${roiEntries}
          ${timeBlock}
          </SPRm-Realtime>
          `;
}

async function downloadAll() {
  if (!parsed || !lastData) return;
  setStatus("Building export…");

  const baseDate  = new Date();
  const startTime = formatTimestamp(baseDate);

  const roiXml = await buildRoiXml(startTime);
  const biXml  = await buildBiXml(startTime);

  const zip = new JSZip();
  addStkFilesToFolder(zip.folder("DATA"), baseDate);
  zip.folder("ROI").file("spr.roi", roiXml);
  zip.folder("TIME").file("data.bi", biXml);

  const blob = await zip.generateAsync({ type: "blob" });
  const a = Object.assign(document.createElement("a"), {
    href:     URL.createObjectURL(blob),
    download: "spr_export.zip"
  });
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(a.href), 2000);

  setStatus("Export ready: spr_export.zip");
}

async function downloadSTK() {
  if (!parsed || !lastData) return;
  const zip = new JSZip();
  addStkFilesToFolder(zip.folder("DATA"), formatTimestamp());
  const blob = await zip.generateAsync({ type: "blob" });
  const a = Object.assign(document.createElement('a'), {
    href:     URL.createObjectURL(blob),
    download: 'DATA.zip'
  });
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(a.href), 2000);
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
on("downloadAll", "click", downloadAll);
on("frame-slider", "input", function () { renderPreview(+this.value); });

/* ── Init ────────────────────────────────────────────────── */
setModelVisibility();
simulate();