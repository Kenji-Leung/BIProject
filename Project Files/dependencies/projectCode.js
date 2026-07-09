/* ═══════════════════════════════════════════════════════════
   SPR Kinetic Simulator + Stack Image Tool  (no graphs)
   ───────────────────────────────────────────────────────────
   Reads the control inputs → simulates the binding model →
   drives the stack-image preview and the .stk/.roi/.bi export.

   Simulation engine: every model is expressed as an ODE derivative
   (makeDeriv) so the optional mass-transport modifier (withTransport)
   can wrap ANY of them. Each concentration is integrated
   independently, which is what the stack image needs (one trace per
   concentration). No plotting of any kind.

   DOM access for graph-era elements (#status, #readouts, ...) is
   guarded so the tool never crashes if those elements are absent.
   ═══════════════════════════════════════════════════════════ */

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
   Each model returns:
     size       - number of state variables
     deriv(y,C) - dy/dt, C = analyte conc the reaction actually sees
     fluxCoef(y)- {a,b} for the net rate analyte is drawn from SOLUTION
                  (J = a*C + b), used only by the mass-transport modifier. */
function makeDeriv(model, Rmax) {
  if (model === "langmuir") {
    const ka = +$("ka").value, kd = +$("kd").value;
    return {
      size: 1,
      deriv: (y, C) => { const R = y[0]; return [ka * C * (Rmax - R) - kd * R]; },
      fluxCoef: (y) => ({ a: ka * (Rmax - y[0]), b: -kd * y[0] })
    };
  }
  if (model === "hetLigand") {
    const ka1 = +$("hetka1").value, kd1 = +$("hetkd1").value,
          ka2 = +$("hetka2").value, kd2 = +$("hetkd2").value, Rmax2 = +$("Rmax2").value;
    return {
      size: 2,
      deriv: (y, C) => { const R1 = y[0], R2 = y[1];
        return [ka1 * C * (Rmax - R1) - kd1 * R1, ka2 * C * (Rmax2 - R2) - kd2 * R2]; },
      fluxCoef: (y) => ({ a: ka1 * (Rmax - y[0]) + ka2 * (Rmax2 - y[1]),
                          b: -(kd1 * y[0] + kd2 * y[1]) })
    };
  }
  if (model === "bivAnalyte") {
    const ka1 = +$("bivka1").value, kd1 = +$("bivkd1").value,
          ka2 = +$("bivka2").value, kd2 = +$("bivkd2").value;
    return {
      size: 2,
      deriv: (y, C) => { const R1 = y[0], R2 = y[1], free = Rmax - R1 - 2 * R2;
        return [2 * ka1 * C * free - kd1 * R1 - ka2 * R1 * free + 2 * kd2 * R2,
                2 * ka2 * R1 * free - 2 * kd2 * R2]; },
      // only the first (solution-coupled) step draws analyte; crosslinking does not
      fluxCoef: (y) => { const free = Rmax - y[0] - 2 * y[1];
        return { a: 2 * ka1 * free, b: -kd1 * y[0] }; }
    };
  }
  // two-state conformational change:  A + B <-> AB <-> AB*
  const ka1 = +$("ka1").value, kd1 = +$("kd1").value,
        ka2 = +$("ka2").value, kd2 = +$("kd2").value;
  return {
    size: 2,
    deriv: (y, C) => { const AB = y[0], ABs = y[1], free = Rmax - AB - ABs;
      return [ka1 * C * free - kd1 * AB - ka2 * AB + kd2 * ABs, ka2 * AB - kd2 * ABs]; },
    fluxCoef: (y) => { const free = Rmax - y[0] - y[1];
      return { a: ka1 * free, b: -kd1 * y[0] }; }
  };
}

/* ── Mass-transport modifier. Wraps ANY model's derivative. Under the
   quasi-steady two-compartment approximation:  kt*(C - Cs) = a*Cs + b,
   so the surface concentration Cs = (kt*C - b)/(kt + a) is fed to the
   reaction in place of bulk C. ─────────────────────────────────────── */
function withTransport(base, kt) {
  return {
    size: base.size,
    deriv: (y, C) => {
      const { a, b } = base.fluxCoef(y);
      let Cs = (kt * C - b) / (kt + a);
      if (!isFinite(Cs) || Cs < 0) Cs = 0;   // guard against pathological inputs
      return base.deriv(y, Cs);
    }
  };
}

/* ── Formatting ──────────────────────────────────────────── */
const parseConcs = str =>
  str.split(/[\s,;]+/).map(Number).filter(v => Number.isFinite(v) && v > 0);

/* ── Model metadata ──────────────────────────────────────── */
const MODEL_HINTS = {
  langmuir:   "Simplest case: one analyte binding one immobilised ligand.",
  hetLigand:  "Two available binding sites with two completely independent dynamics.",
  twostate:   "Binding followed by a conformational change that locks the complex — note the slow dissociation.",
  bivAnalyte: "The analyte may, at sufficient density, bind two membrane receptors simultaneously."
};

/* ── State ───────────────────────────────────────────────── */
let lastData     = null;
let parsed       = null;
let orderedSpots = [];

function simulate() {
  const model   = $("model").value;
  const Rmax    = +$("Rmax").value;
  const concs   = parseConcs($("concSeries").value);
  const tA      = +$("tBase").value;
  const tD      = tA + +$("tAssoc").value;
  const tEnd    = tD + +$("tDissoc").value;
  const noiseOn = $("noiseOn").checked;
  const noiseSd = +$("noiseSd").value || 0;
  const drift   = +$("drift").value   || 0;

  // One frame per second, so nFrames = grid.length drives the stack image.
  const grid = Array.from({ length: Math.round(tEnd) + 1 }, (_, i) => i);

  // Build the model once, optionally wrapped in the mass-transport modifier.
  const base   = makeDeriv(model, Rmax);
  const useMTL = $("mtlOn").checked;
  const engine = useMTL ? withTransport(base, +$("kt").value) : base;

  // Integrate each concentration independently → one trace per concentration.
  const traces = concs.map(Cnm => {
    const C    = Cnm * 1e-9;
    const Cfun = t => (t >= tA && t < tD) ? C : 0;
    let y = simRK4(grid, engine.deriv, new Array(engine.size).fill(0), Cfun);
    if (noiseOn) y = y.map((v, i) => v + noiseSd * gauss() + drift * (grid[i] / tEnd));
    const label = (Cnm >= 1 ? Cnm : Cnm.toPrecision(3)) + " nM";
    return { x: grid, y, name: label };
  });

  lastData = { grid, traces, concs };
  setStatus(`${concs.length} curves · ${grid.length} pts · model: ${model} · MTL ${useMTL ? "on" : "off"}`);
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
   STACK IMAGE — driven entirely by lastData, no file upload
   ══════════════════════════════════════════════════════════ */

const IMG_W = 480, IMG_H = 640, MAX16 = 65535;

const parseConc  = label => { const m = label.match(/[\d.]+/); return m ? +m[0] : 0; };
const sortSpots  = headers =>
  headers.map((h, i) => ({ i, conc: parseConc(h), label: h }))
         .sort((a, b) => a.conc - b.conc);

function updateStackImage() {
  const results = $("results");
  // Nothing to show without at least one concentration curve.
  if (!lastData || !lastData.traces.length) {
    if (results) results.style.display = "none";
    parsed = null; orderedSpots = [];
    return;
  }

  const { grid, traces } = lastData;
  const nFrames = grid.length;
  const nSpots  = traces.length;

  const raw = new Float64Array(nFrames * nSpots);
  traces.forEach((tr, c) => tr.y.forEach((v, r) => { raw[r * nSpots + c] = v; }));

  let globalMin = Infinity, globalMax = -Infinity;
  for (const v of raw) {
    if (v < globalMin) globalMin = v;
    if (v > globalMax) globalMax = v;
  }

  parsed = {
    times: new Float64Array(grid),
    raw, nFrames, nSpots, globalMin, globalMax,
    headers: traces.map(t => t.name)
  };
  orderedSpots = sortSpots(parsed.headers);

  const totalFrames = nFrames * nSpots;
  const slider = $("frame-slider");
  if (slider) {
    slider.max = totalFrames - 1;
    if (+slider.value > totalFrames - 1) slider.value = 0;
  }

  const totalEl = $("total-frames");
  if (totalEl) totalEl.textContent =
    `${totalFrames}  (${nFrames} time points × ${nSpots} concentrations)`;
  if (results) results.style.display = "block";
  renderPreview(slider ? +slider.value : 0);
}

function decodeFrame(globalFrame) {
  return {
    concIdx: Math.floor(globalFrame / parsed.nFrames),
    timeIdx: globalFrame % parsed.nFrames
  };
}

function getBrightness(globalFrame) {
  const { raw, globalMin, globalMax, nSpots } = parsed;
  const { concIdx, timeIdx } = decodeFrame(globalFrame);
  const ru   = raw[timeIdx * nSpots + orderedSpots[concIdx].i];
  const norm = Math.max(0, (ru - globalMin) / (globalMax - globalMin || 1));
  return { brightness16: Math.round(norm * MAX16), norm };
}

/* Returns a full IMG_H × IMG_W Uint16Array for the given frame.
   Currently uniform; swap fill logic here when per-pixel values diverge. */
function getMatrix16(globalFrame) {
  const { brightness16 } = getBrightness(globalFrame);
  return new Uint16Array(IMG_H * IMG_W).fill(brightness16);
}

function renderPreview(globalFrame) {
  if (!parsed || !orderedSpots.length) return;
  const canvas = $("img-canvas");
  if (!canvas) return;
  canvas.width = IMG_W; canvas.height = IMG_H;
  const { norm } = getBrightness(globalFrame);
  const g   = Math.round(norm * 255);
  const ctx = canvas.getContext("2d");
  ctx.fillStyle = `rgb(${g},${g},${g})`;
  ctx.fillRect(0, 0, IMG_W, IMG_H);

  const { concIdx, timeIdx } = decodeFrame(globalFrame);
  const frameVal = $("frame-val");
  if (frameVal) frameVal.textContent = globalFrame;
  const concTag = $("conc-tag");
  if (concTag) concTag.textContent =
    `${orderedSpots[concIdx].label}  —  time point ${timeIdx} / ${parsed.nFrames - 1}`;
}

function findPeakInjectionFrame() {
  // tA / tD come straight from the timing inputs, mirroring simulate().
  const tA = +$("tBase").value;
  const tD = tA + +$("tAssoc").value;

  let bestFrame = 0, bestRU = -Infinity;
  for (let f = 0; f < parsed.nFrames * parsed.nSpots; f++) {
    const { concIdx, timeIdx } = decodeFrame(f);
    const t = parsed.times[timeIdx];
    if (t < tA || t >= tD) continue; // only the injection/association window
    const ru = parsed.raw[timeIdx * parsed.nSpots + orderedSpots[concIdx].i];
    if (ru > bestRU) { bestRU = ru; bestFrame = f; }
  }
  return { frame: bestFrame, ru: bestRU };
}

function u8ToBase64(u8) {
  // chunked to avoid call-stack / argument-count limits on large arrays
  const CHUNK = 0x8000;
  let binary = "";
  for (let i = 0; i < u8.length; i += CHUNK) {
    binary += String.fromCharCode.apply(null, u8.subarray(i, i + CHUNK));
  }
  return btoa(binary);
}

/* Compresses a Uint8Array using raw DEFLATE (no zlib/gzip header).
   Confirmed via Python's zlib.decompress(data, -15) that this is the
   exact format the reference .roi file's SprGray16 Data field uses —
   feeding raw/uncompressed bytes into that field is what was causing
   generated .roi files to fail to load silently. Requires
   CompressionStream('deflate-raw'), supported in all current browsers
   (Chrome/Edge/Firefox/Safari). */
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

/* Decompresses a Uint8Array using raw INFLATE (the inverse of
   deflateRawCompress), via DecompressionStream('deflate-raw'). */
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

/** Encode a typed array (or plain byte buffer) into a base64 string the
 *  same way the reference software's Encode() does: raw-deflate, base64. */
async function encodeCompressedData(typedArrayOrBuffer) {
  const u8 = typedArrayOrBuffer instanceof ArrayBuffer
    ? new Uint8Array(typedArrayOrBuffer)
    : new Uint8Array(typedArrayOrBuffer.buffer, typedArrayOrBuffer.byteOffset, typedArrayOrBuffer.byteLength);
  const compressed = await deflateRawCompress(u8);
  return u8ToBase64(compressed);
}

/** Decode a base64 string produced by encodeCompressedData() back into
 *  raw bytes. Wrap the result in the appropriate typed array to read it. */
async function decodeCompressedData(base64Str) {
  const binary = atob(base64Str);
  const compressed = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i++) compressed[i] = binary.charCodeAt(i);
  return deflateRawDecompress(compressed);
}

async function encodeTimeInput() {
  const frames = parsed.nFrames;
  const concentration = parsed.nSpots;
  const total = frames * concentration;
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

async function encodeResponseInput() {
  const { raw, nFrames, nSpots } = parsed;

  const f32 = new Float32Array(nFrames * nSpots);
  for (let concIdx = 0; concIdx < nSpots; concIdx++) {
    const col = orderedSpots[concIdx].i;
    for (let timeIdx = 0; timeIdx < nFrames; timeIdx++) {
      const ru = raw[timeIdx * nSpots + col];
      f32[concIdx * nFrames + timeIdx] = (ru / 240);
    }
  }
  const dataB64 = await encodeCompressedData(f32);
  return (
    `  <Input>\n` +
    `    <Name>Roi1</Name>\n` +
    `    <Data>${dataB64}</Data>\n` +
    `  </Input>`
  );
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

/* Filesystem-safe name for a concentration's .stk file, e.g. "50 nM" -> spr_stack_50nM.stk */
function stkFileName(concIdx) {
  const tag = orderedSpots[concIdx].label.replace(/[^0-9A-Za-z.]+/g, '');
  return `spr_stack_${tag || concIdx}.stk`;
}

/* ── Builders: return content only, no download ──────────── */

/* Builds ONE .stk file holding just the frames for a single
   concentration (concIdx into the sorted orderedSpots). Timestamps
   are continuous across concentrations via stkTimeOffset(). */
function buildStkBuffer(baseDate = new Date(), concIdx = 0) {
  const FRAME_TYPE_SPR_GRAY16 = 101;
  const nFrames       = parsed.nFrames;
  const bytesPerFrame = IMG_W * IMG_H * 2;
  const timeOffset    = stkTimeOffset(concIdx);

  const startTimeStr  = formatTimestamp(new Date(baseDate.getTime() + timeOffset * 1000));

  const enc = new TextEncoder();
  const concStr  = orderedSpots[concIdx].label;   // this concentration only
  const labelStr = 'SPR simulation';
  const descStr  = `model: ${$('model').value}`;

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

  // Frames for this concentration only
  for (let timeIdx = 0; timeIdx < nFrames; timeIdx++) {
    const f         = concIdx * nFrames + timeIdx;      // global frame index
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

/* Adds one .stk per concentration into the given JSZip folder. */
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

  const winW = Math.round(sprGrayW / 2), winH = Math.round(sprGrayH / 2);
  const winX0 = Math.round((sprGrayW - winW) / 2), winY0 = Math.round((sprGrayH - winH) / 2);
  const winX1 = winX0 + winW,                      winY1 = winY0 + winH;
  const sprWindow = `${winX0}, ${winY0}, ${winX1}, ${winY1}`;

  const topPolygon = `${winX0} ${winY0} ${winX0} ${winY1} ${winX1} ${winY1} ${winX1} ${winY0}`;
  const roiEntries =
    `  <Roi>\n` +
    `    <Polygon>${topPolygon}</Polygon>\n` +
    `    <Sensitivity>1</Sensitivity>\n` +
    `  </Roi>\n`;

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
  addStkFilesToFolder(zip.folder("DATA"), formatTimestamp());   // one .stk per concentration
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

["model","ka","kd","ka1","kd1","ka2","kd2",
 "hetka1","hetkd1","hetka2","hetkd2",
 "bivka1","bivkd1","bivka2","bivkd2",
 "kt","Rmax","Rmax2","concSeries","tBase","tAssoc","tDissoc","noiseSd","drift"]
  .forEach(id => {
    const el = $(id);
    if (!el) { console.warn("Missing element:", id); return; }
    el.addEventListener("input", () => {
      if (id === "model") setModelVisibility();
      simulate();
    });
  });

on("model", "change", () => { setModelVisibility(); simulate(); });

on("noiseOn", "change", () => {
  const on = $("noiseOn").checked;
  const fields = $("noiseFields");
  if (fields) {
    fields.style.opacity       = on ? "1" : ".45";
    fields.style.pointerEvents = on ? "auto" : "none";
  }
  simulate();
});

on("mtlOn", "change", () => {
  const kt = $("ktField");
  if (kt) kt.style.display = $("mtlOn").checked ? "" : "none";
  simulate();
});

on("genDil", "click", genDilution);
on("downloadAll", "click", downloadAll);
on("frame-slider", "input", function () { renderPreview(+this.value); });

/* ── Init ────────────────────────────────────────────────── */
setModelVisibility();
simulate();