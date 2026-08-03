import { state, setStatus, $, fmtConc, IMG_W, IMG_H } from './main.js';
import { getMatrix16, decodeFrame, findPeakInjectionFrame } from './render.js';

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



function formatTimestamp(when = new Date()) {
  const pad = n => String(n).padStart(2, '0');
  return `${pad(when.getMonth()+1)}/${pad(when.getDate())}/${when.getFullYear()} ` +
         `${pad(when.getHours())}:${pad(when.getMinutes())}:${pad(when.getSeconds())}`;
}


async function encodeTimeInput() {
  const { nFrames, nSpots } = state.parsed;
  const grid = state.lastData.grid;
  const f32  = new Float32Array(nFrames * nSpots);

  for (let concIdx = 0; concIdx < nSpots; concIdx++) {
    const timeOffset = stkTimeOffset(concIdx);   // same offset the .stk uses
    for (let timeIdx = 0; timeIdx < nFrames; timeIdx++) {
      f32[concIdx * nFrames + timeIdx] = timeOffset + grid[timeIdx];
    }
  }

  const dataB64 = await encodeCompressedData(f32);
  return (
    `  <Input>\n` +
    `    <Name>Time</Name>\n` +
    `    <Data>${dataB64}</Data>\n` +
    `  </Input>`
  );
}

async function encodeResponseInput() {
  const { regions, nFrames, nSpots } = state.parsed;
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

function stkTimeOffset(concIdx) {
  const g = state.lastData.grid, n = state.parsed.nFrames;
  const step = n > 1 ? (g[1] - g[0]) : 1;
  const span = g[n - 1] - g[0];
  return concIdx * (span + step);
}

function stkFileName(concIdx) {
  const c = state.parsed.concs[concIdx];
  const tag = (c != null) ? fmtConc(c).replace(/[^0-9A-Za-z.]+/g, '') : `spot${concIdx + 1}`;
  return `spr_stack_${tag}.stk`;
}

/* Injection window read straight from the timing inputs: tB is the
   "baseline(s)" field (#tBase), tD is the "assoc.(s)" field (#tAssoc) —
   each offset by concIdx * "dissoc.(s)" (#tDissoc), so later concentrations
   in the serial run get pushed forward by one dissociation length per step.
   Used to place the .stk markers. */
function getInjectionWindow(concIdx = 0) {
  const tBase   = +$("tBase").value;
  const tAssoc  = +$("tAssoc").value;
  const tDissoc = +$("tDissoc").value;
  const offset  = concIdx * tDissoc;
  const tB = tBase + offset;
  const tD = tAssoc + offset;
  return { tB, tD };
}

function buildStkBuffer(baseDate = new Date(), concIdx = 0) {
  const FRAME_TYPE_SPR_GRAY16 = 101;
  const nFrames       = state.parsed.nFrames;
  const bytesPerFrame = IMG_W * IMG_H * 2;
  const timeOffset    = stkTimeOffset(concIdx);

  const startTimeStr  = formatTimestamp(new Date(baseDate.getTime() + timeOffset * 1000));

  const enc = new TextEncoder();
  const c        = state.parsed.concs[concIdx];
  const concStr  = (c != null) ? fmtConc(c) : `spot ${concIdx + 1}`;
  const labelStr = 'SPR simulation';
  const descStr  = `model: ${$('model').value}`;

  const strBytes = s => enc.encode(s);
  const strSize  = s => 1 + strBytes(s).length;

  const headerSize =
    4 + strSize(startTimeStr) + strSize(concStr) + strSize(labelStr) + strSize(descStr) + 4 * 12;
  const frameHeaderSize = 4 + 4 + 4 + 4;
  const markersSize = 6 * 4;   // trailing marker block written below (2 records: int32+float32+int32 each)
  const totalSize = headerSize + nFrames * (frameHeaderSize + bytesPerFrame) + markersSize;

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
    const timestamp = timeOffset + state.lastData.grid[timeIdx];
    const mat       = getMatrix16(f);

    writeInt32(FRAME_TYPE_SPR_GRAY16);
    writeFloat32(timestamp);
    writeInt32(IMG_W);
    writeInt32(IMG_H);
    for (let i = 0; i < mat.length; i++) { view.setUint16(pos, mat[i], LE); pos += 2; }
  }

  // Markers: tB from the "baseline(s)" input, tD from the "assoc.(s)" input,
  // both offset by concIdx * "dissoc.(s)".
  const { tB, tD } = getInjectionWindow(concIdx);
  writeInt32(1000);
  writeFloat32(tB);
  writeInt32(101);
  writeInt32(1000);
  writeFloat32(tD);
  writeInt32(102);

  return buf;
}

function addStkFilesToFolder(folder, startTimeStr = formatTimestamp()) {
  for (let c = 0; c < state.parsed.nSpots; c++) {
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

  // TEMPORARILY DISABLED — roiEntries left empty so the .roi file has no
  // <Roi> polygons. Re-enable by uncommenting this block.
  let roiEntries = "";
  // const regions = state.parsed.regions;
  // regions.forEach(rg => {
  //   let x0, y0, x1, y1;
  //   if (isFinite(rg.x) && isFinite(rg.y) && isFinite(rg.r) && rg.r > 0) {
  //     x0 = Math.max(0, Math.round(rg.x - rg.r));
  //     y0 = Math.max(0, Math.round(rg.y - rg.r));
  //     x1 = Math.min(sprGrayW, Math.round(rg.x + rg.r));
  //     y1 = Math.min(sprGrayH, Math.round(rg.y + rg.r));
  //   } else {
  //     const w = Math.round(sprGrayW / 2), h = Math.round(sprGrayH / 2);
  //     x0 = Math.round((sprGrayW - w) / 2); y0 = Math.round((sprGrayH - h) / 2);
  //     x1 = x0 + w; y1 = y0 + h;
  //   }
  //   const poly = `${x0} ${y0} ${x0} ${y1} ${x1} ${y1} ${x1} ${y0}`;
  //   roiEntries +=
  //     `  <Roi>\n` +
  //     `    <Polygon>${poly}</Polygon>\n` +
  //     `    <Sensitivity>1</Sensitivity>\n` +
  //     `  </Roi>\n`;
  // });

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
  // TEMPORARILY DISABLED — same as buildRoiXml's roiEntries. Re-enable by
  // restoring: const roiEntries = await encodeResponseInput();
  const roiEntries = "";

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
  if (!state.parsed || !state.lastData) return;
  setStatus("Building export…");

  const baseDate  = new Date();
  const startTime = formatTimestamp(baseDate);

  const roiXml = await buildRoiXml(startTime);
  const biXml  = await buildBiXml(startTime);

  const zip = new JSZip();
  addStkFilesToFolder(zip.folder("Stacks"), baseDate);
  zip.folder("ROI").file("spr.roi", roiXml);
  zip.folder("Data").file("data.bi", biXml);

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
  if (!state.parsed || !state.lastData) return;
  const zip = new JSZip();
  addStkFilesToFolder(zip.folder("Stack"), formatTimestamp());
  const blob = await zip.generateAsync({ type: "blob" });
  const a = Object.assign(document.createElement('a'), {
    href:     URL.createObjectURL(blob),
    download: 'DATA.zip'
  });
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(a.href), 2000);
}

document.getElementById("downloadAll").addEventListener("click", downloadAll);