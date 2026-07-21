'use strict';

const DEFAULT_PALETTE = [
  [ 66, 165, 245], // blue
  [239,  83,  80], // red
  [102, 187, 106], // green
  [255, 202,  40], // amber
  [171,  71, 188], // purple
  [ 38, 198, 218], // cyan
  [255, 112,  67], // deep orange
  [141, 110,  99], // brown
];

export class PaintCanvas {
/**
   * @param {HTMLCanvasElement} canvas
   * @param {object} [opts]
   * @param {number} [opts.width=480]   internal pixel width
   * @param {number} [opts.height=640]  internal pixel height
   * @param {number} [opts.brushSize=15]
   * @param {number} [opts.regionCount=1]
   * @param {boolean}[opts.tintActive=true] highlight the active region
   * @param {(idx:number)=>number} [opts.valueProvider] 0..1 brightness per region
   * @param {number[][]} [opts.palette] per-region [r,g,b] colours
   */
    constructor(canvas, opts = {}) {
        if (!canvas) throw new Error('PaintCanvas: a <canvas> element is required');
        this.canvas = canvas;
        this.width  = opts.width  || canvas.width  || 480;
        this.height = opts.height || canvas.height || 640;
        canvas.width  = this.width;
        canvas.height = this.height;
        this.ctx = canvas.getContext('2d');

        this.masks        = new Map();                  // idx -> Uint8Array
        this.regionCount  = Math.max(1, opts.regionCount || 1);
        this.activeRegion = 1;
        this.brushSize    = Math.max(1, opts.brushSize || 15);
        this.mode         = 'paint';                    // 'paint' | 'erase'
        this.tintActive   = opts.tintActive !== false;
        this.palette      = opts.palette || DEFAULT_PALETTE;
        this.valueProvider = opts.valueProvider || null;

        this._listeners = {};
        this._painting  = false;
        this._last      = null;

        this._bindPointer();
        this.render();
    }

  /* ── Mask access ─────────────────────────────────────── */
    getMask(idx) {
        idx = +idx;
        if (!this.masks.has(idx)) this.masks.set(idx, new Uint8Array(this.width * this.height));
        return this.masks.get(idx);
    }

    setMask(idx, arr) {
        this.masks.set(+idx, arr instanceof Uint8Array ? arr : Uint8Array.from(arr));
        this.render();
    }

    clearRegion(idx) {
        this.masks.set(+idx, new Uint8Array(this.width * this.height));
        this.render();
        this._emit('change', { region: +idx });
    }

    clearAll() {
        this.masks.clear();
        this.render();
        this._emit('change', { region: null });
    }

  /* Plain-object export/import (JSON-friendly). */
    exportMasks() {
        const out = {};
        for (const [k, v] of this.masks) out[k] = Array.from(v);
        return out;
    }
    importMasks(obj) {
        this.masks.clear();
        for (const k in obj) this.masks.set(+k, Uint8Array.from(obj[k]));
        this.render();
        this._emit('change', { region: null });
    }

  /* ── Tool settings ───────────────────────────────────── */
    setActiveRegion(idx) { this.activeRegion = Math.max(1, +idx || 1); this.render(); }
    setBrushSize(px)     { this.brushSize    = Math.max(1, +px  || 1); }
    setMode(m)           { this.mode = (m === 'erase') ? 'erase' : 'paint'; }
    setValueProvider(fn) { this.valueProvider = fn || null; this.render(); }

    setRegionCount(n) {
        n = Math.max(1, n | 0);
        this.regionCount = n;
        for (const k of [...this.masks.keys()]) if (k > n) this.masks.delete(k);
        if (this.activeRegion > n) this.activeRegion = n;
        this.render();
    }

  /* ── Events ──────────────────────────────────────────── */
    on(evt, cb)  { (this._listeners[evt] || (this._listeners[evt] = [])).push(cb); return this; }
    off(evt, cb) { this._listeners[evt] = (this._listeners[evt] || []).filter(f => f !== cb); return this; }
    _emit(evt, payload) { (this._listeners[evt] || []).forEach(f => f(payload)); }

  /* ── Painting ────────────────────────────────────────── */
    _bindPointer() {
        const cv = this.canvas;
        cv.style.touchAction = 'none';
        cv.style.cursor = 'crosshair';

        cv.addEventListener('pointerdown', e => {
        this._painting = true; this._last = null;
        try { cv.setPointerCapture(e.pointerId); } catch (_) {}
        this._paintFromEvent(e);
        e.preventDefault();
        });
        cv.addEventListener('pointermove', e => {
        if (this._painting) this._paintFromEvent(e);
        });
        const end = () => { this._painting = false; this._last = null; };
        cv.addEventListener('pointerup', end);
        cv.addEventListener('pointercancel', end);
    }

    _eventToPixel(e) {
        const r = this.canvas.getBoundingClientRect();
        return {
        x: (e.clientX - r.left) * (this.width  / r.width),
        y: (e.clientY - r.top)  * (this.height / r.height)
        };
    }

    _paintFromEvent(e) {
        const { x, y } = this._eventToPixel(e);
        const r    = Math.max(1, this.brushSize);
        const val  = this.mode === 'erase' ? 0 : 1;
        const mask = this.getMask(this.activeRegion);

        if (this._last == null) {
        this._stamp(mask, x, y, r, val);
        } else {
        // interpolate along the drag so fast strokes leave no gaps
        const dx = x - this._last.x, dy = y - this._last.y;
        const dist = Math.hypot(dx, dy);
        const steps = Math.max(1, Math.ceil(dist / Math.max(1, r * 0.5)));
        for (let s = 0; s <= steps; s++) {
            const t = s / steps;
            this._stamp(mask, this._last.x + dx * t, this._last.y + dy * t, r, val);
        }
        }
        this._last = { x, y };
        this.render();
        this._emit('change', { region: this.activeRegion });
    }

  /* Filled brush disk, clamped to image bounds. */
    _stamp(mask, cx, cy, r, val) {
        const W = this.width, H = this.height;
        const x0 = Math.max(0, Math.floor(cx - r)), x1 = Math.min(W - 1, Math.ceil(cx + r));
        const y0 = Math.max(0, Math.floor(cy - r)), y1 = Math.min(H - 1, Math.ceil(cy + r));
        const r2 = r * r;
        for (let py = y0; py <= y1; py++) {
        for (let px = x0; px <= x1; px++) {
            const dx = px - cx, dy = py - cy;
            if (dx * dx + dy * dy <= r2) mask[py * W + px] = val;
        }
        }
    }

  /* Bounding box of a region's painted pixels, or null if empty. */
    maskBBox(idx) {
        const mask = this.masks.get(+idx);
        if (!mask) return null;
        const W = this.width, H = this.height;
        let minX = W, minY = H, maxX = -1, maxY = -1;
        for (let py = 0; py < H; py++) {
        for (let px = 0; px < W; px++) {
            if (mask[py * W + px]) {
            if (px < minX) minX = px; if (px > maxX) maxX = px;
            if (py < minY) minY = py; if (py > maxY) maxY = py;
            }
        }
        }
        return maxX < 0 ? null : { x0: minX, y0: minY, x1: maxX + 1, y1: maxY + 1 };
    }

  /* ── Rendering ───────────────────────────────────────── */
    render() {
        const W = this.width, H = this.height;
        const img = this.ctx.createImageData(W, H);
        const d = img.data;
        for (let i = 0; i < W * H; i++) { const j = i * 4; d[j] = d[j + 1] = d[j + 2] = 0; d[j + 3] = 255; }

        for (let idx = 1; idx <= this.regionCount; idx++) {
        const mask = this.masks.get(idx);
        if (!mask) continue;

        let col;
        if (this.valueProvider) {
            const v = Math.max(0, Math.min(1, this.valueProvider(idx) || 0));
            const g = Math.round(v * 255);
            col = [g, g, g];
        } else {
            col = this.palette[(idx - 1) % this.palette.length];
        }

        const active = (idx === this.activeRegion) && this.tintActive;
        let r = col[0], g = col[1], b = col[2];
        if (active) {                      // blend 45% toward a red highlight
            r = Math.round(r * 0.55 + 255 * 0.45);
            g = Math.round(g * 0.55 +  70 * 0.45);
            b = Math.round(b * 0.55 +  70 * 0.45);
        }

        for (let i = 0; i < mask.length; i++) {
            if (!mask[i]) continue;
            const j = i * 4;
            d[j] = r; d[j + 1] = g; d[j + 2] = b; d[j + 3] = 255;
        }
        }
        this.ctx.putImageData(img, 0, 0);
    }
}

/* Plain-script fallback: expose on window when not imported as a module. */
if (typeof window !== 'undefined') window.PaintCanvas = PaintCanvas;

export default PaintCanvas;