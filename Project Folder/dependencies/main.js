'use strict';

/* ── Shared state ────────────────────────────────────────────
   Mutated (not reassigned) by whichever file produces data —
   simulate() today, potentially a real-data loader later. Anyone
   holding a reference to `state` always sees the latest values. */
export const state = { parsed: null, lastData: null };

/* ── DOM helpers ─────────────────────────────────────────── */
export const $  = id => document.getElementById(id);
export const on = (id, evt, fn) => { const el = $(id); if (el) el.addEventListener(evt, fn); };

export const setStatus = msg => { const el = $("status"); if (el) el.textContent = msg; };

/* ── Cross-file notification ─────────────────────────────────
   Fired whenever state.parsed / state.lastData has been (re)populated.
   Consumers (frame preview, export) subscribe via onDataUpdated()
   instead of the producer calling them directly — avoids a circular
   import between the file that computes data and the file(s) that
   react to it. */
export const DATA_UPDATED_EVENT = 'data-updated';
export const notifyDataUpdated = () => document.dispatchEvent(new CustomEvent(DATA_UPDATED_EVENT));
export const onDataUpdated = fn => document.addEventListener(DATA_UPDATED_EVENT, fn);