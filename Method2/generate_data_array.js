// Very basic version of a different way of generating two-dimentionsal data.
// This prototype does not incorporate user-interface - all parameters hard-coded
// Advantages of this method:
  // Integrates one time per model/set of parameters.
  // Adding any number of randomly generated circles / overlapping circles is fast and computationally inexpensive
  // Adding edge dominance is a simple extension
  // Incorporating drift in kinetic parameters is also reasonable, and can be done inexpensively

// Outline of what this does, assuming one dynamic present:
  // Step 1: Integrate one time (for whole time series, incorporating all concentrations) at RMax = 1. Set aside result Y
  // Step 2: Define a two dimensional, mxn scalar field \in [0,1]. Supplied parameters m and n
    // const s begins as an array of zeros
    // circle mask defines circular regions withing the field where elements equal 1
    // the following loop multiplies the corresponding capacity (this is where we will apply edge dominance and kinetic heterogeneity)
  // Step 3: The final result is the outer product of the scalar field and Y, result is a flat vector of length m*n*t which corresponds to the mxnxT array.

// Next steps
  // Step: Enable model selection.
    //  Note that bivalent analyte model does not allow for the Rmax factorization the same way, and that furthermore may be better represented by a truly spatial model. Must investigate.
  // Step: Random generation of circle placement, number of circles, and capacity.
    // This will consist furthermore of circles that overlap. MUST DECIDE: how to handle circle overlap? Likely not simply additive if this overlap represents vertical stacks of cellular growth.
  // Step: Stronger response on circle edges: This should be quick
  // Step: Heterogeneity of kinetic parameters within the same model/dynamic
    // Must investigate what distribution would be appropriate to use
  // Step: Heterogeneity of multiple models/dynamics within the same data



// Step 1
// assume, for simplicity one model - the Langmuir model
// assume that parameters have been chosen


const model = "langmuir";
const ka = 1e6, kd = 1e-3, Rmax = 1;
const str = "200, 100, 50, 25, 12.5, 6.25";
const tBase = 30, tAssoc = 120, tDiss = 300;

const vadd = (a,b)=>a.map((v,i)=>v+b[i]);
const vscale = (a,s)=>a.map(v=>v*s);
const vsum  = a=>a.reduce((x,y)=>x+y,0);

function simRK4(grid, deriv, y0, Cfun){
  const out=[]; let y=y0.slice(); out.push(vsum(y));
  for (let i=1;i<grid.length;i++){
    const t0=grid[i-1], h=grid[i]-grid[i-1];
    const C0=Cfun(t0), Cm=Cfun(t0+h/2), C1=Cfun(t0+h);
    const k1=deriv(y, C0);
    const k2=deriv(vadd(y,vscale(k1,h/2)), Cm);
    const k3=deriv(vadd(y,vscale(k2,h/2)), Cm);
    const k4=deriv(vadd(y,vscale(k3,h )), C1);
    y = vadd(y, vscale(vadd(vadd(k1,vscale(k2,2)), vadd(vscale(k3,2),k4)), h/6));
    out.push(vsum(y));
  }
  return out;
}

function makeDeriv(model, RmaxArg){
  return {
    size:1,
    deriv:(y,C) => { const R=y[0]; return [ka*C*(RmaxArg - R) - kd*R]; }
  };
}

function parseConcs(str){
  return str.split(/[\s,;]+/).map(s=>parseFloat(s)).filter(v=>Number.isFinite(v)&&v>0);
}

function simulate(){
  const concsNm = parseConcs(str);
  const cyc = tAssoc + tDiss;
  const N = concsNm.length;
  const total = tBase + N*cyc;
  const dt = 1; //Enforcing 1s resolution
  const npts = Math.round(total/dt); // Only distinct from total when dt /neq 1.
  const grid=[];
  for(let i=0;i <= npts; i++) grid.push(+(i*dt).toFixed(4));

  const concsM = concsNm.map(c=>c*1e-9);
  const Cfun = t=>{
    if (t < tBase) return 0;
    let k = Math.floor((t - tBase)/cyc);
    if (k >= N) k = N-1;
    return ((t - tBase) - k*cyc) < tAssoc ? concsM[k] : 0;
  };

  const base = makeDeriv(model, Rmax);
  const Y = simRK4(grid, base.deriv, new Array(base.size).fill(0), Cfun);
  return Y;
}

const Y = simulate();
console.log(Y.length, "points; final RU =", Y[Y.length-1]);

// Step 2: Two sets of hard-coded parameters here
const m = 640, n = 480;   // grid dimensions

// Stores circular region center, radius, and "capacity" as proportion of Rmax
const circles = [
  { ci: 10.5, cj: 20,   r: 5, capacity: 0.8 },
  { ci: 40,   cj: 15.5, r: 8, capacity: 0.5 },
  { ci: 25,   cj: 50,   r: 6, capacity: 1.0 },
];

// Uses flat indices INDEXED BY  i * n + j , m rows, n columns, row-major
function circleMask(m, n, ci, cj, r) {
  const mask = new Float64Array(m * n);   // flat, zeros
  const r2 = r * r;
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      const di = i - ci;
      const dj = j - cj;
      if (di*di + dj*dj <= r2) mask[i * n + j] = 1;   // flat index
    }
  }
  return mask;
}

const s = new Float64Array(m * n);

for (const circle of circles) {
  const mask = circleMask(m, n, circle.ci, circle.cj, circle.r);
  const c = circle.capacity;
  for (let k = 0; k < s.length; k++) {
    s[k] += mask[k] * c;
  }
}

// Step 3: Outer Product ALSO FLAT
// s: flat Float64Array of length m*n (the capacity field)
// Y: array of length T (unit-capacity response over time)
// result: flat Float64Array of length m*n*T

const T = Y.length;
const stack = new Float64Array(s.length * T);

for (let k = 0; k < s.length; k++) {   // each pixel
  const sk = s[k];
  for (let t = 0; t < T; t++) {        // each time point
    stack[k * T + t] = sk * Y[t];
  }
}


// Step 4: Noise

// noise parameters
const D          = 3.5 * Rmax;   // TOTAL accumulated drift (asymptote), not initial rate
const tau        = 500;          // settling time constant, s. Larger = more gradual.
const sigmaOU    = 0.02;         // OU step size — the hard-to-subtract wander
const thetaOU    = 0.005;        // OU mean-reversion rate (1/s). Small = long correlation
const decayOU    = true;         // scale OU jitter by the same exp(-t/tau) envelope
const sigmaPixel = 0.10 * Rmax;  // per-pixel noise: your 5-30% of Rmax range

// Box-Muller - i.e. make it Gaussian
function gauss(){
  let u = 0, v = 0;
  while (u === 0) u = Math.random();
  while (v === 0) v = Math.random();
  return Math.sqrt(-2*Math.log(u)) * Math.cos(2*Math.PI*v);
}

// --- drift vector: one value per frame, shared by ALL pixels ---
const dtGrid = 1;                       // must match the 1 s grid from Step 1
const driftCommon = new Float64Array(T);
let w = 0;                              // OU state — persists across t for temporal correlation
for (let t = 0; t < T; t++) {
  const envelope = decayOU ? Math.exp(-t*dtGrid / tau) : 1;
  // OU update: restoring pull toward 0, plus a random kick
  w += -thetaOU * w * dtGrid + sigmaOU * envelope * Math.sqrt(dtGrid) * gauss();
  driftCommon[t] = D * (1 - Math.exp(-t*dtGrid / tau)) + w;
}

// --- apply drift + noise to the stack ---
for (let k = 0; k < s.length; k++) {
  const sk = s[k];
  for (let t = 0; t < T; t++) {
    stack[k*T + t] = sk * Y[t]              // signal
    + driftCommon[t]          // SAME for every k — common-mode
    + sigmaPixel * gauss();   // fresh draw for every (k,t) — i.i.d.
  }
}



// Appendix
// Two index version
//function circleMask(m, n, ci, cj, r) {
//  const mask = Array.from({length: m}, () => new Array(n).fill(0));
//  const r2 = r * r;
//  for (let i = 0; i < m; i++) {
//    for (let j = 0; j < n; j++) {
//      const di = i - ci;
//      const dj = j - cj;
//      if (di*di + dj*dj <= r2) mask[i][j] = 1;
//    }
//  }
//  return mask;
//}
