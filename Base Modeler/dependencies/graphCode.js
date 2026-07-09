/* ============================================================
   SPR Sensorgram Simulator — single-file reference tool
   ------------------------------------------------------------
   Three conceptual parts:
     (1) CONTROLS  — the HTML inputs above
     (2) MATH      — model functions below
     (3) GLUE      — read inputs, simulate, draw, repeat
   ============================================================ */

/* ---------- small vector helpers (for RK4 multi-state) ---------- */
const vadd = (a,b)=>a.map((v,i)=>v+b[i]);
const vscale = (a,s)=>a.map(v=>v*s);
const vsum  = a=>a.reduce((x,y)=>x+y,0);

/* ---------- (2) MATH: the binding models ---------- */

// 1:1 Langmuir — exact closed-form solution, evaluated on the grid.
function simLangmuir(grid, C, ka, kd, Rmax, tA, tD){
  const kobs = ka*C + kd;
  const Req  = ka*C*Rmax / (ka*C + kd);
  const Rd   = Req*(1 - Math.exp(-kobs*(tD - tA)));   // level at end of assoc.
  return grid.map(t=>{
    if (t < tA) return 0;
    if (t < tD) return Req*(1 - Math.exp(-kobs*(t - tA)));
    return Rd*Math.exp(-kd*(t - tD));
  });
}

function simHetLigand(grid,C,ka1,kd1,ka2,kd2,Rmax,Rmax2,tA,tD){
  const kobs1 = ka1*C + kd1;
  const kobs2 = ka2*C + kd2;
  const Req1  = ka1*C*Rmax / kobs1;
  const Req2  = ka2*C*Rmax2 / kobs2;
  const Rd1   = Req1*(1 - Math.exp(-kobs1*(tD - tA)));
  const Rd2   = Req2*(1 - Math.exp(-kobs2*(tD - tA)));
  return grid.map(t=>{
  if (t<tA) return 0;
  if (t<tD) return Req1*(1 - Math.exp(-kobs1*(t-tA))) + Req2*(1 - Math.exp(-kobs2*(t-tA)));
  return Rd1*Math.exp(-kd1*(t - tD)) + Rd2*Math.exp(-kd2*(t - tD));
  });
}

// Generic fixed-step RK4 integrator for state-space models.
// deriv(y, C) -> dy/dt array;  Cfun(t) -> analyte conc at time t.
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

// 1:1 with mass-transport limitation (quasi-steady surface concentration).
function simMassTransport(grid, C, ka, kd, Rmax, kt, Cfun){
  const deriv=(y,Cc)=>{
    const R=y[0];
    const Csurf=(kt*Cc + kd*R)/(kt + ka*(Rmax - R));
    return [ka*Csurf*(Rmax - R) - kd*R];
  };
  return simRK4(grid, deriv, [0], Cfun);
}

// Two-state conformational change:  A+B <-> AB <-> AB*
function simTwoState(grid, ka1,kd1,ka2,kd2, Rmax, Cfun){
  const deriv=(y,Cc)=>{
    const AB=y[0], ABs=y[1], free=Rmax - AB - ABs;
    return [
      ka1*Cc*free - kd1*AB - ka2*AB + kd2*ABs, // d[AB]/dt
      ka2*AB - kd2*ABs                          // d[AB*]/dt
    ];
  };
  return simRK4(grid, deriv, [0,0], Cfun);
}

function simBivAnalyte(grid, ka1,kd1,ka2,kd2, Rmax, Cfun){
  const deriv=(y,Cc)=>{
    const R1=y[0], R2=y[1], free=Rmax - R1 - 2*R2;
    return [
      2*ka1*Cc*free - kd1*R1 - ka2*R1*free + 2*kd2*R2, // dR1/dt
      2*ka2*R1*free - 2*kd2*R2                          // dR2/dt
    ];
  };
  return simRK4(grid, deriv, [0,0], Cfun);
}

// ---- Serial engine: express every model as an ODE derivative so the whole
// injection sequence can be integrated in ONE continuous pass, letting the
// bound state carry over between injections (no regeneration).
//
// Each model returns:
//   size    - number of state variables
//   deriv(y,C) - dy/dt, where C is the analyte conc the reaction actually sees
//   fluxCoef(y) - {a,b} describing the net rate at which analyte is drawn from
//                 SOLUTION, written as J = a*C + b (linear in C for every model
//                 here). Used only by the mass-transport modifier below. Note
//                 the bivalent model's crosslinking step consumes NO solution
//                 analyte, so only its first step contributes to {a,b}.
function makeDeriv(model, Rmax){
  if(model==="langmuir"){
    const ka=+$("ka").value, kd=+$("kd").value;
    return {size:1,
      deriv:(y,C)=>{ const R=y[0]; return [ka*C*(Rmax-R) - kd*R]; },
      fluxCoef:(y)=>({a: ka*(Rmax-y[0]), b: -kd*y[0]}) };
  }
  if(model==="hetLigand"){
    const ka1=+$("hetka1").value, kd1=+$("hetkd1").value,
          ka2=+$("hetka2").value, kd2=+$("hetkd2").value, Rmax2=+$("Rmax2").value;
    return {size:2,
      deriv:(y,C)=>{ const R1=y[0], R2=y[1];
        return [ka1*C*(Rmax-R1) - kd1*R1, ka2*C*(Rmax2-R2) - kd2*R2]; },
      // both sites draw from the same depleted surface pool
      fluxCoef:(y)=>({a: ka1*(Rmax-y[0]) + ka2*(Rmax2-y[1]),
                      b: -(kd1*y[0] + kd2*y[1])}) };
  }
  if(model==="bivAnalyte"){
    const ka1=+$("bivka1").value, kd1=+$("bivkd1").value,
          ka2=+$("bivka2").value, kd2=+$("bivkd2").value;
    return {size:2,
      deriv:(y,C)=>{ const R1=y[0], R2=y[1], free=Rmax-R1-2*R2;
        return [2*ka1*C*free - kd1*R1 - ka2*R1*free + 2*kd2*R2,
                2*ka2*R1*free - 2*kd2*R2]; },
      // only the first (solution-coupled) step draws analyte; crosslinking does not
      fluxCoef:(y)=>{ const free=Rmax-y[0]-2*y[1];
        return {a: 2*ka1*free, b: -kd1*y[0]}; } };
  }
  // two-state conformational change
  const ka1=+$("ka1").value, kd1=+$("kd1").value,
        ka2=+$("ka2").value, kd2=+$("kd2").value;
  return {size:2,
    deriv:(y,C)=>{ const AB=y[0], ABs=y[1], free=Rmax-AB-ABs;
      return [ka1*C*free - kd1*AB - ka2*AB + kd2*ABs, ka2*AB - kd2*ABs]; },
    // conformational interconversion does not touch solution; only step 1 does
    fluxCoef:(y)=>{ const free=Rmax-y[0]-y[1];
      return {a: ka1*free, b: -kd1*y[0]}; } };
}

// ---- Mass-transport modifier. Wraps ANY model's derivative. Under the
// quasi-steady two-compartment approximation, transport in balances net
// consumption out:  kt*(C - Cs) = J(Cs) = a*Cs + b.  Solving for the surface
// concentration Cs and feeding it to the reaction in place of bulk C:
//        Cs = (kt*C - b) / (kt + a).
// (For 1:1 this reproduces the classic Cs = (kt*C + kd*R)/(kt + ka*(Rmax-R)).)
function withTransport(base, kt){
  return {
    size: base.size,
    deriv:(y,C)=>{
      const {a,b} = base.fluxCoef(y);
      let Cs = (kt*C - b)/(kt + a);
      if(!isFinite(Cs) || Cs < 0) Cs = 0;   // guard against pathological inputs
      return base.deriv(y, Cs);
    }
  };
}

/* ================= Residual analysis: fit a 1:1 Langmuir ==================
   The "data" is one clean single injection from the CURRENT true model (with
   the transport modifier if enabled). We then fit a plain 1:1 Langmuir
   (ka, kd, Rmax) to it by Levenberg–Marquardt and show the residual. When the
   true model isn't 1:1, the best-fit residual carries the mechanism's
   fingerprint. Everything here is self-contained: no libraries. ============ */

// Solve the small dense linear system A x = b (n<=4) by Gaussian elimination
// with partial pivoting. Returns null if singular.
function solveLinear(A, b){
  const n=b.length, M=A.map((row,i)=>row.concat(b[i]));
  for(let c=0;c<n;c++){
    let piv=c; for(let r=c+1;r<n;r++) if(Math.abs(M[r][c])>Math.abs(M[piv][c])) piv=r;
    if(Math.abs(M[piv][c])<1e-300) return null;
    [M[c],M[piv]]=[M[piv],M[c]];
    for(let r=0;r<n;r++){ if(r===c) continue;
      const f=M[r][c]/M[c][c];
      for(let k=c;k<=n;k++) M[r][k]-=f*M[c][k];
    }
  }
  return M.map((row,i)=>row[n]/row[i]);
}

// 1:1 Langmuir response for arbitrary trial parameters, on the given grid.
function langmuirCurve(grid, Cfun, ka, kd, Rmax){
  const deriv=(y,C)=>[ka*C*(Rmax-y[0]) - kd*y[0]];
  return simRK4(grid, deriv, [0], Cfun);
}

// Levenberg–Marquardt fit of a 1:1 Langmuir to Ydata on (grid, Cfun).
// Parameters are fitted in log10 space (keeps them positive & well-scaled).
function fitLangmuir(grid, Cfun, Ydata){
  const model = p => langmuirCurve(grid, Cfun, 10**p[0], 10**p[1], 10**p[2]);
  const resid = Ymod => Ymod.map((v,i)=>v - Ydata[i]);
  const ssr   = r => r.reduce((s,v)=>s+v*v,0);

  // starting guess: generic rates, Rmax from the observed plateau
  const peak = Math.max(1e-6, ...Ydata);
  let p = [Math.log10(1e6), Math.log10(1e-2), Math.log10(peak)];
  let Ymod = model(p), r = resid(Ymod), S = ssr(r);

  let lambda = 1e-2;
  const n = grid.length, np = 3, dp = 1e-4;
  for(let iter=0; iter<80; iter++){
    // finite-difference Jacobian J[i][j] = d resid_i / d p_j
    const cols=[];
    for(let j=0;j<np;j++){
      const pj=p.slice(); pj[j]+=dp;
      const rj=resid(model(pj));
      const col=new Array(n);
      for(let i=0;i<n;i++) col[i]=(rj[i]-r[i])/dp;
      cols.push(col);
    }
    // JtJ (np×np) and Jtr (np)
    const JtJ=[[0,0,0],[0,0,0],[0,0,0]], Jtr=[0,0,0];
    for(let a=0;a<np;a++){
      for(let i=0;i<n;i++) Jtr[a]+=cols[a][i]*r[i];
      for(let b=a;b<np;b++){
        let s=0; for(let i=0;i<n;i++) s+=cols[a][i]*cols[b][i];
        JtJ[a][b]=s; JtJ[b][a]=s;
      }
    }
    // damped normal equations: (JtJ + λ·diag(JtJ)) δ = −Jtr
    const A=JtJ.map((row,a)=>row.map((v,b)=>a===b ? v+lambda*v : v));
    const delta=solveLinear(A, Jtr.map(v=>-v));
    if(!delta){ lambda=Math.min(lambda*4,1e12); continue; }
    const p2=p.map((v,i)=>v+delta[i]);
    const r2=resid(model(p2)), S2=ssr(r2);
    if(S2 < S){
      const rel=(S-S2)/S; p=p2; r=r2; S=S2;
      lambda=Math.max(lambda/3,1e-12);
      if(rel<1e-7) break;               // converged
    } else {
      lambda=Math.min(lambda*4,1e12);
      if(lambda>=1e12) break;           // stuck
    }
  }
  return { ka:10**p[0], kd:10**p[1], Rmax:10**p[2],
           Yfit:model(p), rmse:Math.sqrt(S/n) };
}

// Keep the concentration picker in sync with the series (default = highest).
function syncFitConc(concsNm){
  const sel=$("fitConc");
  const want=concsNm.map(c=>String(c));
  const have=Array.from(sel.options).map(o=>o.value);
  if(JSON.stringify(want)!==JSON.stringify(have)){
    const prev=sel.value;
    sel.innerHTML=concsNm.map(c=>`<option value="${c}">${c} nM</option>`).join("");
    sel.value = want.includes(prev) ? prev : String(Math.max(...concsNm));
  }
}

// Build one clean single-injection curve from the current true model, fit 1:1,
// and draw the fit overlay + residual. Called on every parameter change.
function updateFit(){
  const model=$("model").value, Rmax=parseFloat($("Rmax").value);
  const concsNm=parseConcs($("concSeries").value);
  if(!concsNm.length) return;
  syncFitConc(concsNm);
  const Cnm=parseFloat($("fitConc").value)||Math.max(...concsNm);
  const C=Cnm*1e-9;

  // single-injection timing & grid (lighter grid keeps the live fit snappy)
  const tBase=parseFloat($("tBase").value), tAssoc=parseFloat($("tAssoc").value),
        tDiss=parseFloat($("tDissoc").value), tEnd=tBase+tAssoc+tDiss;
  const npts=600, dt=tEnd/npts;
  const grid=[]; for(let i=0;i<=npts;i++) grid.push(+(i*dt).toFixed(4));
  const Cfun=t=>(t>=tBase && t<tBase+tAssoc)?C:0;

  // "data" = the true model (+ transport if on) at this one concentration
  const base=makeDeriv(model,Rmax);
  const useMTL=$("mtlOn").checked;
  const eng=useMTL ? withTransport(base, parseFloat($("kt").value)) : base;
  const Ydata=simRK4(grid, eng.deriv, new Array(eng.size).fill(0), Cfun);

  // best-fit 1:1 (never transport-limited — that's the "simple model")
  const fit=fitLangmuir(grid, Cfun, Ydata);
  const res=Ydata.map((v,i)=>v-fit.Yfit[i]);

  drawFit(grid, Ydata, fit.Yfit, tBase, tBase+tAssoc);
  drawResid(grid, res, tBase, tBase+tAssoc);

  const isExact = (model==="langmuir" && !useMTL);
  $("fitInfo").innerHTML =
    `best-fit 1:1 · k<sub>a</sub>=${fit.ka.toExponential(2)} · `+
    `k<sub>d</sub>=${fit.kd.toExponential(2)} · R<sub>max</sub>=${fit.Rmax.toFixed(1)} · `+
    `RMSE=${fit.rmse.toFixed(3)} RU` +
    (isExact ? " · (true model is 1:1 → residual is flat)"
             : " · residual structure reflects the true mechanism");
}

function drawFit(grid, Ydata, Yfit, tA, tD){
  const layout={
    margin:{l:60,r:14,t:10,b:20},
    paper_bgcolor:"rgba(0,0,0,0)", plot_bgcolor:"rgba(0,0,0,0)",
    font:{family:"IBM Plex Mono, monospace", size:11, color:"#1b1a17"},
    xaxis:{gridcolor:"#eae4d8", zeroline:false, linecolor:"#c9c2b4", matches:"x2"},
    yaxis:{title:{text:"Response (RU)",font:{size:12}}, gridcolor:"#eae4d8",
           zeroline:true, zerolinecolor:"#ddd6c8", linecolor:"#c9c2b4"},
    legend:{font:{size:10}, bgcolor:"rgba(255,253,248,.7)",
            bordercolor:"#ddd6c8", borderwidth:1, orientation:"h", x:0, y:1.12},
    shapes:[{type:"rect", xref:"x", yref:"paper", x0:tA, x1:tD, y0:0, y1:1,
             fillcolor:"rgba(15,107,102,.06)", line:{width:0}, layer:"below"}]
  };
  const traces=[
    {x:grid, y:Ydata, mode:"lines", type:"scatter", name:"data (true model)",
     line:{color:"#1b1a17", width:2}},
    {x:grid, y:Yfit, mode:"lines", type:"scatter", name:"1:1 fit",
     line:{color:"#b4541f", width:2, dash:"dot"}}
  ];
  Plotly.react("fitPlot", traces, layout,
    {responsive:true, displaylogo:false, displayModeBar:false});
  Plotly.Plots.resize("fitPlot");
}

function drawResid(grid, res, tA, tD){
  const layout={
    margin:{l:60,r:14,t:8,b:40},
    paper_bgcolor:"rgba(0,0,0,0)", plot_bgcolor:"rgba(0,0,0,0)",
    font:{family:"IBM Plex Mono, monospace", size:11, color:"#1b1a17"},
    xaxis:{title:{text:"Time (s)",font:{size:12}}, gridcolor:"#eae4d8",
           zeroline:false, linecolor:"#c9c2b4"},
    yaxis:{title:{text:"residual (RU)",font:{size:12}}, gridcolor:"#eae4d8",
           zeroline:true, zerolinecolor:"#0a4b47", linecolor:"#c9c2b4"},
    shapes:[{type:"rect", xref:"x", yref:"paper", x0:tA, x1:tD, y0:0, y1:1,
             fillcolor:"rgba(15,107,102,.06)", line:{width:0}, layer:"below"}]
  };
  const traces=[{x:grid, y:res, mode:"lines", type:"scatter", name:"residual",
     line:{color:"#0f6b66", width:1.6}}];
  Plotly.react("resid", traces, layout,
    {responsive:true, displaylogo:false, displayModeBar:false});
  Plotly.Plots.resize("resid");
}

/* ---------- viridis colour scale for the concentration family ---------- */
const VIRIDIS = ["#440154","#414487","#2a788e","#22a884","#7ad151","#fde725"];
function hex2rgb(h){return [parseInt(h.slice(1,3),16),parseInt(h.slice(3,5),16),parseInt(h.slice(5,7),16)];}
function viridis(x){
  x=Math.max(0,Math.min(1,x));
  const seg=x*(VIRIDIS.length-1), i=Math.floor(seg), f=seg-i;
  if(i>=VIRIDIS.length-1) return VIRIDIS[VIRIDIS.length-1];
  const a=hex2rgb(VIRIDIS[i]), b=hex2rgb(VIRIDIS[i+1]);
  const c=a.map((v,k)=>Math.round(v+(b[k]-v)*f));
  return `rgb(${c[0]},${c[1]},${c[2]})`;
}

/* ---------- Gaussian noise (Box–Muller) ---------- */
function gauss(){
  let u=0,v=0; while(u===0)u=Math.random(); while(v===0)v=Math.random();
  return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);
}

/* ---------- (3) GLUE: read inputs, simulate, draw ---------- */
const $ = id=>document.getElementById(id);

const MODEL_HINTS = {
  langmuir:"Simplest case: one analyte binding one immobilised ligand. Solved analytically.",
  twostate:"One type of receptor, but binding is followed by a conformational change that -locks- the complex.",
  hetLigand: "Two available binding sites with two completely independent dynamics.",
  bivAnalyte: "The analyte may, with sufficient density, bind two membrane receptors simultaneously."
};

let lastData = null; // cached for CSV export

function parseConcs(str){
  return str.split(/[\s,;]+/).map(s=>parseFloat(s)).filter(v=>Number.isFinite(v)&&v>0);
}

function fmtKD(kd, ka){
  const KD = kd/ka;                 // in molar
  const nM = KD*1e9;
  if (nM < 1)    return [(nM*1000).toPrecision(3), "pM"];
  if (nM < 1000) return [nM.toPrecision(3), "nM"];
  return [(nM/1000).toPrecision(3), "µM"];
}

function simulate(){
  const model = $("model").value;
  const Rmax  = parseFloat($("Rmax").value);
  const concsNm = parseConcs($("concSeries").value);   // injection order = list order

  const tBase  = parseFloat($("tBase").value);
  const tAssoc = parseFloat($("tAssoc").value);
  const tDiss  = parseFloat($("tDissoc").value);
  const cyc    = tAssoc + tDiss;                        // one injection's assoc+dissoc
  const N      = concsNm.length;
  const total  = tBase + N*cyc;                         // full serial timeline

  // One uniform grid over the WHOLE sequence. Keep the per-second resolution
  // identical to the old single-injection grid (1800 pts over one cycle).
  const singleEnd = tBase + cyc;
  const dt = singleEnd/1800;
  const npts = Math.max(1, Math.round(total/dt));
  const grid=[]; for(let i=0;i<=npts;i++) grid.push(+(i*dt).toFixed(4));

  // Staircase analyte profile: 0 during the initial baseline and during every
  // dissociation window; C_k during injection k's association window.
  const concsM = concsNm.map(c=>c*1e-9);               // nM -> M
  const Cfun = t=>{
    if (t < tBase) return 0;
    let k = Math.floor((t - tBase)/cyc);
    if (k >= N) k = N-1;                                // clamp at the far edge
    const inCyc = (t - tBase) - k*cyc;
    return inCyc < tAssoc ? concsM[k] : 0;
  };

  // Build the model, then optionally wrap it in the mass-transport modifier.
  const base = makeDeriv(model, Rmax);
  const useMTL = $("mtlOn").checked;
  const engine = useMTL ? withTransport(base, parseFloat($("kt").value)) : base;

  // Integrate ONCE. Bound state carries across injections automatically because
  // simRK4 never resets y — this is the "no regeneration" behaviour.
  let Y = simRK4(grid, engine.deriv, new Array(engine.size).fill(0), Cfun);

  // Detector noise + linear drift, applied to the single continuous vector.
  const noiseOn = $("noiseOn").checked;
  const noiseSd = parseFloat($("noiseSd").value)||0;
  const drift   = parseFloat($("drift").value)||0;
  if (noiseOn){
    Y = Y.map((v,i)=> v + noiseSd*gauss() + drift*(grid[i]/total));
  }

  // For the figure only: slice the one vector by injection and overlay each
  // segment on a shared local axis [0, cyc]. Each curve begins at the residual
  // left by the previous injection, so incomplete regeneration is visible.
  const traces=[];
  for (let k=0;k<N;k++){
    const start = tBase + k*cyc, end = start + cyc;
    const xs=[], ys=[];
    for (let i=0;i<grid.length;i++){
      if (grid[i] >= start-1e-9 && grid[i] <= end+1e-9){
        xs.push(+(grid[i]-start).toFixed(4));
        ys.push(Y[i]);
      }
    }
    const color = N>1 ? viridis(k/(N-1)) : viridis(0.35);
    const Cnm = concsNm[k];
    traces.push({
      x:xs, y:ys, mode:"lines", type:"scatter",
      name:(Cnm>=1?Cnm:Cnm.toPrecision(3))+" nM",
      line:{color:color, width:2}
    });
  }

  // The authoritative data is the single continuous vector (grid, Y).
  lastData = {serial:true, grid, Y, concsNm, tBase, tAssoc, tDiss, cyc, N};
  drawPlot(traces, 0, tAssoc);                          // shade injection window in local coords
  updateReadouts(model, N);
  updateFit();                                          // refit 1:1 & redraw residual (live)
  $("status").textContent = `${N} injections · serial, no regen · MTL ${useMTL?"on":"off"} · ${grid.length} pts · model: ${model}`;
}

function drawPlot(traces, tA, tD){
  const layout={
    margin:{l:64,r:18,t:14,b:52},
    paper_bgcolor:"rgba(0,0,0,0)",
    plot_bgcolor:"rgba(0,0,0,0)",
    font:{family:"IBM Plex Mono, monospace", size:12, color:"#1b1a17"},
    xaxis:{title:{text:"Time  (s)",font:{size:13}}, gridcolor:"#eae4d8",
           zeroline:false, linecolor:"#c9c2b4"},
    yaxis:{title:{text:"Response  (RU)",font:{size:13}}, gridcolor:"#eae4d8",
           zeroline:true, zerolinecolor:"#ddd6c8", linecolor:"#c9c2b4"},
    legend:{font:{size:11}, bgcolor:"rgba(255,253,248,.7)",
            bordercolor:"#ddd6c8", borderwidth:1},
    shapes:[
      // shade the association injection window
      {type:"rect", xref:"x", yref:"paper", x0:tA, x1:tD, y0:0, y1:1,
       fillcolor:"rgba(15,107,102,.06)", line:{width:0}, layer:"below"}
    ],
    annotations:[
      {x:(tA+tD)/2, y:1, xref:"x", yref:"paper", text:"injection",
       showarrow:false, font:{size:10,color:"#0a4b47"}, yanchor:"bottom"}
    ]
  };
  Plotly.react("plot", traces, layout,
    {responsive:true, displaylogo:false,
     modeBarButtonsToRemove:["select2d","lasso2d","autoScale2d"]});
  Plotly.Plots.resize("plot");
}

function updateReadouts(model, n){
  const box=$("readouts"); box.innerHTML="";
  const add=(k,v,accent)=>{
    const d=document.createElement("div");
    d.className="stat"+(accent?" accent":"");
    d.innerHTML=`<div class="k">${k}</div><div class="v">${v}</div>`;
    box.appendChild(d);
  };

  if (model==="twostate" || model==="bivAnalyte"){ // PLACEHOLDER ONLY for bivAnalyte
    const ka1=parseFloat($("ka1").value), kd1=parseFloat($("kd1").value),
          ka2=parseFloat($("ka2").value), kd2=parseFloat($("kd2").value);
    // apparent overall K_D for the two-state model
    const KD1 = kd1/ka1;                          // molar
    const appKD = KD1 * (kd2/(kd2+ka2));
    const nM = appKD*1e9;
    const [val,unit] = nM<1 ? [(nM*1000).toPrecision(3),"pM"]
                     : nM<1000 ? [nM.toPrecision(3),"nM"]
                     : [(nM/1000).toPrecision(3),"µM"];
    add("apparent K<sub>D</sub>", `${val} <small>${unit}</small>`, true);
    add("k<sub>d1</sub>", `${kd1} <small>s⁻¹</small>`);
  } else if (model==="hetLigand"){
    const ka1=parseFloat($("hetka1").value), kd1=parseFloat($("hetkd1").value),
          ka2=parseFloat($("hetka2").value), kd2=parseFloat($("hetkd2").value);
    // apparent overall K_Ds
    const [val1,unit1] = fmtKD(kd1,ka1);                          // molar
    const [val2,unit2] = fmtKD(kd2,ka2);

    add("K<sub>D1</sub>", `${val1} <small>${unit1}</small>`);
    add("K<sub>D2</sub>", `${val2} <small>${unit2}</small>`);
  } else {
    const ka=parseFloat($("ka").value), kd=parseFloat($("kd").value);
    const [val,unit]=fmtKD(kd,ka);
    add("K<sub>D</sub> = k<sub>d</sub>/k<sub>a</sub>", `${val} <small>${unit}</small>`, true);
    add("k<sub>a</sub>", `${ka.toExponential(1)} <small>M⁻¹s⁻¹</small>`);
    add("k<sub>d</sub>", `${kd.toExponential(1)} <small>s⁻¹</small>`);
  }
  add("Curves", n);
}

/* ---------- CSV / PNG export ---------- */
function exportCsv(){
  if(!lastData) return;
  const {grid,Y,concsNm,tBase,tAssoc,cyc,N}=lastData;
  // which concentration (nM) is being injected at time t; 0 during baseline/dissoc
  const injAt = t=>{
    if (t < tBase) return 0;
    let k = Math.floor((t - tBase)/cyc);
    if (k >= N) k = N-1;
    return (t - tBase - k*cyc) < tAssoc ? concsNm[k] : 0;
  };
  const rows=["time_s\tresponse_RU\tinjected_nM"];
  for(let i=0;i<grid.length;i++){
    rows.push(grid[i].toFixed(3)+"\t"+Y[i].toFixed(4)+"\t"+injAt(grid[i]));
  }
  const blob=new Blob([rows.join("\n")],{type:"text/tab-separated-values"});
  const a=document.createElement("a");
  a.href=URL.createObjectURL(blob);
  a.download="sensorgram_serial.tsv"; a.click();
  URL.revokeObjectURL(a.href);
}
function exportPng(){
  Plotly.downloadImage("plot",
    {format:"png", width:1200, height:760, filename:"sensorgram", scale:2});
}

/* ---------- wire up the UI ---------- */
function setModelVisibility(){
  const m=$("model").value;
  document.querySelectorAll("[data-group]").forEach(el=>{
    const g=el.getAttribute("data-group");
    let show=false;
    if(g==="simple") show=(m==="langmuir");
    if(g==="twostate") show=(m==="twostate");
    if(g==="hetLigand") show=(m==="hetLigand");
    if(g==="bivAnalyte") show=(m==="bivAnalyte");
    el.style.display= show ? "" : "none";
  });
  $("modelHint").textContent = MODEL_HINTS[m];
}

function genDilution(){
  const top=parseFloat($("dilTop").value), f=parseFloat($("dilFactor").value),
        n=Math.max(1,Math.round(parseFloat($("dilN").value)));
  const out=[]; let c=top;
  for(let i=0;i<n;i++){ out.push(+c.toPrecision(4)); c/=f; }
  $("concSeries").value=out.join(", ");
  simulate();
}

// attach listeners
["model","ka","kd","ka1","hetka1","kd1","hetkd1","ka2","hetka2","hetkd2","kd2","kt","Rmax","Rmax2","concSeries",
 "tBase","tAssoc","tDissoc","noiseSd","drift","bivka1","bivka2","bivkd1","bivkd2"].forEach(id=>{
  $(id).addEventListener("input",()=>{ if(id==="model") setModelVisibility(); simulate(); });
});
$("model").addEventListener("change",()=>{ setModelVisibility(); simulate(); });
$("noiseOn").addEventListener("change",()=>{
  const on=$("noiseOn").checked;
  $("noiseFields").style.opacity = on?"1":".45";
  $("noiseFields").style.pointerEvents = on?"auto":"none";
  simulate();
});
$("mtlOn").addEventListener("change",()=>{
  $("ktField").style.display = $("mtlOn").checked ? "" : "none";
  simulate();
});
$("fitConc").addEventListener("change", updateFit);
$("genDil").addEventListener("click", genDilution);
$("exportCsv").addEventListener("click", exportCsv);
$("exportPng").addEventListener("click", exportPng);

// initial render
setModelVisibility();
simulate();
