document.addEventListener("DOMContentLoaded", () => {

'use strict';
const G=9.81;
const TARGET_H=0.43;
const STORAGE='sbls-wt-v4';

const defs={
  d:{s:'d-s',n:'d-n',dv:'d-dv',fmt:v=>v.toFixed(2)+' m',sv:v=>v/10,vs:v=>Math.round(v*10),nmGet:v=>v.toFixed(2),nmSet:v=>v},
  h:{s:'h-s',n:'h-n',dv:'h-dv',fmt:v=>v.toFixed(2)+' m',sv:v=>v/100,vs:v=>Math.round(v*100),nmGet:v=>v.toFixed(2),nmSet:v=>v},
  k:{s:'k-s',n:'k-n',dv:'k-dv',fmt:v=>v.toFixed(0)+' N/m',sv:v=>+v,vs:v=>v,nmGet:v=>v.toFixed(0),nmSet:v=>v},
  x:{s:'x-s',n:'x-n',dv:'x-dv',fmt:v=>v.toFixed(2)+' m',sv:v=>v/100,vs:v=>Math.round(v*100),nmGet:v=>v.toFixed(2),nmSet:v=>v},
  M:{s:'M-s',n:'M-n',dv:'M-dv',fmt:v=>(v*1e3).toFixed(0)+' g',sv:v=>v/1e3,vs:v=>Math.round(v*1e3),nmGet:v=>(v*1e3).toFixed(0),nmSet:v=>v/1e3},
  m:{s:'m-s',n:'m-n',dv:'m-dv',fmt:v=>(v*1e3).toFixed(0)+' g',sv:v=>v/1e3,vs:v=>Math.round(v*1e3),nmGet:v=>(v*1e3).toFixed(0),nmSet:v=>v/1e3},
  e:{s:'e-s',n:'e-n',dv:'e-dv',fmt:v=>v.toFixed(2),sv:v=>v/100,vs:v=>Math.round(v*100),nmGet:v=>v.toFixed(2),nmSet:v=>v},
  e_f:{s:'e_f-s',n:'e_f-n',dv:'e_f-dv',fmt:v=>(v*100).toFixed(1)+'%',sv:v=>v/100,vs:v=>Math.round(v*100),nmGet:v=>(v*100).toFixed(1),nmSet:v=>v/100},
  D_o:{fmt:v=>(v*1e3).toFixed(0)+' mm',value:0.12},
  D_b:{fmt:v=>(v*1e3).toFixed(0)+' mm',value:0.04},
  ang_min:{s:'ang_min-s',n:'ang_min-n',dv:'ang_min-dv',fmt:v=>v.toFixed(1)+'°',sv:v=>+v,vs:v=>v,nmGet:v=>v.toFixed(1),nmSet:v=>v},
  ef_mult:{s:'ef_mult-s',n:'ef_mult-n',dv:null,
    fmt:v=>v>=0?'×'+(1+v).toFixed(2):'×'+(1+v).toFixed(2),
    sv:v=>v/100,vs:v=>Math.round(v*100),
    nmGet:v=>v.toFixed(2),nmSet:v=>v}
};

function gv(){
  const r={};
  for(const[k,d]of Object.entries(defs)){
    if(d.value !== undefined){
      r[k] = d.value;
      continue;
    }
    const el = document.getElementById(d.s);
    if(!el) continue;
    r[k] = d.sv(+el.value);
  }
  return r;
}

function updateEfTag(param){
  const tag=document.getElementById('ef_mult-tag');
  if(!tag)return;
  const mult=1+param;
  tag.textContent='×'+mult.toFixed(2);
  tag.className='ef-tag '+(param>0.001?'pos':param<-0.001?'neg':'zero');
}

function calcShot(deg, p){
  const th = deg * Math.PI / 180;
  const {k, x, M, m, e, h} = p;

  // Spring energy → plunger speed
  const KE = 0.5 * k * x * x - M * G * x * Math.sin(th);
  if(KE <= 0) return null;

  const vp = Math.sqrt(2 * KE / M);
  // Inelastic collision: plunger hits ball
  const vb = M * (1 + e) * vp / (M + m);

  const ux = vb * Math.cos(th);
  const uy = vb * Math.sin(th);

  // Solve y(t) = h + uy*t - 0.5*G*t² = TARGET_H
  // 0.5*G*t² - uy*t + (TARGET_H - h) = 0
  const a_ = 0.5 * G;
  const b_ = -uy;
  const c_ = TARGET_H - h;          // NOTE: corrected sign vs. original
  const disc = b_*b_ - 4*a_*c_;
  if(disc < 0) return null;

  const sqrtDisc = Math.sqrt(disc);
  const t1 = (-b_ - sqrtDisc) / (2 * a_);
  const t2 = (-b_ + sqrtDisc) / (2 * a_);

  // vy at each root
  const vy1 = uy - G * t1;
  const vy2 = uy - G * t2;

  // Prefer the root where ball is DESCENDING (vy < 0) and t > 0
  let tFinal = null;
  if(t2 > 1e-9 && vy2 <= 0) {
    tFinal = t2;                     // almost always the correct one
  } else if(t1 > 1e-9 && vy1 <= 0) {
    tFinal = t1;
  } else {
    // Launch is below TARGET_H: ball crosses on the way up — accept the
    // positive root even if vy > 0 (no descending crossing exists)
    const candidates = [t1, t2].filter(t => t > 1e-9);
    if(!candidates.length) return null;
    tFinal = Math.max(...candidates);
  }

  const xAtTarget = ux * tFinal;
  const vy_at     = uy - G * tFinal;
  const descentDeg = Math.atan2(-vy_at, ux) * 180 / Math.PI;

  // Time of flight to ground (y=0) — for display only
  const tof = (uy + Math.sqrt(uy*uy + 2*G*h)) / G;

  // Peak
  const tp   = uy / G;
  const ymax = tp > 0 ? h + uy*tp - 0.5*G*tp*tp : h;

  return {xAtTarget, descentDeg, vb, tof, ymax,
          ux, uy, x0:0, y0:h, tTarget:tFinal, vy_at};
}

function findAngle(targetDist, p){
  let prev = null;
  for(let i = 1; i <= 895; i++){
    const deg  = i * 0.1;
    const shot = calcShot(deg, p);
    if(!shot){ prev = null; continue; }
    if(prev !== null){
      const dA = prev.xAtTarget - targetDist;
      const dB = shot.xAtTarget - targetDist;
      if(dA * dB <= 0){
        const sol = Math.abs(dB) <= Math.abs(dA)
          ? {theta:deg,  ...shot}
          : {theta:deg-0.1, ...prev};
        if(sol.descentDeg >= p.ang_min) return sol;
      }
    }
    prev = shot; prev._deg = deg;
  }
  return null;
}

function pathPoints(theta, p, steps=300){
  const shot = calcShot(theta, p);
  if(!shot) return [];

  const {ux, uy, y0, tTarget} = shot;
  const pts = [];

  for(let i = 0; i <= steps; i++){
    const t  = (i / steps) * tTarget;
    const px = ux * t;
    const py = y0 + uy * t - 0.5 * G * t * t;
    pts.push({x: px, y: py});
  }
  // Snap last point exactly to target height
  pts[pts.length - 1] = {x: shot.xAtTarget, y: TARGET_H};
  return pts;
}

function saveStorage(){
  try{
    const r={};
    for(const[k,d]of Object.entries(defs)){
      const sl=document.getElementById(d.s);
      if(sl) r[k]=sl.value;
    }
    localStorage.setItem(STORAGE,JSON.stringify(r));
  }catch(_){}
}

function loadStorage(){
  try{
    const raw=JSON.parse(localStorage.getItem(STORAGE)||'null');
    if(!raw) return;
    for(const[k,d]of Object.entries(defs)){
      if(raw[k]==null) continue;
      const sl=document.getElementById(d.s);
      const nm=document.getElementById(d.n);
      if(!sl||!nm) continue;
      sl.value=raw[k];
      const param=d.sv(+raw[k]);
      nm.value=d.nmGet(param);
      if(d.dv){const dv=document.getElementById(d.dv);if(dv)dv.textContent=d.fmt(param);}
      if(k==='ef_mult') updateEfTag(param);
    }
  }catch(_){}
}

function renderResults(solNom, solMin, solMax, solShift, p){
  const row      = document.getElementById('results-row');
  const statusEl = document.getElementById('status-txt');
  const angleEl  = document.getElementById('calc-angle');
  const Esp = 0.5 * p.k * p.x * p.x;

  if(p.D_b >= p.D_o){
    angleEl.textContent='--'; statusEl.textContent='Ball too large';
    row.innerHTML=`<div class="alert-wrap"><div class="alert-box"><div class="alert-icon">!</div><div class="alert-body"><h4>ลูกใหญ่กว่าช่องเปิด</h4><p>⌀ลูก <b>${(p.D_b*1e3).toFixed(0)} mm</b> ≥ ⌀ช่อง <b>${(p.D_o*1e3).toFixed(0)} mm</b></p></div></div></div>`;
    return;
  }
  if(!solNom){
    angleEl.textContent='--'; statusEl.textContent='No solution';
    row.innerHTML=`<div class="alert-wrap"><div class="alert-box"><div class="alert-icon">!</div><div class="alert-body"><h4>ไม่พบมุมยิงที่เหมาะสม</h4><p>E<sub>spring</sub> ${Esp.toFixed(3)} J ไม่เพียงพอ หรือ descent &lt; ${p.ang_min.toFixed(1)}°<br>ลอง ↑ k · ↑ x · ↑ h · ↓ ระยะ</p></div></div></div>`;
    return;
  }

  angleEl.textContent = solNom.theta.toFixed(1)+'°';
  const ef_disp = p.ef_mult!==0 ? ` · ef×${(1+p.ef_mult).toFixed(2)}` : '';
  statusEl.textContent = `θ=${solNom.theta.toFixed(1)}° · descent ${solNom.descentDeg.toFixed(1)}°${ef_disp}`;

  const hMin=!!solMin, hMax=!!solMax;
  let bandStr = hMin&&hMax
    ? `${solMin.theta.toFixed(1)}° — ${solMax.theta.toFixed(1)}°`
    : hMin ? `≥ ${solMin.theta.toFixed(1)}°`
    : hMax ? `≤ ${solMax.theta.toFixed(1)}°`
    : `${solNom.theta.toFixed(1)}° (±0)`;
  const descOk = solNom.descentDeg >= p.ang_min;

  row.innerHTML=`
  <div class="rcard nom">
    <div class="wm">${solNom.theta.toFixed(0)}°</div>
    <div class="rbadge nom">Launch Angle</div>
    <div class="rangle">${solNom.theta.toFixed(1)}<span class="unit">°</span></div>
    <div class="rband">Band <b>${bandStr}</b> &nbsp;(±${(p.e_f*100).toFixed(1)}%)</div>
    <div class="rmetrics">
      <div class="rm"><div class="rm-lbl">Ball v₀</div><div class="rm-val">${solNom.vb.toFixed(2)}<span class="rm-unit">m/s</span></div></div>
      <div class="rm"><div class="rm-lbl">Air time</div><div class="rm-val">${solNom.tof.toFixed(3)}<span class="rm-unit">s</span></div></div>
      <div class="rm"><div class="rm-lbl">Peak h</div><div class="rm-val">${solNom.ymax.toFixed(2)}<span class="rm-unit">m</span></div></div>
      <div class="rm"><div class="rm-lbl">Descent</div><div class="rm-val ${descOk?'ok':'warn'}">${solNom.descentDeg.toFixed(1)}<span class="rm-unit">°</span></div></div>
    </div>
  </div>
  <div class="rcard geo">
    <div class="wm geo-wm">GEO</div>
    <div class="rbadge geo">Fixed Geometry</div>
    <div class="geo-spec-row">
      <div class="geo-spec-item">
        <div class="geo-spec-val">430<span class="geo-spec-unit">mm</span></div>
        <div class="geo-spec-lbl">Target H</div>
      </div>
      <div class="geo-spec-sep"></div>
      <div class="geo-spec-item">
        <div class="geo-spec-val">${(p.D_o*1e3).toFixed(0)}<span class="geo-spec-unit">mm</span></div>
        <div class="geo-spec-lbl">Opening ⌀</div>
      </div>
      <div class="geo-spec-sep"></div>
      <div class="geo-spec-item">
        <div class="geo-spec-val">${(p.D_b*1e3).toFixed(0)}<span class="geo-spec-unit">mm</span></div>
        <div class="geo-spec-lbl">Ball ⌀</div>
      </div>
      <div class="geo-spec-sep"></div>
      <div class="geo-spec-item">
        <div class="geo-spec-val ${(p.D_o-p.D_b)*1e3 <= 20 ? 'warn' : ''}">${((p.D_o-p.D_b)*1e3).toFixed(0)}<span class="geo-spec-unit">mm</span></div>
        <div class="geo-spec-lbl">Clearance</div>
      </div>
    </div>
    <div class="rmetrics" style="border-top:1px solid var(--stroke);padding-top:10px;margin-top:2px">
      <div class="rm">
        <div class="rm-lbl">Spring E</div>
        <div class="rm-val">${Esp.toFixed(3)}<span class="rm-unit">J</span></div>
      </div>
      <div class="rm">
        <div class="rm-lbl">ef→ θ</div>
        <div class="rm-val" style="color:${p.ef_mult!==0?'var(--tiger)':'var(--ink-light)'}">${solShift&&p.ef_mult!==0?solShift.theta.toFixed(1)+'°':'—'}</div>
      </div>
      <div class="rm">
        <div class="rm-lbl">ef_mult</div>
        <div class="rm-val" style="color:${p.ef_mult!==0?'var(--tiger)':'var(--ink-light)'}">×${(1+p.ef_mult).toFixed(2)}</div>
      </div>
      <div class="rm">
        <div class="rm-lbl">Gap ratio</div>
        <div class="rm-val">${((p.D_o-p.D_b)/p.D_o*100).toFixed(0)}<span class="rm-unit">%</span></div>
      </div>
    </div>
  </div>`;
}

function drawTraj(solNom, solMin, solMax, solShift, p){
  const cv   = document.getElementById('traj');
  const wrap = cv.parentElement;
  const W    = wrap.clientWidth;
  const H    = wrap.clientHeight;
  cv.width   = W * devicePixelRatio;
  cv.height  = H * devicePixelRatio;
  cv.style.width  = W + 'px';
  cv.style.height = H + 'px';
  const ctx = cv.getContext('2d');
  ctx.scale(devicePixelRatio, devicePixelRatio);
  ctx.clearRect(0, 0, W, H);

  const pad = {l:52, r:24, t:48, b:38};
  const pw = W - pad.l - pad.r;
  const ph = H - pad.t - pad.b;
  if(pw <= 0 || ph <= 0) return;

  function peakOf(sol){
    if(!sol) return 0;
    const tp = sol.uy / G;
    return tp > 0 ? sol.y0 + sol.uy*tp - 0.5*G*tp*tp : sol.y0;
  }

  let maxX = Math.max(p.d * 1.35, 3);
  let maxY = Math.max(p.h * 1.4, TARGET_H * 1.5, 1);
  [solNom, solMin, solMax, solShift].forEach(s=>{
    if(s) maxY = Math.max(maxY, peakOf(s) * 1.18);
  });
  maxY = Math.max(maxY, p.h * 1.15);

  const toX = rx => pad.l + (rx / maxX) * pw;
  const toY = ry => pad.t + ph - (ry / maxY) * ph;

  // Grid
  ctx.strokeStyle='#f0ede8'; ctx.lineWidth=1;
  for(let i=0;i<=6;i++){const cx=pad.l+pw*i/6;ctx.beginPath();ctx.moveTo(cx,pad.t);ctx.lineTo(cx,pad.t+ph);ctx.stroke();}
  for(let i=0;i<=4;i++){const cy=pad.t+ph*i/4;ctx.beginPath();ctx.moveTo(pad.l,cy);ctx.lineTo(pad.l+pw,cy);ctx.stroke();}

  ctx.strokeStyle='#d8d4cc'; ctx.lineWidth=1;
  ctx.beginPath();ctx.moveTo(pad.l,pad.t);ctx.lineTo(pad.l,pad.t+ph);ctx.stroke();
  ctx.beginPath();ctx.moveTo(pad.l,pad.t+ph);ctx.lineTo(pad.l+pw,pad.t+ph);ctx.stroke();

  ctx.fillStyle='#b0aca4'; ctx.font='9px "JetBrains Mono",monospace'; ctx.textAlign='center';
  for(let i=0;i<=6;i++) ctx.fillText((maxX*i/6).toFixed(1)+'m', pad.l+pw*i/6, pad.t+ph+20);
  ctx.textAlign='right';
  for(let i=0;i<=4;i++) ctx.fillText((maxY*(4-i)/4).toFixed(1), pad.l-7, pad.t+ph*i/4+3);

  // Ground line
  ctx.strokeStyle='#c8c4bc'; ctx.lineWidth=1.5;
  ctx.beginPath();ctx.moveTo(pad.l,toY(0));ctx.lineTo(pad.l+pw,toY(0));ctx.stroke();

  // Target opening
  const scale_ppx = pw / maxX;
  const txPx = toX(p.d), tyPx = toY(TARGET_H);
  const openPx = p.D_o * scale_ppx, ballPx = p.D_b * scale_ppx;

  ctx.save();ctx.setLineDash([5,5]);ctx.strokeStyle='#c42020';ctx.lineWidth=1;ctx.globalAlpha=.2;
  ctx.beginPath();ctx.moveTo(txPx,pad.t);ctx.lineTo(txPx,pad.t+ph);ctx.stroke();ctx.restore();
  ctx.save();ctx.strokeStyle='#c42020';ctx.lineWidth=5;ctx.lineCap='round';
  ctx.beginPath();ctx.moveTo(txPx-openPx/2,tyPx);ctx.lineTo(txPx+openPx/2,tyPx);ctx.stroke();
  ctx.lineWidth=2;ctx.globalAlpha=.5;
  ctx.beginPath();
  ctx.moveTo(txPx-openPx/2,tyPx);ctx.lineTo(txPx-openPx/2,tyPx+18);
  ctx.moveTo(txPx+openPx/2,tyPx);ctx.lineTo(txPx+openPx/2,tyPx+18);
  ctx.stroke();ctx.restore();
  ctx.fillStyle='#c42020';ctx.font='9px "JetBrains Mono",monospace';ctx.textAlign='center';
  ctx.fillText(`${p.d.toFixed(1)}m  ⌀${(p.D_o*1e3).toFixed(0)}mm`, txPx, pad.t+ph+20);

  // Pivot post
  ctx.strokeStyle='#c8c4bc';ctx.lineWidth=1.5;
  ctx.beginPath();ctx.moveTo(pad.l,toY(0));ctx.lineTo(pad.l,toY(p.h));ctx.stroke();
  ctx.fillStyle='#9a9690';ctx.beginPath();ctx.arc(pad.l,toY(p.h),5,0,Math.PI*2);ctx.fill();
  ctx.fillStyle='#fff';ctx.beginPath();ctx.arc(pad.l,toY(p.h),2,0,Math.PI*2);ctx.fill();

  if(!solNom) return;

  // ── Clip to plot area ──────────────────────────────────────────────────────
  ctx.save();
  ctx.beginPath();ctx.rect(pad.l-1, pad.t-1, pw+2, ph+2);ctx.clip();

  const pMin   = solMin   ? pathPoints(solMin.theta,   p) : null;
  const pMax   = solMax   ? pathPoints(solMax.theta,   p) : null;
  const pNom   = pathPoints(solNom.theta, p);
  const pShift = (solShift && p.ef_mult !== 0) ? pathPoints(solShift.theta, p) : null;

  // Band fill
  if(pMin && pMax){
    ctx.save();ctx.beginPath();
    pMin.forEach((pt,i)=>i===0?ctx.moveTo(toX(pt.x),toY(pt.y)):ctx.lineTo(toX(pt.x),toY(pt.y)));
    ctx.lineTo(toX(pMax[pMax.length-1].x), toY(pMax[pMax.length-1].y));
    for(let i=pMax.length-2;i>=0;i--) ctx.lineTo(toX(pMax[i].x),toY(pMax[i].y));
    ctx.closePath();ctx.fillStyle='rgba(26,95,196,0.10)';ctx.fill();ctx.restore();

    [pMin, pMax].forEach(pts=>{
      ctx.save();ctx.strokeStyle='#1a5fc4';ctx.lineWidth=1;ctx.globalAlpha=.18;ctx.setLineDash([4,4]);
      ctx.beginPath();pts.forEach((pt,i)=>i===0?ctx.moveTo(toX(pt.x),toY(pt.y)):ctx.lineTo(toX(pt.x),toY(pt.y)));
      ctx.stroke();ctx.restore();
    });
  }

  // Shifted arc
  if(pShift){
    ctx.save();ctx.strokeStyle='#c85a12';ctx.lineWidth=2;ctx.lineJoin='round';ctx.setLineDash([6,4]);ctx.globalAlpha=.8;
    ctx.beginPath();pShift.forEach((pt,i)=>i===0?ctx.moveTo(toX(pt.x),toY(pt.y)):ctx.lineTo(toX(pt.x),toY(pt.y)));
    ctx.stroke();ctx.restore();
    ctx.fillStyle='#c85a12';ctx.strokeStyle='#fff';ctx.lineWidth=1.5;
    ctx.beginPath();ctx.arc(toX(pShift[pShift.length-1].x),toY(pShift[pShift.length-1].y),3.5,0,Math.PI*2);ctx.fill();ctx.stroke();
  }

  // Nominal arc
  ctx.save();ctx.strokeStyle='#1a5fc4';ctx.lineWidth=2.5;ctx.lineJoin='round';
  ctx.beginPath();pNom.forEach((pt,i)=>i===0?ctx.moveTo(toX(pt.x),toY(pt.y)):ctx.lineTo(toX(pt.x),toY(pt.y)));
  ctx.stroke();ctx.restore();

  ctx.restore(); // end clip

  // Angle arc annotation
  const lx = pad.l, ly = toY(p.h);
  ctx.save();ctx.setLineDash([4,5]);ctx.strokeStyle='#1a5fc4';ctx.lineWidth=1.5;ctx.globalAlpha=.35;
  ctx.beginPath();ctx.moveTo(lx,ly);ctx.lineTo(toX(solNom.x0),toY(solNom.y0));ctx.stroke();ctx.restore();

  const th = solNom.theta * Math.PI / 180;
  ctx.save();ctx.strokeStyle='#1a5fc4';ctx.lineWidth=1.5;ctx.globalAlpha=.25;
  ctx.beginPath();ctx.arc(lx,ly,26,-th,0);ctx.stroke();ctx.restore();
  const midA = -th/2;
  ctx.fillStyle='#1a5fc4';ctx.font='bold 10px "JetBrains Mono",monospace';ctx.textAlign='left';
  ctx.fillText(solNom.theta.toFixed(1)+'°', lx+32*Math.cos(midA)+2, ly-32*Math.sin(midA)+4);

  // Launch dot
  ctx.fillStyle='#1a5fc4';ctx.strokeStyle='#fff';ctx.lineWidth=2;
  ctx.beginPath();ctx.arc(toX(solNom.x0),toY(solNom.y0),5,0,Math.PI*2);ctx.fill();ctx.stroke();

  // Peak dot
  const tp = solNom.uy / G;
  if(tp > 0){
    ctx.fillStyle='#1a5fc4';ctx.globalAlpha=.3;
    ctx.beginPath();ctx.arc(toX(solNom.x0+solNom.ux*tp),toY(solNom.ymax),4,0,Math.PI*2);ctx.fill();ctx.globalAlpha=1;
  }

  // Descent arrow at target
  const ax = toX(solNom.xAtTarget), ay = toY(TARGET_H);
  const spd = Math.sqrt(solNom.ux*solNom.ux + solNom.vy_at*solNom.vy_at);
  const arrowLen = 30;
  const dxN = solNom.ux / spd, dyN = -solNom.vy_at / spd;
  const bx = ax - dxN*arrowLen, by = ay - dyN*arrowLen;
  ctx.save();ctx.strokeStyle='#c85a12';ctx.lineWidth=2;ctx.globalAlpha=.9;
  ctx.beginPath();ctx.moveTo(bx,by);ctx.lineTo(ax,ay);ctx.stroke();
  const ang = Math.atan2(ay-by, ax-bx);
  ctx.fillStyle='#c85a12';
  ctx.beginPath();ctx.moveTo(ax,ay);ctx.lineTo(ax-8*Math.cos(ang-0.38),ay-8*Math.sin(ang-0.38));ctx.lineTo(ax-8*Math.cos(ang+0.38),ay-8*Math.sin(ang+0.38));ctx.closePath();ctx.fill();
  ctx.font='bold 9px "JetBrains Mono",monospace';ctx.textAlign='right';
  ctx.fillText(solNom.descentDeg.toFixed(1)+'° ↓', ax-3, by-4);
  ctx.restore();

  // Ball outline at target
  if(ballPx > 3){
    ctx.save();ctx.globalAlpha=.18;ctx.strokeStyle='#1a5fc4';ctx.lineWidth=1.5;
    ctx.beginPath();ctx.arc(ax, ay-ballPx/2, ballPx/2, 0, Math.PI*2);ctx.stroke();ctx.restore();
  }
  ctx.fillStyle='#1a5fc4';ctx.strokeStyle='#fff';ctx.lineWidth=2;
  ctx.beginPath();ctx.arc(ax,ay,4,0,Math.PI*2);ctx.fill();ctx.stroke();
}

function calculate(){
  const p = gv();
  saveStorage();
  const delta  = p.d * p.e_f;
  const dMin   = Math.max(0.1, p.d - delta);
  const dMax   = p.d + delta;
  const dShift = p.d * (1 + p.ef_mult);
  const solNom   = findAngle(p.d,   p);
  const solMin   = findAngle(dMin,  p);
  const solMax   = findAngle(dMax,  p);
  const solShift = p.ef_mult !== 0 ? findAngle(Math.max(0.1, dShift), p) : null;
  renderResults(solNom, solMin, solMax, solShift, p);
  drawTraj(solNom, solMin, solMax, solShift, p);
}

// Wire up all controls
for(const[key,d] of Object.entries(defs)){
  const sl = document.getElementById(d.s);
  const nm = document.getElementById(d.n);
  if(!sl||!nm) continue;
  sl.addEventListener('input', ()=>{
    const param = d.sv(+sl.value);
    nm.value = d.nmGet(param);
    if(d.dv){const dv=document.getElementById(d.dv);if(dv)dv.textContent=d.fmt(param);}
    if(key==='ef_mult') updateEfTag(param);
    calculate();
  });
  nm.addEventListener('input', ()=>{
    const raw = nm.value;
    if(raw===''||isNaN(+raw)) return;
    const param = d.nmSet(+raw);
    sl.value = d.vs(param);
    if(d.dv){const dv=document.getElementById(d.dv);if(dv)dv.textContent=d.fmt(param);}
    if(key==='ef_mult') updateEfTag(param);
    calculate();
  });
}

// Step-size buttons
const scaleBtns = document.querySelectorAll('.scale-btn');
let currentStep = 0.1;

function updateSliderSteps(){
  for(const[key,d] of Object.entries(defs)){
    const sl = document.getElementById(d.s);
    if(!sl) continue;
    const useIntStep = (key==='M'||key==='m');
    if(useIntStep) continue;
    // Convert display-space step → slider-space delta
    const paramDelta = d.nmSet ? d.nmSet(currentStep) - d.nmSet(0) : currentStep;
    const sliderDelta = Math.max(1, Math.abs(Math.round(d.vs(paramDelta) - d.vs(0))));
    sl.step = sliderDelta;
  }
}

scaleBtns.forEach(btn=>{
  btn.addEventListener('click', ()=>{
    scaleBtns.forEach(b=>b.classList.remove('active'));
    btn.classList.add('active');
    currentStep = parseFloat(btn.dataset.step);
    updateSliderSteps();
  });
});
document.querySelector('.scale-btn[data-step="0.1"]').classList.add('active');
updateSliderSteps();

// ± buttons
document.querySelectorAll('.pm-btn').forEach(btn=>{
  btn.addEventListener('click', ()=>{
    const field = btn.dataset.field;
    const dir   = +btn.dataset.dir;
    const d     = defs[field];
    const sl    = document.getElementById(d.s);
    const nm    = document.getElementById(d.n);
    if(!sl||!nm) return;
    const useIntStep = (field==='M'||field==='m'||field==='D_o'||field==='D_b');
    const step = useIntStep ? 1 : currentStep;
    // Read directly from number input (not slider) to preserve sub-step precision
    let nmNum = parseFloat(nm.value);
    // Safe float addition — avoid 0.1+0.2 drift
    nmNum = Math.round((nmNum + dir * step) * 1e9) / 1e9;
    const nmMin = parseFloat(nm.min), nmMax = parseFloat(nm.max);
    nmNum = Math.max(nmMin, Math.min(nmMax, nmNum));
    const newParam = d.nmSet(nmNum);
    nm.value = d.nmGet(newParam);
    sl.value = d.vs(newParam);
    if(d.dv){const dv=document.getElementById(d.dv);if(dv)dv.textContent=d.fmt(newParam);}
    if(field==='ef_mult') updateEfTag(newParam);
    calculate();
  });
});

loadStorage();
calculate();
window.addEventListener('resize', calculate);
});