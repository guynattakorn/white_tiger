/**
 * WHITE TIGER · Projectile Calculator — main.js v6.0
 *
 * Architecture:
 *   ConfigManager     — Constants (physics, render, fields)
 *   PhysicsEngine     — Projectile math (no DOM)
 *   TimingEngine      — Timing math (no DOM)
 *   AngleController   — Angle state & sidebar UI
 *   CanvasRenderer    — Canvas 2D trajectory drawing
 *   UIRenderer        — DOM state switches & status updates
 *   InputController   — Input events & value management
 *   TimingController  — Timing UI in right panel
 *   HistoryController — Saved presets with localStorage
 *   SimulationApp     — Mediator that wires all modules
 */

'use strict';

/* ── 1. ConfigManager ── */
class ConfigManager {
  static PHYSICS = Object.freeze({
    G: 9.81,          // m/s²
    Y_DIFF: 0.05,     // m — height difference (target - launch)
    Y_LAUNCH: 0.38,   // m — launch height (380 mm)
    Y_TARGET: 0.43,   // m — target height (430 mm)
    THETA_MIN: 46,    // ° — valid angle lower bound
    THETA_MAX: 73,    // ° — valid angle upper bound
  });

  static RENDER = Object.freeze({
    TRAJ_COLOR: '#3b82f6',
    TRAJ_BAND_COLOR: 'rgba(59,130,246,0.07)',
    TARGET_COLOR: '#ef4444',
    GRID_COLOR: '#f0f0f0',
    AXIS_COLOR: '#e0e0e0',
    TEXT_COLOR: '#888888',
    TEXT_DARK: '#555555',
    PAD: { t: 48, b: 42, l: 58, r: 28 },
  });

  static FIELDS = Object.freeze({
    sx: { id: 'in-sx', min: 0, max: 2680, step: 0.1, dec: 1 },
  });

  static TIMING_FIELDS = Object.freeze({
    trig_angle: { id: 'in-trig-angle', min: 0, max: 360, step: 1, dec: 1 },
    rpm:        { id: 'in-rpm',        min: 1, max: 500, step: 1, dec: 1 },
    n:          { id: 'in-n',          min: 0, max: 100, step: 1, dec: 0 },
    t_delay:    { id: 'in-t-delay',    min: 0, max: 100, step: 0.001, dec: 3 },
  });

  static DEFAULTS = Object.freeze({
    sx:         163.5,
    trig_angle: 180,
    rpm:        10,
    n:          0,
    t_delay:    0,
  });

  static Y_MEASURE_DEFAULT = 80;
  static ANGLE_STEP = 0.1;
}


/* ── 1.5 DirtyTracker ── */
class DirtyTracker {
  constructor() {
    this._baseline = new Map();
    this._dirtyIds = new Set();
    this._isDirty = false;
    this._hasCalculated = false;

    this._btnCalc = document.getElementById('btn-calc');
    this._statusDot = document.getElementById('status-dot');
    this._statusTxt = document.getElementById('status-txt');
  }

  static EXEMPT_IDS = new Set(['in-y-measure', 'in-sx']);

  snapshot(inputIds) {
    for (const id of inputIds) {
      const el = document.getElementById(id);
      if (el) this._baseline.set(id, el.value);
    }
    this._dirtyIds.clear();
    this._isDirty = false;
    this._hasCalculated = true;
    this._updateGlobalUI();
    this._updateFieldUI();
  }

  check(inputId) {
    if (DirtyTracker.EXEMPT_IDS.has(inputId)) return;
    const el = document.getElementById(inputId);
    if (!el) return;
    const baseVal = this._baseline.get(inputId);
    if (baseVal === undefined) {
      if (!this._hasCalculated) return;
      this._dirtyIds.add(inputId);
    } else if (el.value !== baseVal) {
      this._dirtyIds.add(inputId);
    } else {
      this._dirtyIds.delete(inputId);
    }
    this._isDirty = this._dirtyIds.size > 0;
    this._updateGlobalUI();
    this._updateFieldUI();
  }

  get isDirty() { return this._isDirty; }

  clear() {
    this._dirtyIds.clear();
    this._isDirty = false;
    this._hasCalculated = false;
    this._baseline.clear();
    this._updateGlobalUI();
    this._updateFieldUI();
  }

  markResultsStale(stale) {
    const cards = document.querySelectorAll('#result-primary, #result-timing, #result-electrical');
    cards.forEach(c => c.classList.toggle('result--stale', stale));
  }

  _updateGlobalUI() {
    if (this._btnCalc) {
      this._btnCalc.classList.toggle('btn-calc--dirty', this._isDirty);
    }
    if (this._statusDot && this._isDirty) {
      this._statusDot.className = 'status-dot amber-pulse';
    }
    if (this._statusTxt && this._isDirty) {
      this._statusTxt.textContent = 'ค่ามีการเปลี่ยนแปลง — กด CALCULATE เพื่ออัปเดต';
    }
    this.markResultsStale(this._isDirty && this._hasCalculated);
  }

  _updateFieldUI() {
    for (const [id] of this._baseline) {
      const input = document.getElementById(id);
      if (!input) continue;
      const prow = input.closest('.prow');
      const isDirty = this._dirtyIds.has(id);
      if (isDirty) {
        input.classList.add('input-dirty');
        if (prow) prow.classList.add('prow--dirty');
      } else {
        input.classList.remove('input-dirty');
        if (prow) prow.classList.remove('prow--dirty');
      }
    }
  }
}


/* ── 2. PhysicsEngine ── */
class PhysicsEngine {
  constructor(cfg = ConfigManager.PHYSICS) {
    this._c = cfg;
  }

  /**
   * Solve for launch angle using empirical quadratic:
   *   sx = 86 + 8.22·θ − 0.103·θ²
   */
  solve(sx_cm) {
    const a = -0.103;
    const b = 8.22;
    const c = 86 - sx_cm;
    const discriminant = b * b - 4 * a * c;
    if (discriminant < 0) return null;
    const theta1 = (-b + Math.sqrt(discriminant)) / (2 * a);
    const theta2 = (-b - Math.sqrt(discriminant)) / (2 * a);
    const theta = Math.max(theta1, theta2);
    const v0 = this.calculateV0ForTarget(theta, sx_cm / 100);
    return Object.freeze({ theta, v0, x: sx_cm / 100 });
  }

  calculateV0ForTarget(thetaDeg, x_m) {
    const r = thetaDeg * Math.PI / 180;
    const cosSq = Math.pow(Math.cos(r), 2);
    const num = this._c.G * x_m * x_m;
    const den = 2 * cosSq * (x_m * Math.tan(r) - this._c.Y_DIFF);
    if (den <= 0) return 5;
    return Math.sqrt(num / den);
  }

  trajectoryPoints(thetaDeg, v0) {
    const pts = [];
    const { ux, uy } = this._vel(thetaDeg, v0);
    const dt = 0.005;
    let t = 0;
    while (true) {
      const x = ux * t;
      const y = this._c.Y_LAUNCH + uy * t - 0.5 * this._c.G * t * t;
      pts.push({ x, y });
      if ((y < 0 && t > 0) || t > 10) break;
      t += dt;
    }
    return pts;
  }

  landingX(thetaDeg, v0) {
    const { ux, uy } = this._vel(thetaDeg, v0);
    const a = 0.5 * this._c.G;
    const b = -uy;
    const c = this._c.Y_DIFF;
    const d = b * b - 4 * a * c;
    if (d < 0) return null;
    const t = (-b + Math.sqrt(d)) / (2 * a);
    return ux * t;
  }

  /** Compute flight metrics for the metrics bar. */
  computeMetrics(thetaDeg, v0) {
    const { ux, uy } = this._vel(thetaDeg, v0);

    // Time of flight to Y_TARGET
    const a = 0.5 * this._c.G;
    const b = -uy;
    const c = this._c.Y_DIFF;
    const d = b * b - 4 * a * c;
    const tof = d >= 0 ? (-b + Math.sqrt(d)) / (2 * a) : 0;

    // Peak height
    const tPeak = uy / this._c.G;
    const ymax = this._c.Y_LAUNCH + uy * tPeak - 0.5 * this._c.G * tPeak * tPeak;

    // Descent angle at target
    const vyAtTarget = uy - this._c.G * tof;
    const descentDeg = Math.abs(Math.atan2(-vyAtTarget, ux) * 180 / Math.PI);

    return { v0, tof, ymax, descentDeg };
  }

  _vel(thetaDeg, v0) {
    const r = thetaDeg * Math.PI / 180;
    return { ux: v0 * Math.cos(r), uy: v0 * Math.sin(r) };
  }
}


/* ── 3. AngleController ── */
class AngleController {
  constructor(step = ConfigManager.ANGLE_STEP) {
    this._step = step;
    this._target = null;
    this._current = null;
    this._onAdjust = null;
    this._els = {
      group: document.getElementById('angle-result-group'),
      target: document.getElementById('target-angle-val'),
      display: document.getElementById('ca-display'),
      plus: document.getElementById('ca-plus'),
      minus: document.getElementById('ca-minus'),
      panel: document.getElementById('delta-panel'),
      icon: document.getElementById('dp-icon'),
      val: document.getElementById('dp-val'),
      sub: document.getElementById('dp-sublabel'),
    };
  }

  setTarget(theta) {
    this._target = theta;
    this._current = this._snap(theta);
    this._renderTarget();
    this._renderCurrent();
    this._renderDelta();
  }

  onAdjust(cb) { this._onAdjust = cb; }
  get current() { return this._current; }
  get target() { return this._target; }

  bindButtons() {
    this._els.plus?.addEventListener('click', () => this._adjust(+1));
    this._els.minus?.addEventListener('click', () => this._adjust(-1));
  }

  show() { if (this._els.group) this._els.group.style.display = ''; }
  hide() { if (this._els.group) this._els.group.style.display = 'none'; }

  _adjust(dir) {
    if (this._current === null) return;
    this._current = Math.max(10, Math.min(89,
      this._snap(this._current + dir * this._step)
    ));
    this._renderCurrent();
    this._renderDelta();
    this._onAdjust?.(this._current);
  }

  _snap(v) { return Math.round(v / this._step) * this._step; }

  _renderTarget() {
    if (this._els.target) this._els.target.textContent = this._target?.toFixed(2) ?? '--';
  }

  _renderCurrent() {
    if (this._els.display) this._els.display.textContent = this._current?.toFixed(1) ?? '--';
  }

  _renderDelta() {
    if (!this._els.panel || this._current === null || this._target === null) return;
    const delta = this._target - this._current;
    const absDelta = Math.abs(delta);
    const { panel, icon, val, sub } = this._els;
    panel.classList.remove('on-target', 'need-up', 'need-down');
    if (absDelta <= 0.05) {
      panel.classList.add('on-target');
      icon.textContent = '✓'; val.textContent = 'ON TARGET'; sub.textContent = 'ถึงเป้าหมายแล้ว!';
    } else if (delta > 0) {
      panel.classList.add('need-up');
      icon.textContent = '↑'; val.textContent = `${absDelta.toFixed(1)}°`; sub.textContent = 'เพิ่มมุมขึ้น';
    } else {
      panel.classList.add('need-down');
      icon.textContent = '↓'; val.textContent = `${absDelta.toFixed(1)}°`; sub.textContent = 'ลดมุมลง';
    }
  }
}


/* ── 4. CanvasRenderer — Professional graph matching reference ── */
class CanvasRenderer {
  constructor(canvasId = 'traj', rcfg = ConfigManager.RENDER, pcfg = ConfigManager.PHYSICS) {
    this._cv = document.getElementById(canvasId);
    this._rcfg = rcfg;
    this._pcfg = pcfg;
  }

  draw({ targetPts, targetAngle, currentPts, currentAngle, currentLandX, xTarget, bandPtsUpper, bandPtsLower }) {
    this._resize();
    const ctx = this._cv.getContext('2d');
    const W = this._cv.width;
    const H = this._cv.height;
    ctx.clearRect(0, 0, W, H);

    // White background
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, W, H);

    const allPts = [...targetPts, ...currentPts];
    if (!allPts.length) return;

    const { bounds, toX, toY } = this._coords(allPts, xTarget);
    const { pad } = bounds;

    this._drawGrid(ctx, bounds);
    this._drawAxes(ctx, bounds);
    this._drawAxisLabels(ctx, bounds);

    // Reference lines
    this._drawRefLine(ctx, bounds, toY, this._pcfg.Y_TARGET, '#ef4444', 'target y = 430mm');
    this._drawRefLine(ctx, bounds, toY, this._pcfg.Y_LAUNCH, '#d0d0d0', 'launch y = 380mm');

    const isOnTarget = Math.abs(currentAngle - targetAngle) <= 0.05;

    // Trajectory band (shaded area ±2° around target)
    if (bandPtsUpper && bandPtsLower) {
      this._drawBand(ctx, bandPtsUpper, bandPtsLower, toX, toY, this._rcfg.TRAJ_BAND_COLOR);
    }

    // Target trajectory (dashed blue)
    if (!isOnTarget) {
      this._drawPath(ctx, targetPts, toX, toY, this._rcfg.TRAJ_COLOR, 1.5, [6, 4], 0.35);
    }

    // Current trajectory (solid blue)
    this._drawPath(ctx, currentPts, toX, toY, this._rcfg.TRAJ_COLOR, 2.5, null, 1);

    // Target marker (red vertical line + dot)
    this._drawTargetMarker(ctx, toX, toY, xTarget, this._pcfg.Y_TARGET);

    // Landing dot
    if (currentLandX !== null) {
      this._drawLandingDot(ctx, toX, toY, currentLandX, this._pcfg.Y_TARGET);
    }

    // Launch dot
    this._drawLaunchDot(ctx, toX, toY);

    // Angle arc
    this._drawAngleArc(ctx, toX, toY, targetAngle);

    // Legend
    this._drawLegend(isOnTarget, targetAngle, currentAngle);
  }

  clear() {
    this._resize();
    const ctx = this._cv.getContext('2d');
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, this._cv.width, this._cv.height);
  }

  _resize() {
    const p = this._cv.parentElement;
    this._cv.width = p.clientWidth;
    this._cv.height = p.clientHeight;
  }

  _coords(pts, xTarget) {
    const pad = { ...this._rcfg.PAD };
    const W = this._cv.width - pad.l - pad.r;
    const H = this._cv.height - pad.t - pad.b;
    const xs = pts.map(p => p.x);
    const ys = pts.map(p => p.y);
    const xMax = Math.max(...xs, xTarget) * 1.15;
    const yMin = 0;
    const yMax = Math.max(...ys, this._pcfg.Y_TARGET) * 1.2;
    const toX = v => pad.l + (v / xMax) * W;
    const toY = v => pad.t + H - ((v - yMin) / (yMax - yMin)) * H;
    return { bounds: { pad, W, H, xMax, yMin, yMax }, toX, toY };
  }

  _drawGrid(ctx, { pad, W, H, xMax, yMin, yMax }) {
    ctx.save();
    ctx.strokeStyle = this._rcfg.GRID_COLOR;
    ctx.lineWidth = 1;

    // Compute nice grid intervals
    const xStep = this._niceStep(xMax, 6);
    const yStep = this._niceStep(yMax - yMin, 5);

    // Vertical grid
    for (let xv = xStep; xv < xMax; xv += xStep) {
      const gx = pad.l + (xv / xMax) * W;
      ctx.beginPath(); ctx.moveTo(gx, pad.t); ctx.lineTo(gx, pad.t + H); ctx.stroke();
    }

    // Horizontal grid
    for (let yv = yMin + yStep; yv < yMax; yv += yStep) {
      const gy = pad.t + H - ((yv - yMin) / (yMax - yMin)) * H;
      ctx.beginPath(); ctx.moveTo(pad.l, gy); ctx.lineTo(pad.l + W, gy); ctx.stroke();
    }
    ctx.restore();
  }

  _drawAxes(ctx, { pad, W, H }) {
    ctx.save();
    ctx.strokeStyle = this._rcfg.AXIS_COLOR;
    ctx.lineWidth = 1.5;
    // Y axis
    ctx.beginPath(); ctx.moveTo(pad.l, pad.t); ctx.lineTo(pad.l, pad.t + H); ctx.stroke();
    // X axis
    ctx.beginPath(); ctx.moveTo(pad.l, pad.t + H); ctx.lineTo(pad.l + W, pad.t + H); ctx.stroke();
    ctx.restore();
  }

  _drawAxisLabels(ctx, { pad, W, H, xMax, yMin, yMax }) {
    ctx.save();
    ctx.fillStyle = this._rcfg.TEXT_COLOR;
    ctx.font = '10px "IBM Plex Mono", monospace';

    const xStep = this._niceStep(xMax, 6);
    const yStep = this._niceStep(yMax - yMin, 5);

    // X labels
    ctx.textAlign = 'center';
    for (let xv = 0; xv <= xMax; xv += xStep) {
      const gx = pad.l + (xv / xMax) * W;
      const label = xv >= 1 ? xv.toFixed(1) + 'm' : (xv * 100).toFixed(0) + 'cm';
      ctx.fillText(label, gx, pad.t + H + 18);
    }

    // Y labels
    ctx.textAlign = 'right';
    for (let yv = yMin; yv <= yMax + 0.001; yv += yStep) {
      const gy = pad.t + H - ((yv - yMin) / (yMax - yMin)) * H;
      ctx.fillText(yv.toFixed(1), pad.l - 8, gy + 4);
    }
    ctx.restore();
  }

  /** Compute a "nice" step size for grid lines. */
  _niceStep(range, targetSteps) {
    const rough = range / targetSteps;
    const mag = Math.pow(10, Math.floor(Math.log10(rough)));
    const residual = rough / mag;
    let nice;
    if (residual <= 1.5) nice = 1;
    else if (residual <= 3.5) nice = 2;
    else if (residual <= 7.5) nice = 5;
    else nice = 10;
    return nice * mag;
  }

  _drawRefLine(ctx, { pad, W }, toY, yVal, color, label) {
    const gy = toY(yVal);
    ctx.save();
    ctx.strokeStyle = color;
    ctx.lineWidth = 1;
    ctx.setLineDash([5, 3]);
    ctx.globalAlpha = 0.6;
    ctx.beginPath(); ctx.moveTo(pad.l, gy); ctx.lineTo(pad.l + W, gy); ctx.stroke();
    ctx.setLineDash([]);

    // Label
    ctx.globalAlpha = 0.7;
    ctx.fillStyle = color;
    ctx.font = '9px "IBM Plex Mono", monospace';
    ctx.textAlign = 'left';
    ctx.fillText(label, pad.l + 6, gy - 5);
    ctx.restore();
  }

  /** Draw shaded band between upper and lower trajectory bounds. */
  _drawBand(ctx, upper, lower, toX, toY, color) {
    if (!upper.length || !lower.length) return;
    ctx.save();
    ctx.fillStyle = color;
    ctx.beginPath();
    ctx.moveTo(toX(upper[0].x), toY(upper[0].y));
    for (let i = 1; i < upper.length; i++) {
      ctx.lineTo(toX(upper[i].x), toY(upper[i].y));
    }
    for (let i = lower.length - 1; i >= 0; i--) {
      ctx.lineTo(toX(lower[i].x), toY(lower[i].y));
    }
    ctx.closePath();
    ctx.fill();
    ctx.restore();
  }

  _drawPath(ctx, pts, toX, toY, color, width, dash, alpha) {
    if (!pts.length) return;
    ctx.save();
    ctx.globalAlpha = alpha;
    ctx.strokeStyle = color;
    ctx.lineWidth = width;
    ctx.lineJoin = 'round';
    ctx.lineCap = 'round';
    if (dash) ctx.setLineDash(dash);
    ctx.beginPath();
    ctx.moveTo(toX(pts[0].x), toY(pts[0].y));
    for (let i = 1; i < pts.length; i++) {
      ctx.lineTo(toX(pts[i].x), toY(pts[i].y));
    }
    ctx.stroke();
    ctx.restore();
  }

  /** Red target marker — vertical dashed line + crosshair dot. */
  _drawTargetMarker(ctx, toX, toY, x, y) {
    const px = toX(x), py = toY(y);
    const { pad, H } = this._lastBounds || {};

    ctx.save();
    // Vertical dashed line
    ctx.strokeStyle = this._rcfg.TARGET_COLOR;
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 3]);
    ctx.globalAlpha = 0.4;
    ctx.beginPath();
    ctx.moveTo(px, toY(0));
    ctx.lineTo(px, py - 12);
    ctx.stroke();
    ctx.setLineDash([]);

    // Crosshair ring
    ctx.globalAlpha = 0.3;
    ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.arc(px, py, 10, 0, Math.PI * 2); ctx.stroke();

    // Center dot
    ctx.globalAlpha = 1;
    ctx.fillStyle = this._rcfg.TARGET_COLOR;
    ctx.beginPath(); ctx.arc(px, py, 4, 0, Math.PI * 2); ctx.fill();

    // Label
    ctx.globalAlpha = 0.85;
    ctx.font = 'bold 9px "IBM Plex Mono", monospace';
    ctx.textAlign = 'left';
    ctx.fillText(`Target`, px + 14, py - 6);
    ctx.globalAlpha = 0.6;
    ctx.font = '9px "IBM Plex Mono", monospace';
    ctx.fillText(`(${(x*100).toFixed(0)}cm, ${(y*1000).toFixed(0)}mm)`, px + 14, py + 6);
    ctx.restore();
  }

  _drawLandingDot(ctx, toX, toY, landX, yTarget) {
    const px = toX(landX), py = toY(yTarget);
    ctx.save();
    ctx.fillStyle = this._rcfg.TRAJ_COLOR;
    ctx.strokeStyle = '#ffffff';
    ctx.lineWidth = 2.5;
    ctx.beginPath(); ctx.arc(px, py, 5, 0, Math.PI * 2);
    ctx.fill(); ctx.stroke();

    // Distance label
    ctx.fillStyle = this._rcfg.TEXT_DARK;
    ctx.font = '9px "IBM Plex Mono", monospace';
    ctx.textAlign = 'center';
    ctx.fillText(`${(landX*100).toFixed(0)}cm`, px, py - 10);
    ctx.restore();
  }

  _drawLaunchDot(ctx, toX, toY) {
    const px = toX(0), py = toY(this._pcfg.Y_LAUNCH);
    ctx.save();
    ctx.fillStyle = '#555555';
    ctx.strokeStyle = '#ffffff';
    ctx.lineWidth = 2;
    ctx.beginPath(); ctx.arc(px, py, 4, 0, Math.PI * 2);
    ctx.fill(); ctx.stroke();
    ctx.restore();
  }

  _drawAngleArc(ctx, toX, toY, thetaDeg) {
    const thRad = thetaDeg * Math.PI / 180;
    const ox = toX(0), oy = toY(this._pcfg.Y_LAUNCH);
    const R = 28;
    ctx.save();
    ctx.strokeStyle = this._rcfg.TRAJ_COLOR;
    ctx.lineWidth = 1.5;
    ctx.globalAlpha = 0.3;
    ctx.beginPath(); ctx.arc(ox, oy, R, -thRad, 0); ctx.stroke();

    // Angle label
    ctx.globalAlpha = 1;
    ctx.fillStyle = this._rcfg.TEXT_DARK;
    ctx.font = 'bold 11px "IBM Plex Mono", monospace';
    ctx.textAlign = 'left';
    const midA = -thRad / 2;
    const LR = R + 12;
    ctx.fillText(thetaDeg.toFixed(1) + '°', ox + LR * Math.cos(midA) + 2, oy + LR * Math.sin(midA) + 4);
    ctx.restore();
  }

  _drawLegend(isOnTarget, targetAngle, currentAngle) {
    const el = document.getElementById('legend-row');
    if (!el) return;
    const blue = this._rcfg.TRAJ_COLOR;
    if (isOnTarget) {
      el.innerHTML = `
        <div class="leg">
          <div class="leg-swatch" style="background:${blue}"></div>
          <span>θ = ${currentAngle.toFixed(2)}° · On Target ✓</span>
        </div>
        <div class="leg">
          <div class="leg-swatch" style="background:rgba(59,130,246,.15);border:1px solid rgba(59,130,246,.3)"></div>
          <span>Band</span>
        </div>`;
    } else {
      el.innerHTML = `
        <div class="leg">
          <div class="leg-swatch" style="background:${blue}"></div>
          <span>Nominal ${targetAngle.toFixed(2)}°</span>
        </div>
        <div class="leg">
          <div class="leg-swatch" style="background:rgba(59,130,246,.15);border:1px solid rgba(59,130,246,.3)"></div>
          <span>Band</span>
        </div>
        <div class="leg">
          <div class="leg-swatch" style="background:${blue};opacity:.5"></div>
          <span>Current ${currentAngle.toFixed(1)}°</span>
        </div>
        <div class="leg">
          <div class="leg-swatch" style="background:${this._rcfg.TARGET_COLOR}"></div>
          <span>Target</span>
        </div>`;
    }
  }
}


/* ── 5. UIRenderer ── */
class UIRenderer {
  constructor() {
    this._els = {
      empty:      document.getElementById('empty-state'),
      noSol:      document.getElementById('no-sol-wrap'),
      graphWrap:  document.getElementById('graph-wrap'),
      rangeWarn:  document.getElementById('range-warn'),
      statusTxt:  document.getElementById('status-txt'),
      statusDot:  document.getElementById('status-dot'),
      barWidget:  document.getElementById('bar-angle-widget'),
      toast:      document.getElementById('reset-toast'),
      metricsBar: document.getElementById('metrics-bar'),
      rpIdle:     document.getElementById('rp-idle'),
      rpSolved:   document.getElementById('rp-solved'),
      rpAngleVal: document.getElementById('rp-angle-val'),
    };
  }

  showIdle() {
    this._show(this._els.empty);
    this._hide(this._els.noSol);
    this._hide(this._els.graphWrap);
    this._hide(this._els.barWidget);
    this._hide(this._els.metricsBar);
    this._els.rangeWarn.classList.remove('show');
    this._setStatus('Ready — กรอกค่าแล้วกด CALCULATE', 'green');
    this._resetRightPanel();
    this._resetMetrics();
  }

  showNoSolution() {
    this._hide(this._els.empty);
    this._show(this._els.noSol);
    this._els.graphWrap.style.display = 'flex';
    this._hide(this._els.barWidget);
    this._hide(this._els.metricsBar);
    this._els.rangeWarn.classList.add('show');
    this._setStatus('ไม่พบคำตอบในช่วง 46°–73°', 'red');
    this._resetRightPanel();
    this._resetMetrics();
  }

  showSolved(targetAngle, metrics) {
    this._hide(this._els.empty);
    this._hide(this._els.noSol);
    this._els.graphWrap.style.display = 'flex';
    this._show(this._els.barWidget);
    this._show(this._els.metricsBar);
    this._els.rangeWarn.classList.remove('show');
    this._setStatus(`θ = ${targetAngle.toFixed(2)}° · descent ${metrics.descentDeg.toFixed(1)}°`, 'green');
    this._updateRightPanel(targetAngle);
    this._updateMetrics(metrics);
  }

  updateLiveStatus(currentAngle) {
    if (this._els.barWidget) {
      const valEl = document.getElementById('baw-val');
      if (valEl) valEl.textContent = currentAngle.toFixed(1);
    }
  }

  showToast() {
    const t = this._els.toast;
    t.classList.add('show');
    setTimeout(() => t.classList.remove('show'), 2200);
  }

  _updateRightPanel(theta) {
    this._hide(this._els.rpIdle);
    this._show(this._els.rpSolved);
    if (this._els.rpAngleVal) this._els.rpAngleVal.textContent = theta.toFixed(1);
  }

  _resetRightPanel() {
    this._show(this._els.rpIdle);
    this._hide(this._els.rpSolved);
  }

  _updateMetrics(m) {
    this._setEl('met-v0', m.v0.toFixed(2));
    this._setEl('met-tof', m.tof.toFixed(3));
    this._setEl('met-peak', m.ymax.toFixed(2));
    this._setEl('met-descent', m.descentDeg.toFixed(1));
  }

  _resetMetrics() {
    ['met-v0', 'met-tof', 'met-peak', 'met-descent'].forEach(id => this._setEl(id, '--'));
  }

  _show(el) { if (el) el.style.display = ''; }
  _hide(el) { if (el) el.style.display = 'none'; }
  _setEl(id, text) { const el = document.getElementById(id); if (el) el.textContent = text; }

  _setStatus(text, dotState) {
    const el = this._els.statusTxt;
    if (el) el.textContent = text;
    if (this._els.statusDot) {
      this._els.statusDot.className = `status-dot ${dotState}`;
    }
  }
}


/* ── 6. InputController ── */
class InputController {
  constructor() {
    this._projectileFields = ConfigManager.FIELDS;
    this._timingFields = ConfigManager.TIMING_FIELDS;
    this._defaults = ConfigManager.DEFAULTS;
    this._allFields = { ...this._projectileFields, ...this._timingFields };
    this._dirtyTracker = null;
  }

  setDirtyTracker(tracker) { this._dirtyTracker = tracker; }

  bindEvents() {
    for (const key of Object.keys(this._allFields)) {
      const f = this._allFields[key];
      const el = document.getElementById(f.id);
      if (!el) continue;
      el.addEventListener('change', () => {
        let v = parseFloat(el.value);
        if (isNaN(v)) v = f.min;
        v = Math.max(f.min, Math.min(f.max, v));
        el.value = v.toFixed(f.dec);
        if (this._dirtyTracker) this._dirtyTracker.check(f.id);
      });
      el.addEventListener('input', () => {
        if (this._dirtyTracker) this._dirtyTracker.check(f.id);
      });
    }

    // Target Calculator sync (y_measure ↔ sx)
    const inYMeasure = document.getElementById('in-y-measure');
    const inSx = document.getElementById('in-sx');
    if (inYMeasure && inSx) {
      let isUpdating = false;
      const updateFromY = () => {
        if (isUpdating) return;
        isUpdating = true;
        const y = parseFloat(inYMeasure.value) || 0;
        inSx.value = (250 - 6.5 - y).toFixed(1);
        isUpdating = false;
      };
      const updateFromSx = () => {
        if (isUpdating) return;
        isUpdating = true;
        const sx = parseFloat(inSx.value) || 0;
        inYMeasure.value = (250 - 6.5 - sx).toFixed(1);
        isUpdating = false;
      };
      inYMeasure.addEventListener('input', updateFromY);
      inYMeasure.addEventListener('change', updateFromY);
      inSx.addEventListener('input', updateFromSx);
      inSx.addEventListener('change', updateFromSx);
      updateFromY();
    }
  }

  onSxAutoChanged(callback) { this._onSxAutoChangedCallback = callback; }

  getSx() {
    const el = document.getElementById('in-sx');
    return el ? parseFloat(el.value) || 0 : 0;
  }

  resetToDefaults() {
    for (const [key, val] of Object.entries(this._defaults)) {
      const field = this._allFields[key];
      if (!field) continue;
      const el = document.getElementById(field.id);
      if (el) el.value = typeof val === 'number' ? val.toFixed(field.dec) : val;
    }
    const yEl = document.getElementById('in-y-measure');
    if (yEl) yEl.value = ConfigManager.Y_MEASURE_DEFAULT;
  }
}


/* ── 7. TimingEngine ── */
class TimingEngine {
  constructor(G = 9.81) { this._G = G; }

  /**
   * R = t_delay / (1.10251450676983 × 470000 × 0.000001)
   * V = (12.45 × 3.2) / R
   */
  compute(t_delay) {
    const k = 1.10251450676983 * 470 * 1000 * 0.000001;
    const r = t_delay / k;
    const v = r !== 0 ? (12.45 * 3.2) / r : 0;
    return { ok: true, t_delay, r_total: r, volt: v };
  }
}


/* ── 8. TimingController ── */
class TimingController {
  constructor(engine, inputCtrl) {
    this._engine = engine;
    this._input = inputCtrl;
    this._lastTheta = null;
    this._lastV0 = null;
    this._lastSy = null;
    this._els = {
      idle:       document.getElementById('rt-idle'),
      solved:     document.getElementById('rt-solved'),
      error:      document.getElementById('rt-error'),
      badge:      document.getElementById('rt-badge'),
      delayCard:  document.getElementById('rt-delay-card'),
      elecIdle:   document.getElementById('elec-idle'),
      elecSolved: document.getElementById('elec-solved'),
    };
    this._bindTimingInputs();
  }

  _bindTimingInputs() {}

  render(theta, v0, s_y) {
    this._lastTheta = theta;
    this._lastV0 = v0;
    this._lastSy = s_y;

    const { idle, solved, error, badge, delayCard } = this._els;
    const trig = parseFloat(document.getElementById('in-trig-angle')?.value) || 180;
    const rpm = parseFloat(document.getElementById('in-rpm')?.value) || 10;
    const n = parseFloat(document.getElementById('in-n')?.value) || 0;

    let t_delay = 0;
    if (rpm !== 0) {
      t_delay = ((Math.PI / 180) * (360 - trig + 2 * Math.PI * n)) / (rpm * (Math.PI / 30)) - (0.00954 * theta + 0.335);
    }

    const res = this._engine.compute(t_delay);

    this._hide(idle);
    this._hide(error);
    if (solved) solved.style.display = '';
    if (badge) { badge.textContent = 'OK'; badge.className = 'rt-status-badge rt-badge--ok'; }

    const valEl = document.getElementById('rt-delay-val');
    if (valEl) valEl.textContent = res.t_delay.toFixed(3);

    const { elecIdle, elecSolved } = this._els;
    this._hide(elecIdle);
    if (elecSolved) elecSolved.style.display = '';

    const rEl = document.getElementById('rt-r-val');
    if (rEl) rEl.innerHTML = `${res.r_total.toFixed(2)} <span class="elec-item-unit">kΩ</span>`;
    const vEl = document.getElementById('rt-v-val');
    if (vEl) vEl.innerHTML = `${res.volt.toFixed(2)} <span class="elec-item-unit">V</span>`;

    if (delayCard) {
      delayCard.classList.remove('rt-pop');
      void delayCard.offsetWidth;
      delayCard.classList.add('rt-pop');
    }
  }

  reset() {
    this._lastTheta = null;
    const { idle, solved, error, badge, elecIdle, elecSolved } = this._els;
    if (idle) idle.style.display = '';
    if (solved) solved.style.display = 'none';
    if (error) error.style.display = 'none';
    if (badge) { badge.textContent = '—'; badge.className = 'rt-status-badge'; }
    if (elecIdle) elecIdle.style.display = '';
    if (elecSolved) elecSolved.style.display = 'none';
  }

  _hide(el) { if (el) el.style.display = 'none'; }
}


/* ── 8.5 HistoryController ── */
class HistoryController {
  constructor(app) {
    this._app = app;
    this._history = [];
    this._els = { list: document.getElementById('history-list') };
    this._loadFromStorage();
    this.render();
  }

  _loadFromStorage() {
    try {
      const data = localStorage.getItem('wt_history');
      if (data) this._history = JSON.parse(data);
    } catch { /* ignore */ }
  }

  _saveToStorage() {
    try { localStorage.setItem('wt_history', JSON.stringify(this._history)); }
    catch { /* ignore */ }
  }

  saveCurrent() {
    if (!this._app._lastSol) return;
    const record = {
      id: Date.now().toString(),
      date: new Date().toLocaleTimeString('th-TH', {hour: '2-digit', minute:'2-digit'}),
      y: document.getElementById('in-y-measure')?.value || 0,
      sx: document.getElementById('in-sx')?.value || 0,
      theta: this._app._lastSol.theta.toFixed(1),
      tdelay: document.getElementById('rt-delay-val')?.textContent || '--'
    };
    this._history.unshift(record);
    if (this._history.length > 10) this._history.pop();
    this._saveToStorage();
    this.render();
  }

  loadRecord(id) {
    const record = this._history.find(r => r.id === id);
    if (!record) return;
    const yEl = document.getElementById('in-y-measure');
    if (yEl) { yEl.value = record.y; yEl.dispatchEvent(new Event('change')); }
    setTimeout(() => this._app.calculate(), 50);
  }

  deleteRecord(id) {
    this._history = this._history.filter(r => r.id !== id);
    this._saveToStorage();
    this.render();
  }

  render() {
    if (!this._els.list) return;
    if (this._history.length === 0) {
      this._els.list.innerHTML = '<div class="history-empty">ไม่มีประวัติการบันทึก</div>';
      return;
    }
    this._els.list.innerHTML = this._history.map(r => `
      <div class="hist-item">
        <div class="hist-top">
          <div>
            <div class="hist-title">ระยะ ${r.y} cm</div>
            <div class="hist-sub">เวลาบันทึก: ${r.date}</div>
          </div>
          <div class="hist-actions">
            <button class="hist-btn btn-load" data-id="${r.id}" title="Load"><svg width="12" height="12" viewBox="0 0 12 12" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round"><path d="M2 6h8M6 2l4 4-4 4"/></svg></button>
            <button class="hist-btn btn-del" data-id="${r.id}" title="Delete"><svg width="12" height="12" viewBox="0 0 12 12" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round"><path d="M3 3l6 6M9 3L3 9"/></svg></button>
          </div>
        </div>
        <div class="hist-metrics">
          <span class="hm-badge">θ = ${r.theta}°</span>
          <span class="hm-badge">t = ${r.tdelay}s</span>
        </div>
      </div>
    `).join('');

    this._els.list.querySelectorAll('.btn-load').forEach(btn => {
      btn.addEventListener('click', () => this.loadRecord(btn.dataset.id));
    });
    this._els.list.querySelectorAll('.btn-del').forEach(btn => {
      btn.addEventListener('click', () => this.deleteRecord(btn.dataset.id));
    });
  }
}


/* ── 9. SimulationApp ── */
class SimulationApp {
  constructor() {
    this._physics = new PhysicsEngine(ConfigManager.PHYSICS);
    this._angle = new AngleController(ConfigManager.ANGLE_STEP);
    this._canvas = new CanvasRenderer('traj', ConfigManager.RENDER, ConfigManager.PHYSICS);
    this._ui = new UIRenderer();
    this._input = new InputController();
    this._dirty = new DirtyTracker();
    this._timing = new TimingController(
      new TimingEngine(ConfigManager.PHYSICS.G),
      this._input
    );
    this._lastSol = null;
    this._lastTargetPts = null;
    this._lastMetrics = null;
    this._history = new HistoryController(this);
  }

  get _allInputIds() {
    return [
      ...Object.values(ConfigManager.FIELDS).map(f => f.id),
      ...Object.values(ConfigManager.TIMING_FIELDS).map(f => f.id),
    ];
  }

  init() {
    this._input.setDirtyTracker(this._dirty);
    this._input.bindEvents();
    this._angle.bindButtons();
    this._wireButtons();
    this._wireAngleAdjust();
    this._wireResize();
    this._ui.showIdle();
  }

  calculate() {
    const sx_cm = this._input.getSx();
    const sol = this._physics.solve(sx_cm);
    this._dirty.snapshot(this._allInputIds);

    if (!sol) {
      this._angle.hide();
      this._ui.showNoSolution();
      this._timing.reset();
      const previewAngle = (ConfigManager.PHYSICS.THETA_MIN + ConfigManager.PHYSICS.THETA_MAX) / 2;
      const previewV0 = this._physics.calculateV0ForTarget(previewAngle, sx_cm / 100);
      const previewPts = this._physics.trajectoryPoints(previewAngle, previewV0);
      const previewLandX = this._physics.landingX(previewAngle, previewV0);
      this._canvas.draw({
        targetPts: previewPts, targetAngle: previewAngle,
        currentPts: previewPts, currentAngle: previewAngle,
        currentLandX: previewLandX, xTarget: sx_cm / 100,
      });
      return;
    }

    this._lastSol = sol;
    this._lastTargetPts = this._physics.trajectoryPoints(sol.theta, sol.v0);
    this._lastMetrics = this._physics.computeMetrics(sol.theta, sol.v0);

    this._angle.setTarget(sol.theta);
    this._angle.show();
    this._ui.showSolved(sol.theta, this._lastMetrics);
    this._ui.updateLiveStatus(sol.theta);
    this._timing.render(sol.theta, sol.v0, ConfigManager.PHYSICS.Y_DIFF);
    this._redraw(sol.theta);
    this._history.saveCurrent();
  }

  reset() {
    this._lastSol = null;
    this._lastTargetPts = null;
    this._lastMetrics = null;
    this._input.resetToDefaults();
    this._angle.hide();
    this._ui.showIdle();
    this._canvas.clear();
    this._timing.reset();
    this._dirty.clear();
    this._dirty.markResultsStale(false);
    this._ui.showToast();
  }

  _wireButtons() {
    document.getElementById('btn-calc')?.addEventListener('click', () => this.calculate());
    document.getElementById('btn-reset')?.addEventListener('click', () => this.reset());
  }

  _wireAngleAdjust() {
    this._angle.onAdjust((currentAngle) => {
      if (!this._lastSol) return;
      requestAnimationFrame(() => this._redraw(currentAngle));
    });
  }

  _wireResize() {
    window.addEventListener('resize', () => {
      if (!this._lastSol) return;
      requestAnimationFrame(() => this._redraw(this._angle.current));
    });
  }

  /** Redraw canvas with trajectory band for the given angle. */
  _redraw(currentAngle) {
    const sol = this._lastSol;
    if (!sol || !this._lastTargetPts) return;

    const currentPts = this._physics.trajectoryPoints(currentAngle, sol.v0);
    const currentLandX = this._physics.landingX(currentAngle, sol.v0);

    // Generate trajectory band (±2° around target angle)
    const bandOffset = 2;
    const bandUpper = this._physics.trajectoryPoints(sol.theta + bandOffset, sol.v0);
    const bandLower = this._physics.trajectoryPoints(sol.theta - bandOffset, sol.v0);

    this._canvas.draw({
      targetPts: this._lastTargetPts,
      targetAngle: sol.theta,
      currentPts,
      currentAngle,
      currentLandX,
      xTarget: sol.x,
      bandPtsUpper: bandUpper,
      bandPtsLower: bandLower,
    });

    this._ui.updateLiveStatus(currentAngle);
  }
}


/* ── Bootstrap ── */
document.addEventListener('DOMContentLoaded', () => {
  const app = new SimulationApp();
  app.init();
});