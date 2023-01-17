var framecount = 0;
var prev_BPMTri = 0;
// factors to compute slowed BPMTri LFOs for, as rationals
var slow_BPMTri = [
  // [1, 8],
  // [1, 4],
  // [1, 3],
  // [1, 2], //double speed
  [1], // should match syn_BPMTri
  // [3, 2],
  [2], //half speed (should match syn_BPMTri2)
  // [3], //third speed etc
  [4]
];
var phi_BPMTri = Array.apply(null, Array(slow_BPMTri.length)).map(function() { return 0 });

function clamp(x, lo, hi) {
  return Math.max(Math.min(x, hi), lo)
}

function scale(x, lo, hi) {
  return x * (hi - lo) + lo;
}

function fract(x) {
  return x - Math.floor(x);
}

function getSpectrum() {
  var highs = inputs.syn_HighHits*inputs.syn_Intensity;
  var highmids = inputs.syn_MidHighLevel;
  var mids = inputs.syn_MidLevel;
  var bass = inputs.syn_BassHits*inputs.syn_Intensity;
  return [bass, mids, highmids, highs];
}

function update(dt) {
  getSpectrum().forEach(function(el, i) {
    uniforms['spectrum_' + i] = el
  })

  var BPMTri = (inputs.syn_BPMTri / inputs.syn_BPMConfidence) || 0;
  var delta_BPMTri = fract(BPMTri - prev_BPMTri);
  slow_BPMTri.forEach(function(r, i) {
    var num = r[0];
    var denom = r[1] || 1;
    phi_BPMTri[i] = fract(phi_BPMTri[i] + delta_BPMTri * denom / num);
    uniforms['BPMTri' + num + (denom == 1 ? '' : '_' + denom)] = phi_BPMTri[i];
  });

  prev_BPMTri = BPMTri;

  framecount++;
}