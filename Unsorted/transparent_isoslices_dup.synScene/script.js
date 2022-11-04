var t = 0;
var tMorph = 0;

function update(dt) {
  t = t + inputs.rate_in + inputs.syn_BassLevel*inputs.syn_BassLevel*inputs.bass_rate;
  tMorph = tMorph + inputs.syn_BassLevel*inputs.syn_BassLevel*inputs.bass_morph*inputs.bass_morph;

  uniforms.script_time = t;
  uniforms.morph_time = tMorph;

}