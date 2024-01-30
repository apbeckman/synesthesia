function Timer() {
  this.time = 0.0;
}

Timer.prototype.updateTime = function (rate, val, dt) {
  this.time = this.time + rate * dt * val;
};

var bassTimevar = new Timer();

var zoomBassTimevar = new Timer();

var glitchBassTimevar = new Counter();

var midTimevar = new Timer();

var midHighTimevar = new Timer();

var nRateMid = new Timer();

var fade = 0.0;

var teet = 0.0;

function Counter() {
  this.count = 0.0;
}

Counter.prototype.updateCount = function (val) {
  this.count += val * 15.0;
};

var bassCount = new Counter();

function update(dt) {
  try {
    bassTimevar.updateTime(
      0.2 * inputs.syn_BassPresence,
      (inputs.syn_BassLevel * 2.0 +
        inputs.syn_BassPresence +
        inputs.syn_BassHits * 2.0) *
        3.0,
      dt
    );

    // setControl("meta/playback_speed", 0.5+inputs.syn_BassLevel*0.5);

    bassCount.updateCount(inputs.syn_BassHits * 0.075);

    uniforms.bass_count = bassCount.count;

    uniforms.bass_time = bassTimevar.time + inputs.TIME * 0.1;

    uniforms.show_media =
      inputs.beat_flash > 0.5
        ? Math.floor(inputs.syn_OnBeat * 10.5) * 0.1
        : inputs.level_media > 0.5
        ? inputs.syn_MidLevel * inputs.syn_Presence * 0.25
        : inputs._show_media > 0.5
        ? 0.2
        : inputs.hint_media > 0.5
        ? inputs.syn_Presence * inputs.syn_Intensity * 0.125
        : 0.0; //+(inputs.syn_BassLevel + inputs.syn_BassPresence)/2.0*0.0625 + inputs.syn_MidLevel*0.03 : 0.0;

    uniforms.traceMix = 0.975; // + 0.0125*inputs._flat;// (1.0 - fade) + 0.001*fade;

    uniforms.tS = 0.05 * inputs.syn_BassLevel + 0.05 * inputs.syn_BassPresence;
  } catch (e) {
    print(e);
  }
}
