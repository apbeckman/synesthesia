function Timer () {
    this.time = 0.0;
    this.updateTime = function(rate, val, dt) {
        this.time = this.time+rate*dt*val;
    }
}

var time_bass = 0;
var time_highs = 0;
var time_levels = 0;
var bassTimevar = new Timer();
var timevar = new Timer();
var levelTimevar = new Timer();

function update(dt) {
    var rate_in = inputs.rate_in;
    var fly_speed = inputs.fly_speed;
    /*if (inputs.zoom >= 1 || inputs.zoom <= -1) {
        rate_in += (inputs.zoom * 0.1);
    }*/
    time_bass += 0.05*(inputs.syn_BassLevel + inputs.syn_BassLevel*inputs.syn_BassLevel)*rate_in;
    time_highs += 0.1*(inputs.syn_HighHits + inputs.syn_HighHits*inputs.syn_HighHits)*rate_in;
    time_levels += 0.05*(inputs.syn_Level + inputs.syn_Level*inputs.syn_Level)*fly_speed;

    bassTimevar.updateTime(0.8, inputs.reactive_time ? (inputs.syn_BassLevel*2.0+inputs.syn_BassPresence+inputs.syn_BassHits*2.0)*rate_in : rate_in, dt);
    timevar.updateTime(0.4, rate_in, dt);
    levelTimevar.updateTime(0.8, inputs.reactive_time ? (inputs.syn_Level*2.0+inputs.syn_Presence+inputs.syn_Hits*2.0+inputs.syn_BassHits*2.0)*rate_in : rate_in, dt);

    uniforms.time_bass = time_bass;
    uniforms.time_highs = time_highs;
    uniforms.time_levels = time_levels;
    uniforms.script_time = timevar.time;
    uniforms.bass_time = bassTimevar.time;
    uniforms.fly_time = levelTimevar.time;
}

function transition() {

}
