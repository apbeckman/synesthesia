function Timer () {
    this.time = 0.0;
    this.updateTime = function(rate, val, dt) {
        this.time = this.time+rate*dt*val;
    }
}

var time_bass = 0;
var time_highs = 0;
var time_fly = 0;
var time_morph = 0;
var bassTimevar = new Timer();
var timevar = new Timer();
var flyTimeVar = new Timer();
var morphTimeVar = new Timer();

function update(dt) {
    var fly_speed = inputs.fly_speed;
    var rate_in = inputs.rate_in;
    var morph_rate = inputs.morph_rate;
    var react = 0;

    //if (inputs.reactive_time > 0.5) react = 1;

    time_bass += 0.05*(inputs.syn_BassLevel + inputs.syn_BassLevel*inputs.syn_BassLevel)*rate_in;
    time_highs += 0.1*(inputs.syn_HighHits + inputs.syn_HighHits*inputs.syn_HighHits)*rate_in;
    time_fly += 0.05*(inputs.Level + inputs.Level*inputs.syn_Level)*fly_speed;
    time_morph += 0.05*(inputs.syn_BassLevel + inputs.syn_BassLevel*inputs.syn_BassLevel)*morph_rate;

    bassTimevar.updateTime(0.8, inputs.reactive_time ? (inputs.syn_BassLevel*2.0+inputs.syn_BassPresence+inputs.syn_BassHits*2.0)*rate_in : rate_in, dt);
    timevar.updateTime(0.4, rate_in, dt);
    flyTimeVar.updateTime(0.4, inputs.reactive_time ?
    (inputs.syn_Level*2.0 + inputs.syn_Presence + inputs.syn_Hits*1.0 + inputs.syn_BassHits*2.0)*fly_speed : fly_speed, dt);
    morphTimeVar.updateTime(0.8, inputs.reactive_morph ? (inputs.syn_BassLevel*2.0+inputs.syn_BassPresence+inputs.syn_BassHits*2.0)*morph_rate : 0, dt);

    uniforms.time_bass = time_bass;
    uniforms.time_highs = time_highs;
    uniforms.time_fly = time_fly;
    uniforms.script_time = timevar.time;
    uniforms.bass_time = bassTimevar.time;
    uniforms.fly_time = flyTimeVar.time;
    uniforms.morph_time = morphTimeVar.time;
}

function transition() {

}
