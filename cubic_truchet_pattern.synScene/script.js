function Timer () {
    this.time = 0.0;
    this.updateTime = function(rate, val, dt) {
        this.time = this.time+rate*dt*val;
    }
}

var time_fly = 0;
var time_morph = 0;
var flyTimeVar = new Timer();
var morphTimeVar = new Timer();

function update(dt) {
    var fly_speed = inputs.Fly_Speed;
    var morph_rate = inputs.morph_rate;

    //if (inputs.reactive_time > 0.5) react = 1;

    time_fly += 0.05*(inputs.Level + inputs.Level*inputs.syn_Level)*fly_speed;
    time_morph += 0.05*(inputs.syn_BassLevel + inputs.syn_BassLevel*inputs.syn_BassLevel)*morph_rate;

    flyTimeVar.updateTime(0.4, inputs.reactive_time ?
    (inputs.syn_Level*2.0 + inputs.syn_Presence + inputs.syn_Hits*1.0 + inputs.syn_BassHits*2.0)*fly_speed : fly_speed, dt);
    morphTimeVar.updateTime(0.8, inputs.reactive_morph ? (inputs.syn_BassLevel*2.0+inputs.syn_BassPresence+inputs.syn_BassHits*2.0)*morph_rate : 0, dt);

    uniforms.time_fly = time_fly;
    uniforms.fly_time = flyTimeVar.time;
    uniforms.morph_time = morphTimeVar.time;
}

function transition() {

}
