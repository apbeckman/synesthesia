function Timer () {
    this.time = 0.0;
    this.updateTime = function(rate, val, dt) {
        this.time = this.time+rate*dt*val;
    }
}
function BPMCounter () {
    this.time = 0.0;
    this.count = 0.0;
    this.timeWithinBeat = 0.0;
    this.didIncrement = 0.0;
  }
  
  BPMCounter.prototype.updateTime = function(bpm, dt) {
    this.didIncrement = 0.0;
    var amountToStepThroughBeat = bpm*dt/60.0;
    this.time = this.time+amountToStepThroughBeat;
    if(this.count != Math.floor(this.time)){
      this.count = Math.floor(this.time);
      this.didIncrement = 1.0;
    };
    this.timeWithinBeat = this.time-this.count;
  }
  
var time_bass = 0;
var time_highs = 0;
var time_levels = 0;
var bassTimevar = new Timer();
var timevar = new Timer();
var levelTimevar = new Timer();
var bpmcount = new BPMCounter();
//var cPos = new CameraPos();
//var cLook = new CameraLook();

var decimator = 0;
var tAtLast0 = 0;
var bpmTime = 0;
var bassT = 0.0;
var highT = 0.0;
var midT = 0.0;
var highhits = 0.0;

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

    bassT = bassT + Math.pow(inputs.syn_BassLevel*0.35+inputs.syn_MidLevel*0.15,2.0);
    uniforms.script_bass_time = bassT;
    midT = midT + Math.pow(inputs.syn_MidLevel*0.35+inputs.syn_MidHighLevel*0.15,2.0);
    uniforms.script_bass_time = bassT;
  
    highT = highT + Math.pow(inputs.syn_MidHighLevel*0.35+inputs.syn_HighLevel*0.15,2.0);
    uniforms.script_high_time = highT;
    uniforms.smooth_basstime = bpmcount.time*0.85+(bassT*0.35);
    uniforms.smooth_midtime = bpmcount.time*0.85+(midT*0.35);
    uniforms.smooth_hightime = bpmcount.time*0.85+(highT*0.35);
    uniforms.smoothTime = (uniforms.smooth_basstime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.125);
    uniforms.smoothTimeB = (uniforms.smooth_hightime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.125);
    uniforms.smoothTimeC = (uniforms.smooth_midtime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.25);
  
  
    uniforms.time_bass = time_bass;
    uniforms.time_highs = time_highs;
    uniforms.time_levels = time_levels;
    uniforms.script_time = timevar.time;
    uniforms.bass_time = bassTimevar.time;
    uniforms.fly_time = levelTimevar.time;
}

function transition() {

}
