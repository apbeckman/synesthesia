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
  


var bpmcount = new BPMCounter();

var bpmTime = 0;
var bassT = 0.0;
var highT = 0.0;
var midT = 0.0;
var highhits = 0.0;
var basshits = 0.0;
var time_bass = 0;
var time_highs = 0;
var time_fly = 0;
var time_morph = 0;
var bassTimevar = new Timer();
var timevar = new Timer();
var flyTimeVar = new Timer();
var morphTimeVar = new Timer();
var hueTimeVar = new Timer();
var lightTimeVar = new Timer();

function update(dt) {
    var fly_speed = inputs.fly_speed;
    var rate_in = inputs.rate_in;
    var morph_rate = inputs.morph_rate;    var react = 0;

    //if (inputs.reactive_time > 0.5) react = 1;
    uniforms.highhits = 0.65*Math.pow(inputs.syn_HighLevel*0.25 + (inputs.syn_Level*0.25+inputs.syn_HighHits*0.35), 2.0)*inputs.syn_Intensity;
    uniforms.basshits = 0.65*Math.pow(inputs.syn_BassLevel*0.25 + (inputs.syn_Level*0.25+inputs.syn_BassHits*0.35), 2.0)*inputs.syn_Intensity;
  
    bassT = bassT + Math.pow(inputs.syn_BassLevel*0.35+inputs.syn_MidLevel*0.15+syn_Intensity*0.2,2.0);
    uniforms.script_bass_time = bassT;
    midT = midT + Math.pow(inputs.syn_MidLevel*0.25+inputs.syn_MidHighLevel*0.25+syn_Intensity*0.2,2.0);
    uniforms.script_bass_time = bassT;
  
    highT = highT + Math.pow(inputs.syn_MidHighLevel*0.25+inputs.syn_HighLevel*0.25,2.0);
    uniforms.script_high_time = highT;
    uniforms.smooth_basstime = bpmcount.time*0.85+(bassT*0.35);
    uniforms.smooth_midtime = bpmcount.time*0.85+(midT*0.35);
    uniforms.smooth_hightime = bpmcount.time*0.85+(highT*0.35);
    uniforms.smoothTime = (uniforms.smooth_basstime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.2);
    uniforms.smoothTimeB = (uniforms.smooth_hightime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.2);
    uniforms.smoothTimeC = (uniforms.smooth_midtime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.2);
  
    time_bass += 0.05*(bassT)*rate_in;
    time_highs += 0.1*(inputs.syn_HighHits + inputs.syn_HighHits*inputs.syn_HighHits)*rate_in;
    time_fly += 0.05*(inputs.Level + inputs.Level*inputs.syn_Level)*fly_speed;
    time_morph += 0.05*(midT+0.1 )*morph_rate;

    bassTimevar.updateTime(0.8, inputs.reactive_time ? (basshits)*rate_in : rate_in, dt);
    timevar.updateTime(0.4, rate_in, dt);
    flyTimeVar.updateTime(0.6, inputs.reactive_time ? (inputs.syn_Level*2.0 + inputs.syn_Presence + inputs.syn_Hits*1.0 + inputs.syn_BassHits*2.0)*fly_speed : fly_speed, dt);
    morphTimeVar.updateTime(0.5, inputs.reactive_morph ? (inputs.syn_BassLevel*2.0+inputs.syn_BassPresence+inputs.syn_BassHits*2.0)*morph_rate : 0, dt);
    hueTimeVar.updateTime(0.3, inputs.reactive_hue ? (inputs.syn_BassLevel*2.0+inputs.syn_BassPresence+inputs.syn_BassHits*4.0)*inputs.hue_rate : 0, dt);
    lightTimeVar.updateTime(0.8, inputs.reactive_lights ? (inputs.syn_BassLevel*2.0+inputs.syn_BassPresence+inputs.syn_BassHits*8.0)*inputs.light_rate : 0, dt);

    uniforms.time_bass = time_bass;
    uniforms.time_highs = time_highs;
    uniforms.time_fly = time_fly;
    uniforms.script_time = timevar.time;
    uniforms.bass_time = bassTimevar.time;
    uniforms.fly_time = flyTimeVar.time;
    uniforms.morph_time = morphTimeVar.time;
    uniforms.hue_time = hueTimeVar.time;
    uniforms.light_time = lightTimeVar.time;
}

function transition() {

}
