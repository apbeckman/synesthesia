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


function TimeVar () {
  this.rate = 0.0;
  this.time = 0.0;
}

TimeVar.prototype.updateTime = function(val, dt) {
  this.time = this.time+val*this.rate*(dt*60);
}

function Accumulator () {
  this.value = 0.0;
  this.decay = 0.9;
}

Accumulator.prototype.update = function(val, dt){
  this.value += (val - this.value) * (this.decay * dt * 60.0);
  if (this.value > 1.0){
    this.value = 1.0;
  }
}
var bpmcount = new BPMCounter();

var timevar = new TimeVar();
var accumulator = new Accumulator();
var bpmTime = 0;
var bassT = 0.0;
var highT = 0.0;
var midT = 0.0;

var decimator = 0;
function update(dt) {
  var bpm = inputs.syn_BPM/4.0;
  timevar.rate = inputs.speed;
  accumulator.update(Math.pow(inputs.syn_Level*0.125+syn_BassLevel*0.35+syn_MidLevel*0.25,2.0), dt);
  timevar.updateTime(accumulator.value, dt);
  //uniforms.script_time = timevar.time;
  uniforms.script_time = bpmcount.time;



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

  if(decimator%50==0){
    // console.log(timevar.time);
  }
  decimator ++;

}
function transition() {
}