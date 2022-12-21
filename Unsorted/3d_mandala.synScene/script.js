function Timer () {
  this.time = 0.0;
}

Timer.prototype.updateTime = function(rate, val, dt) {
  this.time = this.time+rate*dt*val;
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

function SmoothCounter () {
  this.oldCount = 0.0;
  this.isGoing = 0.0;
  this.currentValue = 0.0;
}

SmoothCounter.prototype.update = function(dt, newCount, speed) {
  this.currentValue = this.currentValue+(newCount-this.currentValue)*speed;
}

// var bpmcount = new BPMCounter();
// var downrot = new SmoothCounter();
// var rightrot = new SmoothCounter();
var bpmcount = new BPMCounter();

var bassTimevar = new Timer();
var midTimevar = new Timer();
var timevar = new Timer();
var spinTimevar = new Timer();

var decimator = 0;
var t = 0;
var bpmTime = 0;
var bassT = 0.0;
var highT = 0.0;
var midT = 0.0;
var highhits = 0.0;
var basshits = 0.0;
var midhits = 0.0;
var bass = 0.0;
var hits = 0.0;

function update(dt) {
  var bpm = inputs.syn_BPM/4.0;
  bpmcount.updateTime(bpm, dt);

  try {
    bassTimevar.updateTime(0.125, Math.pow((Math.max(0.75, inputs.syn_BassLevel)+inputs.syn_MidLevel*0.75+syn_Intensity*0.25), 2.0)*(inputs.rate_in+inputs.syn_Intensity), dt);

    midTimevar.updateTime(0.1, Math.pow((inputs.syn_BassLevel*0.785+inputs.syn_MidHighLevel*0.9+inputs.syn_MidLevel*1.5+syn_Intensity), 2.0)*(inputs.Pattern_Rot), dt);

    spinTimevar.updateTime(.1,  Math.pow(0.5+inputs.syn_BassLevel*0.25+inputs.syn_MidHighLevel*1.25+inputs.syn_MidLevel*1.75+syn_Intensity,2.0)*(inputs.Pattern_Rot), dt);

    timevar.updateTime(0.4, inputs.rate_in, dt);

    uniforms.script_time = timevar.time;
    uniforms.bass_time = bassTimevar.time;
    uniforms.mid_time = midTimevar.time;
    uniforms.spin_time = spinTimevar.time;

    uniforms.highhits = Math.pow( (inputs.syn_HighLevel*0.5 + inputs.syn_Hits*0.125+inputs.syn_HighHits*0.375)*inputs.syn_Intensity, 2.0);
    uniforms.basshits = Math.pow( (inputs.syn_BassLevel*0.5 + inputs.syn_Level*0.125+inputs.syn_BassHits*0.375)*inputs.syn_Intensity, 2.0);
    uniforms.midhits = Math.pow( (inputs.syn_MidLevel*0.5 + inputs.syn_Level*0.125+inputs.syn_MidHits*0.375)*inputs.syn_Intensity, 2.0);
    uniforms.hits = Math.pow( (inputs.syn_HighLevel*0.5 + inputs.syn_Level*0.125+inputs.syn_Hits*0.375)*inputs.syn_Intensity, 2.0);
  
    uniforms.low = Math.pow( (inputs.syn_BassLevel*0.5 + inputs.syn_MidLevel*0.375+inputs.syn_Level*0.125)*inputs.syn_Intensity, 2.0);
  
    bassT = bassT + Math.pow(inputs.syn_BassLevel*0.45+inputs.syn_MidLevel*0.35+syn_Intensity*0.2,2.0);
    uniforms.script_bass_time = bassT;
    midT = midT + Math.pow(inputs.syn_MidLevel*0.45+inputs.syn_MidHighLevel*0.35+syn_Intensity*0.2,2.0);
    uniforms.script_bass_time = bassT;
  
    highT = highT + Math.pow(inputs.syn_MidHighLevel*0.45+inputs.syn_HighLevel*0.35+syn_Intensity*0.2,2.0);
    uniforms.script_high_time = highT;
    uniforms.smooth_basstime = bpmcount.time*0.85+(bassT*0.35);
    uniforms.smooth_midtime = bpmcount.time*0.85+(midT*0.35);
    uniforms.smooth_hightime = bpmcount.time*0.85+(highT*0.35);
    uniforms.smoothTime = (uniforms.smooth_basstime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.2);
    uniforms.smoothTimeB = (uniforms.smooth_hightime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.2);
    uniforms.smoothTimeC = (uniforms.smooth_midtime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.2);
    
    // downrot.update(dt, inputs.down_rot, 0.1);

    // uniforms.down_rot_scr = downrot.currentValue;

  }  

  // decimator++;
  // if (decimator%10==0){
  //   print(t);
  // }

    catch (e){
    print(e);
  }

}