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

function CameraPos () {
  this.x = 0.0;
  this.y = 0.0;
  this.z = 0.0;
}

function CameraLook () {
  this.x = 0.0;
  this.y = 0.0;
  this.z = 0.0;
}


var bpmcount = new BPMCounter();
var cPos = new CameraPos();
var cLook = new CameraLook();

var decimator = 0;
var tAtLast0 = 0;
var bpmTime = 0;
var bassT = 0.0;
var highT = 0.0;
var midT = 0.0;
var highhits = 0.0;
var basshits = 0.0;
var midhits = 0.0;

function update(dt) {

  var bpm = inputs.syn_BPM/4.0;
  bpmcount.updateTime(bpm, dt);

  uniforms.script_time = bpmcount.time;

  // if (bpmcount.didIncrement == 1.0){
  //   tAtLast0 = bpmTime;
  // }
  // bpmTime = tAtLast0;
  // bpmTime += (1. - Math.exp(-bpmcount.timeWithinBeat*3.))*inputs.amount_to_step;
  // uniforms.script_time = bpmTime;
  uniforms.highhits = 0.75*Math.pow( (inputs.syn_HighLevel*0.45 + inputs.syn_Level*0.125+inputs.syn_HighHits*0.25), 2.0)*inputs.syn_Intensity;
  uniforms.basshits = 0.75*Math.pow( (inputs.syn_BassLevel*0.45 + inputs.syn_Level*0.125+inputs.syn_BassHits*0.25), 2.0)*inputs.syn_Intensity;
  uniforms.midhits = 0.75*Math.pow( (inputs.syn_MidLevel*0.45 + inputs.syn_Level*0.125+inputs.syn_MidHits*0.25), 2.0)*inputs.syn_Intensity;
  uniforms.level = 
  bassT = bassT + Math.pow(inputs.syn_BassLevel*0.35+inputs.syn_MidLevel*0.15,2.0);
  uniforms.bassT = bassT;
  uniforms.script_bass_time = bassT;
  midT = midT + Math.pow(inputs.syn_MidLevel*0.25+inputs.syn_MidHighLevel*0.25,2.0);
  uniforms.script_bass_time = bassT;

  highT = highT + Math.pow(inputs.syn_MidHighLevel*0.25+inputs.syn_HighLevel*0.25,2.0);
  uniforms.script_high_time = highT;
  uniforms.smooth_basstime = bpmcount.time*0.85+(bassT*0.35);
  uniforms.smooth_midtime = bpmcount.time*0.85+(midT*0.35);
  uniforms.smooth_hightime = bpmcount.time*0.85+(highT*0.35);
  uniforms.smoothTime = (uniforms.smooth_basstime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.2);
  uniforms.smoothTimeB = (uniforms.smooth_hightime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.2);
  uniforms.smoothTimeC = (uniforms.smooth_midtime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.2);
  
  // cPos.x = cPos.x + (-0.5+Math.random(inputs.syn_BeatTime)*inputs.syn_OnBeat);
  // decimator++;
  // if (decimator%10==0){
  //   print(t);
  // }

}
function transition() {
  //log(5);
}
