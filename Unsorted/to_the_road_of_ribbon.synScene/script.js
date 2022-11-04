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
//var tAtLast0 = 0;
var bpmTime = 0;
var bassT = 0.0;
var highT = 0.0;
var midT = 0.0;
function update(dt) {

  var bpm = inputs.syn_BPM/4.0;
  bpmcount.updateTime(bpm, dt);

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
  uniforms.smoothTime = (uniforms.smooth_basstime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.2);
  uniforms.smoothTimeB = (uniforms.smooth_hightime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.2);
  uniforms.smoothTimeC = (uniforms.smooth_midtime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.2);


  // cPos.x = cPos.x + (-0.5+Math.random(inputs.syn_BeatTime)*inputs.syn_OnBeat);
  uniforms.camera_pos = cPos;
  uniforms.camera_look = cLook;


  // decimator++;
  // if (decimator%10==0){
  //   print(t);
  // }

}
function transition() {
  //log(5);
}
