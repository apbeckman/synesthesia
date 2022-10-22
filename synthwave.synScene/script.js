function Slicer () {
  this.laggyX = 0.0;
  this.laggyY = 0.0;

  this.x = 0;
  this.y = 0;
  this.angle = 0;
  this.motionAmount = 1.0;


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

var bassT = 0.0;
var highT = 0.0;
var midT = 0.0


function normalize(x, y) {
  var magnitude = Math.sqrt(Math.pow(x,2)+Math.pow(y,2));
  x = x/magnitude;
  y = y/magnitude;
}

Slicer.prototype.updateLaggyMotion = function() {
  this.laggyX = this.laggyX-(this.laggyX-inputs.manual_slice.x)*0.1;
  this.laggyY = this.laggyY-(this.laggyY-inputs.manual_slice.y)*0.1;

}


var decimator = 0;
function update(dt) {

  bassT = bassT + Math.pow(inputs.syn_BassLevel*0.35+inputs.syn_MidLevel*0.15,2.0);
  uniforms.script_bass_time = bassT;
  midT = midT + Math.pow(inputs.syn_MidLevel*0.25+inputs.syn_MidHighLevel*0.25,2.0);
  uniforms.script_bass_time = bassT;

  highT = highT + Math.pow(inputs.syn_MidHighLevel*0.25+inputs.syn_HighLevel*0.25,2.0);
  uniforms.script_high_time = highT;
  uniforms.smooth_basstime = bpmcount.time*0.85+(bassT*0.35);
  uniforms.smooth_midtime = bpmcount.time*0.85+(midT*0.35);
  uniforms.smooth_hightime = bpmcount.time*0.85+(highT*0.35);
  uniforms.smoothTime = (uniforms.smooth_basstime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.125);
  uniforms.smoothTimeB = (uniforms.smooth_hightime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.125);
  uniforms.smoothTimeC = (uniforms.smooth_midtime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.25);

  if(decimator%50==0){
    // console.log(slicer.x);
    // console.log(slicer.y);
    // console.log("next");


    // console.log(slicer.angle);

    // console.log(slicer.motionAmount);

  }
  decimator ++;

}
