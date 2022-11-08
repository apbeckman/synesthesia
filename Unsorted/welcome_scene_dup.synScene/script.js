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

var decimator = 0;
var tAtLast0 = 0;
var bpmTime = 0;
var bassT = 0.0;

function update(dt) {

  var bpm = inputs.syn_BPM;
  bpmcount.updateTime(bpm, dt*Math.pow(inputs.syn_FadeInOut,2.0));

  uniforms.script_time = (bpmcount.count + Math.pow(bpmcount.timeWithinBeat, 0.2));

  // if (bpmcount.didIncrement == 1.0){
  //   tAtLast0 = bpmTime;
  // }
  // bpmTime = tAtLast0;
  // bpmTime += (1. - Math.exp(-bpmcount.timeWithinBeat*3.))*inputs.amount_to_step;
  // uniforms.script_time = bpmTime;

  bassT = bassT + Math.pow(inputs.syn_BassLevel*0.5,2.0);
  uniforms.script_bass_time = bassT;



  // decimator++;
  // if (decimator%10==0){
  //   print(t);
  // }

}
function transition() {
  //log(5);
}
