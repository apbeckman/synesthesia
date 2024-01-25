function BPMCounter () {
  this.time = 0.0;
  this.count = 0.0;
  this.timeWithinBeat = 0.0;
  this.didIncrement = 0.0;
}

BPMCounter.prototype.updateTime = function(bpm, dt) {
  this.didIncrement = 0.0;
  var amountToStepThroughBeat = 0.25*bpm*dt/60.0;
  this.time = this.time+amountToStepThroughBeat;
  if(this.count !== Math.floor(this.time)){
    this.count = Math.floor(this.time);
    this.didIncrement = 1.0;
  }
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

var bpmcount = new BPMCounter();
var firepulse = new SmoothCounter();
// var rightrot = new SmoothCounter();

var decimator = 0;
var t = 0;

function update(dt) {
  bpmcount.updateTime(inputs.syn_BPM, dt);
  uniforms.bpmTime = bpmcount.time;
  firepulse.update(dt, inputs.fire_pulse, 0.01);

  uniforms.fire_pulse_scr = firepulse.currentValue;

  // decimator++;
  // if (decimator%20==0){
  //   print(uniforms.fire_pulse_scr);
  // }

}
function transition() {
  //log(5);
}
