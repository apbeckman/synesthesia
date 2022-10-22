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
var timevar = new Timer();

var decimator = 0;
var t = 0;
var bpmTime = 0;
var bassT = 0.0;
var highT = 0.0;
var midT = 0.0;

function update(dt) {
  var bpm = inputs.syn_BPM/4.0;
  bpmcount.updateTime(bpm, dt);

  try {
    bassTimevar.updateTime(0.1, (inputs.syn_BassLevel*2.0+inputs.syn_BassPresence+inputs.syn_BassHits*2.0)*inputs.rate_in, dt);
    timevar.updateTime(0.4, inputs.rate_in, dt);

    uniforms.script_time = timevar.time;
    uniforms.bass_time = bassTimevar.time;

    uniforms.smooth_basstime = bpmcount.time*0.425 +timevar.time*0.425+(bassT*0.35);
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