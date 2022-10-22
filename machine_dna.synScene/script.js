function ApproachRandomSpot () {
  this.goal = 0.0;
  this.currentSpot = 0.0;
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

ApproachRandomSpot.prototype.update = function(amt) {
  this.currentSpot = this.currentSpot + (this.goal-this.currentSpot)*amt;
}

var mover = new ApproachRandomSpot();

var decimator = 0;
var time = 0.0;
var rate = 0.0;
var bpmcount = new BPMCounter();
var bassT = 0.0;
var highT = 0.0;
var midT = 0.0;

function update(dt) {

  time += 0.5*((inputs.syn_BassLevel*inputs.syn_BassLevel*0.1 + inputs.syn_BassPresence*0.2 + 0.01)*inputs.bass_time*1.5);

  uniforms.script_time = time;
  var bpm = inputs.syn_BPM/4.0;
  bpmcount.updateTime(bpm, dt);

  //uniforms.script_time = bpmcount.time;



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
  uniforms.smoothTimeC = (uniforms.smooth_midtime * 0.75 + inputs.TIME * 0.125 + inputs.syn_Time * 0.125);


  // decimator++;
  // if (decimator%10==0){
  //   print(t);
  // }
}

function transition() {
  //log(5);
}
