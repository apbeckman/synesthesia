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
function Timer () {
  this.time = 0.0;
}

Timer.prototype.updateTime = function(rate, val, dt) {
  this.time = this.time+rate*dt*val;
}
function SmoothCounter () {
  this.oldCount = 0.0;
  this.isGoing = 0.0;
  this.currentValue = 0.0;
}

SmoothCounter.prototype.update = function(dt, newCount, speed) {
  this.currentValue = this.currentValue+(newCount-this.currentValue)*speed;
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

var bassTimevar = new Timer();
var bpmcount = new BPMCounter();
var cPos = new CameraPos();
var cLook = new CameraLook();
var clicked = false;
var decimator = 0;
var tAtLast0 = 0;
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



    if (_click.x > 0.5) {
        setControl('manual_position',[_muv.x*2-1, _muv.y*2-1])
        setControl('manual_impulse',1);
        clicked = true;
    } else if (clicked) {
        setControl('manual_impulse',0);
        clicked = false;
    }

  var bpm = inputs.syn_BPM/4.0;
  bpmcount.updateTime(bpm, dt);
  bassTimevar.updateTime(1.,  Math.pow(0.5+inputs.syn_BassLevel+inputs.syn_MidLevel*0.75+syn_Intensity*0.75,2.0)*(inputs.rate_in), dt);
  //timevar.updateTime(0.4, inputs.rate_in, dt);
  //timevar.updateTime(0.4, inputs.rate_in, dt);


  uniforms.script_time = bpmcount.time;
  uniforms.bass_time = bassTimevar.time;

  // if (bpmcount.didIncrement == 1.0){
  //   tAtLast0 = bpmTime;
  // }
  // bpmTime = tAtLast0;
  // bpmTime += (1. - Math.exp(-bpmcount.timeWithinBeat*3.))*inputs.amount_to_step;
  // uniforms.script_time = bpmTime;
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
  
  // cPos.x = cPos.x + (-0.5+Math.random(inputs.syn_BeatTime)*inputs.syn_OnBeat);
  // decimator++;
  // if (decimator%10==0){
  //   print(t);
  // }

}
function transition() {
  //log(5);
}
