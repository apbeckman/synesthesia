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

// var bpmcount = new BPMCounter();
// var downrot = new SmoothCounter();
// var rightrot = new SmoothCounter();

var bassTimevar = new Timer();
var timevar = new Timer();

var decimator = 0;
var t = 0;

function update(dt) {

  try {
    bassTimevar.updateTime(0.1, (inputs.syn_BassLevel*2.0+inputs.syn_BassPresence+inputs.syn_BassHits*2.0)*inputs.rate_in, dt);
    timevar.updateTime(0.4, inputs.rate_in, dt);

    uniforms.script_time = timevar.time;
    uniforms.bass_time = bassTimevar.time;

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