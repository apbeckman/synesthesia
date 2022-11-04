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

var time_fly = 0;
var bassTimevar = new Timer();
var timevar = new Timer();
var flyTimeVar = new Timer();

var decimator = 0;
var t = 0;

function update(dt) {
  var fly_speed = inputs.Fly_Speed;
  try {
    flyTimeVar.updateTime(1.0,
    (inputs.syn_Level*2.0 + inputs.syn_Presence + inputs.syn_Hits*1.0 + inputs.syn_BassHits*4.0)*fly_speed, dt);

    uniforms.fly_time = flyTimeVar.time;

  }

  // decimator++;
  // if (decimator%10==0){
  //   print(t);
  // }

    catch (e){
    print(e);
  }

}
