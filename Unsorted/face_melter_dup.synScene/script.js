function Timer () {

  this.time = 0.0;

}

Timer.prototype.updateTime = function(rate, val, dt) {

  this.time = this.time+rate*dt*val;

}

var timeVar = new Timer();
var bassTimevar = new Timer();
var teet = 0.0;

function update(dt) {

  try {

    if (inputs.teeter_time > 0.5) {
      teet = (Math.floor(Math.floor(inputs.syn_ToggleOnBeat + 0.5) + 0.5) > 0.5 ? -1.0 : 1.0);
    }
    else {
      teet = 1.0;
    }
    
    bassTimevar.updateTime(0.5, (inputs.syn_BassLevel*2.0+inputs.syn_BassPresence+inputs.syn_BassHits*2.0)*inputs.rate*teet*0.5, dt);
    timeVar.updateTime(0.5, inputs.rate*teet, dt);
    uniforms.bass_time = bassTimevar.time;
    uniforms.dynamic_time = (timeVar.time + bassTimevar.time);
    uniforms.fbRotate = inputs._hold > 0.5 ? 0.0 : inputs.syn_BassPresence*0.0005*teet;
    uniforms.zoom = inputs._zoom*0.01;
    uniforms.traceMix = 0.99;

  } 

  catch (e){
    print(e);
  }
}