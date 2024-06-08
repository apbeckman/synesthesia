//Vial Spiller - Psybernautics - 2022
function Timer () {

  this.time = 0.0;

}

Timer.prototype.updateTime = function(rate, val, dt) {

  this.time = this.time+rate*dt*val;
  

}
var bassTimevar = new Timer();
var fade = 0.0;
var teet = 0.0;
var gone = 0.0;

function Counter () {

  this.count = 0.0;

}

Counter.prototype.updateCount = function(val) {

  this.count += val*15.0;

}

var bassCount = new Counter();
var media = 0.0;
function update(dt) {

  try {
    bassTimevar.updateTime(0.1+0.1*_size+0.025*inputs.syn_BassPresence, (inputs.syn_BassLevel*2.0+inputs.syn_BassPresence+inputs.syn_BassHits*2.0), dt);
    


    bassCount.updateCount(syn_BassHits*(0.075+_size*0.005));
    uniforms.bass_count = bassCount.count;

    uniforms.bass_time = bassTimevar.time;

    media = inputs.flash_feedback > 0.001 ? Math.abs(0.0-inputs.flash_feedback) : (inputs.mosh_on > 0.5 ?  ( inputs.drum_hits > 0.5 ? (syn_BassHits>0.5 ? 1.0 : 1.0-syn_Level*0.125) : inputs.syn_Presence*(1.2-inputs.syn_BassHits*0.25)) : inputs._show_feedback);
    uniforms.show_media = alt_fb > 0.5 || alt_auto > 0.5 && mosh_on < 0.5 ? 0.05*syn_Hits : (1.0 - (media > 1.0 ? 1.0 : media));
    uniforms.cs = inputs._cartoon_style;
    uniforms.csf = inputs._cartoon_style;
    fade = (1.0 - inputs.syn_BassPresence)*(0.000625);
    uniforms.traceMix = 1.0;
    uniforms.tS = 0.1*inputs.syn_Presence; 
    uniforms.alt_av = syn_Presence;
  } 

  catch (e){
    print(e);
  }
}