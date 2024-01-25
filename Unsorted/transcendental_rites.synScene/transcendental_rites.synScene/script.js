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

  this.currentValue = this.currentValue+newCount;

}

var stripeTimeVar = new Timer();
var timevar = new Timer();
var bassTimevar = new Timer();
var pMidTimevar = new Timer();
var pBassTimevar = new Timer();
var pMidHighTimevar = new Timer();
var dynTime = 0.0;

var timevarGeo = new Timer();
var bassTimevarGeo = new Timer();
var dynTimeGeo = 0.0;
var randomOnBeat = 0.0;

function update(dt) {

  try {

    stripeTimeVar.updateTime(0.2, (inputs.syn_MidHits*3.0)*inputs.beat_stripe, dt);
    
    uniforms.beatStripe = stripeTimeVar.time;

    bassTimevar.updateTime(0.1, (inputs.syn_BassLevel*2.0+inputs.syn_BassPresence+inputs.syn_BassHits*2.0)*inputs.spin_rate, dt);

    
    randomOnBeat = inputs.syn_RandomOnBeat;
    timevar.updateTime(0.4, inputs.spin_rate, dt);

    dynTime = inputs.reactive_time > 0.5 ? bassTimevar.time : timevar.time;

    uniforms.script_time = timevar.time;

    uniforms.bass_time = bassTimevar.time;

    uniforms.dynamic_time = dynTime;


    bassTimevarGeo.updateTime(0.1, (inputs.syn_BassLevel*2.0+inputs.syn_BassPresence+inputs.syn_BassHits*2.0)*inputs.geo_rate, dt);

    timevarGeo.updateTime(0.4, inputs.geo_rate, dt);

    dynTimeGeo = inputs.geo_reactive_time > 0.5 ? bassTimevarGeo.time : timevarGeo.time;

    uniforms.geo_script_time = timevarGeo.time;

    uniforms.geo_bass_time = bassTimevarGeo.time;

    uniforms.geo_dynamic_time = dynTimeGeo;
    
    pMidTimevar.updateTime(0.1, (inputs.syn_MidLevel*2.0+inputs.syn_MidPresence+inputs.syn_MidHits*2.0)*5.0, dt);

    pBassTimevar.updateTime(0.1, (inputs.syn_MidHighLevel*2.0+inputs.syn_MidHighPresence+inputs.syn_MidHighHits*2.0)*-5.0, dt);

    pMidHighTimevar.updateTime(0.1, (inputs.syn_BassLevel*2.0+inputs.syn_BassPresence+inputs.syn_BassHits*2.0)*5.0, dt);


    uniforms.pBass = pBassTimevar.time;
    uniforms.pMid = pMidTimevar.time;
    uniforms.pMidHigh = pMidHighTimevar.time;

    uniforms.abs_sin_time_8th_speed = Math.abs(Math.sin(dynTime*0.125));


    uniforms.edgeMix = inputs.audio_edge > 0.5 ? inputs.syn_BassHits*0.5 : inputs.edge_mix;

    uniforms.traceMix = (inputs.feedback_on - inputs.fb_fade) + 0.07*inputs.fb_fade; 
    uniforms.trace_zoom = inputs.trace_on > 0.0 || inputs.keep_trace > 0.0 ? inputs.zoom_trace : 0.0;
    uniforms.inTube = inputs.rand_mode > 0.5 ? Math.floor(inputs.syn_ToggleOnBeat+0.5)*39.0 : inputs.in_da_tube;
    uniforms.fov_correct = inputs.in_da_tube > 38.0 ? 1.0 : 0.0;

  } 
  catch (e){

    console.log(JSON.stringify(e));

  }







}