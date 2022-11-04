// function ApproachRandomSpot () {
//   this.goal = 0.0;
//   this.currentSpot = 0.0;
// }

// ApproachRandomSpot.prototype.update = function(amt) {
//   this.currentSpot = this.currentSpot + (this.goal-this.currentSpot)*amt;
// }

// var mover = new ApproachRandomSpot();

var decimator = 0;
var time_forward = 0;
var time_color = 0;
var time_bass = 0;
var time_highs = 0;

function update(dt) {

  // if (inputs.syn_OnBeat > 0.9){
  //   mover.goal = Math.random();
  // }
  // // mover.goal = 0.5+0.5*Math.sin(inputs.TIME);
  // mover.update(0.02);
  // uniforms.random_spot = mover.currentSpot;

  time_forward += 0.005*(inputs.syn_BassLevel*inputs.syn_BassLevel*1.5 + inputs.syn_Presence*0.3)*inputs.rate_forward;
  time_color += 0.01*(inputs.syn_HighHits + inputs.syn_Level*inputs.syn_Level)*inputs.rate_color;
  time_bass += 0.05*(inputs.syn_BassLevel + inputs.syn_BassLevel*inputs.syn_BassLevel)*inputs.smooth_rotate;
  time_highs += 0.1*(inputs.syn_HighHits + inputs.syn_HighHits*inputs.syn_HighHits)*inputs.rate_color;

  uniforms.time_forward = time_forward;
  uniforms.time_color = time_color;
  uniforms.time_bass = time_bass;
  uniforms.time_highs = time_highs;

  // decimator++;
  // if (decimator%10==0){
  //   print(time_bass);
  // }
}

function transition() {
  //log(5);
}
