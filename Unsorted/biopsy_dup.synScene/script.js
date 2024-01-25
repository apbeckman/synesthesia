function Slicer () {
  this.laggyX = 0.0;
  this.laggyY = 0.0;

  this.x = 0;
  this.y = 0;
  this.angle = 0;
  this.motionAmount = 1.0;


}

function Painter () {
  this.x = 0;
  this.y = 0;
  this.laggyX = 0.0;
  this.laggyY = 0.0;

}

function normalize(x, y) {
  var magnitude = Math.sqrt(Math.pow(x,2)+Math.pow(y,2));
  x = x/magnitude;
  y = y/magnitude;
}

Slicer.prototype.updateLaggyMotion = function() {
  this.laggyX = this.laggyX-(this.laggyX-inputs.manual_slice.x)*0.1;
  this.laggyY = this.laggyY-(this.laggyY-inputs.manual_slice.y)*0.1;

}

Slicer.prototype.getMagDirMotion = function() {
  var tempX = this.laggyX - inputs.manual_slice.x;
  var tempY = this.laggyY - inputs.manual_slice.y;
  var magnitude = Math.sqrt(Math.pow(tempX,2)+Math.pow(tempY,2));

  this.x = tempX/magnitude;
  this.y = tempY/magnitude;

  this.angle = Math.atan(this.y/this.x);

  this.motionAmount = Math.max(Math.abs(inputs.manual_slice.x-this.laggyX),Math.abs(inputs.manual_slice.y-this.laggyY));

}

var slicer = new Slicer();
var painter = new Painter();

var decimator = 0;
function update(dt) {
  slicer.updateLaggyMotion();
  slicer.getMagDirMotion();

  uniforms.slicerDirection = slicer;
  uniforms.slicerAngle = slicer.angle;
  uniforms.slicerMotionAmt = slicer.motionAmount;

  if(decimator%50==0){
    // console.log(slicer.x);
    // console.log(slicer.y);
    // console.log("next");


    // console.log(slicer.angle);

    // console.log(slicer.motionAmount);

  }
  decimator ++;

}
function transition() {
}