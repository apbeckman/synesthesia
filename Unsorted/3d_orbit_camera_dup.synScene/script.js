/* LastClick: 
*  store an x and y value from each of the mouse buttons
*  used to calculate click drag distance
*/
function LastClick() {
    this.x = -100
    this.y = -100
}
LastClick.prototype.click = function() {
    this.x = _muvc.x
    this.y = _muvc.y
}
LastClick.prototype.release = function() {
    this.x = -100
    this.y = -100
}
var last_clickL = new LastClick()
var last_clickR = new LastClick()
/* Smooth3:
*  a generic vec3 interface with smoothing applied
*  used for camera position and zoom
*/
function Smooth3() {
    this.x = 0
    this.y = 0
    this.z = 0
    this.goal_x = 0
    this.goal_y = 0
    this.goal_z = 0
}
// helper function used for random values from -2PI to 2PI
function rand2PIsigned() {return Math.random() * Math.PI * 2 - Math.PI;}
// smoothly move to the x,y,z input values
Smooth3.prototype.smooth = function(x,y,z) {
    // accept an array as input for convenience
    if (typeof x == 'object') {
        z = x[2]
        y = x[1]
        x = x[0]
    }
    // fallback to random numbers from -2pi to 2pi
    x = typeof x == 'number' ? x : rand2PIsigned()
    y = typeof y == 'number' ? y : rand2PIsigned()
    z = typeof z == 'number' ? z : rand2PIsigned()
    this.goal_x = x
    this.goal_y = y
    this.goal_z = z
}
// instantly move the x,y,z input values
Smooth3.prototype.instant = function(x,y,z) {
    // accept an array as input for convenience
    if (typeof x == 'object') {
        z = x[2]
        y = x[1]
        x = x[0]
    }
    // fallback to random numbers from -2pi to 2pi
    x = typeof x == 'number' ? x : rand2PIsigned()
    y = typeof y == 'number' ? y : rand2PIsigned()
    z = typeof z == 'number' ? z : rand2PIsigned()
    this.x = x
    this.y = y
    this.z = z
    this.smooth(x,y,z)
}
// process the smoothed movement, applied each frame
Smooth3.prototype.tick = function() {
    this.x += (this.goal_x - this.x) * 0.5
    this.y += (this.goal_y - this.y) * 0.5
    this.z += (this.goal_z - this.z) * 0.5
}
// create our camera Smooth3's and send them to the GLSL as vec3 uniforms
var cam = new Smooth3()
uniforms.cam = cam
// set up some global variables
var deltaL = {x:0,y:0}
var deltaR = {x:0,y:0}
var deltaM = {x:0,y:0}
var drag_timeout = 0;
/* setup():
*
*  Synesthesia runs this function once at the start of the scene.
*  Here we're using it to set the default camera position and register some 
*  callbacks - javascript functions to call when actions happen in Synesthesia.
*/
function setup() {
    // copy default values from control panel so we don't have to define them twice
    cam.instant(cam_xy.x,cam_xy.y,cam_zoom)
    // call the viewFront() function every time the view_front bang control is pressed
    onOffToOn('view_front','viewFront')
    // and the others
    onOffToOn('view_top','viewTop')
    onOffToOn('view_right','viewRight')
    onOffToOn('view_back','viewBack')
    onOffToOn('view_left','viewLeft')
    onOffToOn('view_bottom','viewBottom')
}
// Here are the callback functions we register in setup
function viewFront() {setControl('cam_xy',[0,0])}
function viewTop() {setControl('cam_xy',[0,Math.PI/2])}
function viewBottom() {setControl('cam_xy',[0,-Math.PI/2])}
function viewRight() {setControl('cam_xy',[Math.PI/2,0])}
function viewLeft() {setControl('cam_xy',[-Math.PI/2,0])}
function viewBack() {setControl('cam_xy',[-Math.PI,0])}
/* update():
*
*  This is the main loop that Synesthesia will run before every frame is rendered.
*  Here we will process mouse events and control panel events and translate 
*  those into smoothed camera position and rotation variables. We'll then be able
*  to access those in main.glsl to use in our raymarcher.
*/
function update(dt) {
    // get click status as booleans for convenience
    var clickL = _click.x > 0.5
    var clickR = _click.y > 0.5
    // LEFT CLICK  
    if (clickL && last_clickL.x < -50)
        last_clickL.click()
    else if (!clickL && last_clickL.x > -50)
        last_clickL.release()
    if (clickL) { // smoothed drag
        deltaL.x = _muvc.x - last_clickL.x
        deltaL.y = _muvc.y - last_clickL.y
    } else if (deltaL.x!=0||deltaL.y!=0) { // apply momentum decay on release
        deltaL.x *= 0.95
        deltaL.y *= 0.95
        // zero out very low momentum to stop unnecessary calculations
        if (Math.abs(deltaL.x)<0.01) deltaL.x=0;
        if (Math.abs(deltaL.y)<0.01) deltaL.y=0;
    }
    // apply the left click drag information
    // left click x:
    cam.goal_x += deltaL.x * 0.08
    // left click y:
    cam.goal_y += deltaL.y * 0.08
    
    // RIGHT CLICK
    if (clickR && last_clickR.x < -50)
        last_clickR.click()
    else if (!clickR && last_clickR.x > -50)
        last_clickR.release()
    if (clickR) {
        deltaR.y = _muvc.y - last_clickR.y
    } else if (deltaR.x!=0||deltaR.y!=0) { // apply momentum decay on release
        deltaR.y *= 0.95
        // zero out very low momentum to stop unnecessary calculations
        if (Math.abs(deltaR.y)<0.01) deltaR.y=0
    }
    
    cam.goal_z += deltaR.y * 0.08
    
    // limit and wrap XY to range -PI to PI
    if (cam.x < -Math.PI) {
        cam.x += 2*Math.PI
        cam.goal_x += 2*Math.PI
    }
    if (cam.x > Math.PI) {
        cam.x -= 2*Math.PI
        cam.goal_x -= 2*Math.PI
    }
    // clamp Y from negative half pi to positive half pi
    // limits vertical rotation from top to bottom pole
    cam.goal_y = Math.min(Math.max(cam.goal_y, -Math.PI/2), Math.PI/2)
    
    // figure out if we're using mouse or control panel controls currently
    var dragging = deltaL.x!=0||deltaL.y!=0||deltaR.y!=0
    if (dragging) drag_timeout = 20
    if (drag_timeout > 0) {
        // if we were dragging within the last 20 frames, process the drag data
        if (!dragging) drag_timeout--
        setControl('cam_xy',[cam.x,cam.y]);
        setControl('cam_zoom',cam.z);
    } else {
        // otherwise get updates from the Synesthesia scene controls
        cam.smooth(
            (inputs.cam_xy.x + cam.x*19) / 20.,
            (inputs.cam_xy.y + cam.y*19) / 20.,
            (inputs.cam_zoom + cam.z*19) / 20.
        )
    }
    // update smooth3s
    cam.tick()
}