function LastMouse() {

    this.x = 0

    this.y = 0

    this.lastx = 0

    this.lasty = 0

}



LastMouse.prototype.tick = function() {

    this.x = this.lastx

    this.y = this.lasty

    this.lastx = _muvc.x

    this.lasty = _muvc.y

}



lastmouse = new LastMouse

uniforms.lastmouse = lastmouse



var timer = 0;



function update() {

    timer++

    if (timer > 60) {

        // console.log(_muvc.x + " " + _muvc.y + " " + lastmouse.x + " " + lastmouse.y)

        timer = 0

    }

    lastmouse.tick()

}