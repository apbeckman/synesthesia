function CirclePos() {
    this.x = 0
    this.y = 0
    this.goal_x = 0
    this.goal_y = 0
}

// CirclePos.prototype.tick = function(x, y) {
//     this.x += (this.goal_x - this.x) * 0.4
//     this.y += (this.goal_y - this.y) * 0.4
// }

// CirclePos.prototype.setGoal = function(x, y) {
//     this.goal_x = x
//     this.goal_y = y
// }

CirclePos.prototype.move = function(x, y) {
    this.x = x
    this.y = y
}

var circle1 = new CirclePos()
uniforms.circle1 = circle1

function setup() {
    // onChange etc here
}

function update() {
    // circle1.tick()
    if (_click.x > 0.) {
        circle1.move(_muvc.x, _muvc.y)
        // circle1.setGoal(_muvc.x, _muvc.y)
    }
}