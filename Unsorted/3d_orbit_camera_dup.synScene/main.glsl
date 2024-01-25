/* SHAPE FUNCTIONS:
* 
*  These generate Signed Distance Functions, or SDFs, of basic 3 dimensional shapes.
*  You can combine these together in map() to compose a scene using different boolean operations.
*  You can learn more about SDFs and find examples at this link:
*  https://iquilezles.org/articles/distfunctions/
*/
float sdBox(vec3 p, vec3 s) {
    p = abs(p)-s;
	return length(max(p, 0.))+min(max(p.x, max(p.y, p.z)), 0.);
}

float sdSphere(vec3 p, float s) {
    return length(p) - s;
}

/* MAP:
*  
*  This is where we combine our SDFs to create a 3D scene.
*  map() receives a vec3 position in space and returns a float distance to the nearest surface.
*  We break this out into its own function so we can call it in the raymarcher,
*  the getNorm() function, and anywhere else we want to know about our scene.
*
*  You can use this scene as an SDF editing playground by editing the map function.
*  Use the controls to view your work as you go, and then either copy the contents
*  to another scene or edit this scene to be more stylized.
*/
float map(in vec3 p) {
    float d = 0.;
    if (shape < 0.5) { 
        // shape 0: Sphere
        float sphere1 = sdSphere(p,0.25); // a sphere of size 0.25 at the origin
        d = sphere1; // add a sphere to the scene. since it's the first object we can just set d equal to our sdf
    } else if (shape < 1.5) { 
        // shape 1: Box Minus Sphere
        float sphere1 = sdSphere(p,shape_mod/8.+.2); // set sphere size based on shape_mod user control
        d = sphere1; // add sphere to the scene
        float box1 = sdBox(p,vec3(0.2)); // a box of size 0.2 at the origin
        d = max(-d,box1); // boolean subtract the existing scene (the sphere) and the box
    } else if (shape < 2.5) {
        // shape 2: Gyroid
        float sphere1 = sdSphere(p,0.25); // start with a sphere, this will be our outer bounds
        d = sphere1;
        // create a gyroid
        float scale = shape_mod*29. + .01;
        float gyroid = abs(dot(sin(p*scale), cos(p.zxy*scale))) / scale - .05;
        // 3D boolean AND operation with sphere
        d = max(d, -gyroid);
    }
    
    return d;
}

/* GET NORMAL:
*
*  getNorm accepts a 3d point on a surface in space (returned from the raymarch function),
*  samples 3 nearby points, and then returns a 3d vector perpendicular to the surface
*  Useful in shading and lighting. You shouldn't need to edit this function.
*/
vec3 getNorm(vec3 p) {
	float d = map(p);
    vec2 e = vec2(.001, 0);
    vec3 n = d - vec3(
        map(p-e.xyy),
        map(p-e.yxy),
        map(p-e.yyx));
    return normalize(n);
}

/* RAYMARCH:
*
*  This function starts at a position ro and direction rd in scene space,
*  and runs a loop where it calls the map function to get the nearest surface in the scene,
*  then moves forward that distance (in the rd direction), until it hits a surface or runs 
*  out of steps, and gives us the final position in space and the depth of the ray.
*  You can modify the for loop to optimize your scene, create optical effects, add wild glitches,
*  and more, but in most cases you can leave that loop alone.
*  
*  We're also using the distance, surface normals, and scene position to calculate
*  shading and lighting, and returning the color results as our final output.
*  We've provided a few basic views but you'll probably want to edit the 
*  colorizing section to create your own custom stylized look.
*/

const int max_steps = 100;

vec4 raymarch(in vec3 ro, in vec3 rd) {
    // initialize variables
    vec3 col = vec3(0.);
    vec3 p = ro;
    float d = 0.;
    int steps = 0;
    bool hit = false;
    // iteratively step forward, looking for a surface
    for (steps; steps<max_steps; steps++) {
        p = ro + rd * d;
        float dist_step = map(p);
        d += dist_step;
        if (dist_step < 0.001) {
            hit = true;
            break;
        }
    }
    
    // Colorizing: now that we have a position in space, we can use that to calculate
    // lighting and shading, and then use that to colorize our scene.    
    if (hit) { // color the objects
        // the direction perpendicular to the surface at the point p
        vec3 norm = getNorm(p);
        if (shader_mode < 0.5) { 
            // shader 0: normals
            col = abs(norm);
        } else if (shader_mode < 1.5) { 
            // shader 1: depth
            col = 1. - vec3(d*d*0.5,d*d,d*.5);
        } else if (shader_mode < 2.5) { 
            // shader 2: world space grid

            // evenly spaced color grids
            col = smoothstep(0.9,0.91,fract(abs(p)*10.));
            // white center lines
            col += smoothstep(0.01,0.0,abs(p.x));
            col += smoothstep(0.01,0.0,abs(p.y));
            col += smoothstep(0.01,0.0,abs(p.z));
        } else if (shader_mode < 3.5) { 
            // shader 3: basic lighting
            vec3 objcol = vec3(1.); // white surface
            vec3 lightcol = vec3(0.9,0.85,0.6); // warm light
            vec3 lightpos = vec3(2.,2.,-0.5); // place light at top right of default view, slightly in front of object
            vec3 lightdir = normalize(lightpos - p); // calculate direction to the light
            vec3 diff = max(dot(norm,lightdir),0.) * lightcol;
            vec3 amb = lightcol * .2; // ambient light, multiplier determines strength
            col = (amb + diff) * objcol;
            /* Learn more about lighting here: 
            *    https://learnopengl.com/Lighting/Basic-Lighting
            *  Or if you prefer video tutorials: 
            *    https://www.youtube.com/watch?v=mL8U8tIiRRg 
            *    "Healthbars, SDFs & Lighting • Shaders for Game Devs [Part 2]" by Freya Holmér
            */
        }
    } else { // didn't hit object, set color of the background
        col = vec3(0.03,0.02,0.01);
    }
    
    return vec4(col, 1.);
}

/* MAIN VIEW:
*
*  This function receives all of the camera position data from the UI and calls
*  the raymarch function from that camera viewpoint. This is used to create the
*  primary view, as well as the top right view in grid view.
*/
vec4 mainView(in vec2 uv) {
    
    // rotate & move cam based on controls
    vec3 ro = vec3(0.,0.,cam_zoom);
    vec3 rd = normalize(vec3(uv, 1.0));
    
    ro.yz = _rotate(ro.yz, cam_xy.y);
    ro.xz = _rotate(ro.xz, cam_xy.x);
    
    rd.yz = _rotate(rd.yz, cam_xy.y);
    rd.xz = _rotate(rd.xz, cam_xy.x);
    
    return raymarch(ro, rd);
}

/* STATIC VIEW
*
*  This function is called to raymarch the scene from 3 static angles in grid view.
*/
vec4 staticView(in vec3 viewdir, in vec2 uv) {
    // set camera position
    vec3 ro = vec3(
        viewdir.x,
        viewdir.y,
        viewdir.z
    );
    uv *= -cam_zoom;
    // orthogonal camera - spaced out ray origins, parallel rays
    if (viewdir.y < -0.5) { // top left
        ro.x += uv.x;
        ro.z -= uv.y;
    }
    if (viewdir.z < -0.5) { // bottom left
        ro.x += uv.x;
        ro.y += uv.y;
    }
    if (viewdir.x < -0.5) { // bottom right
        ro.z -= uv.x;
        ro.y += uv.y;
    }
    vec3 rd = normalize(vec3(0.,0.,1.));
    // rotate to match viewdir
    rd.xz = _rotate(rd.xz, -viewdir.x * PI / 2.);
    rd.yz = _rotate(rd.yz, -viewdir.y * PI / 2.);
    
    return raymarch(ro, rd);
}

/* RENDER MAIN
*
*  This is the main render loop that Synesthesia looks for, whatever you return 
*  here is shown on the screen. Here we're using it as a wrapper to implement "grid_view".
*  If grid_view is on we divide the screen into quadrants to view 4 different angles,
*  otherwise we just call mainView for the full screen. 
*/
vec4 renderMain(void) {
    
    vec2 aspect = RENDERSIZE.xy/RENDERSIZE.x;
    
    if (grid_view < 0.5) return mainView(_uvc);
    
    vec2 uvc2 = _uvc * 2.;
    
    if (_uv.x<0.5) // left
        if (_uv.y>0.5) // top left - top view
            return staticView(vec3(0,-1.,0),uvc2 + aspect * vec2(1.,-1.));
        else // bottom left - front view
            return staticView(vec3(0,0,-1.),uvc2 + aspect);
    else // right
        if (_uv.y>0.5) // top right - variable
            return mainView(uvc2 - aspect);
        else // bottom right - right view
            return staticView(vec3(-1.,0,0),uvc2 + aspect * vec2(-1.,1.));
}