

			//******** Common Code Begins ********


// shared constants
const float delay = 15.;
const float max_dist = 200.;

// snippets
#define R RENDERSIZE
#define T(uv) texture(BuffA, uv).a
#define N(x,y,z) normalize(vec3(x,y,z))
#define ss(a,b,t) smoothstep(a,b,t)
#define repeat(p,r) (mod(p,r)-r/2.)
mat2 rot(float a) { return mat2(cos(a),-sin(a),sin(a),cos(a)); }

// Victor Shepardson + Inigo Quilez 
// https://www.shadertoy.com/view/XlXcW4
const uint k = 1103515245U;  // GLIB C
vec3 hash( uvec3 x )
{
    x = ((x>>8U)^x.yzx)*k;
    x = ((x>>8U)^x.yzx)*k;
    x = ((x>>8U)^x.yzx)*k;
    return vec3(x)*(1.0/float(0xffffffffU));
}

// Dave Hoskins
// https://www.shadertoy.com/view/4djSRW
float hash13(vec3 p3)
{
	p3  = fract(p3 * .1031);
    p3 += dot(p3, p3.zyx + 31.32);
    return fract((p3.x + p3.y) * p3.z);
}
vec3 hash31(float p)
{
   vec3 p3 = fract(vec3(p) * vec3(.1031, .1030, .0973));
   p3 += dot(p3, p3.yzx+33.33);
   return fract((p3.xxy+p3.yzz)*p3.zyx); 
}
vec2 hash21(float p)
{
	vec3 p3 = fract(vec3(p) * vec3(.1031, .1030, .0973));
	p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.xx+p3.yz)*p3.zy);
}

			//******** BuffA Code Begins ********


// Frooty by Leon Denise 2023-01-16

// this frame buffer draws the depth buffer with pixel data (steps, timeline, id, depth)
// it draws shape only if it is closer that previous depth
// we can then draw trails of shapes and relight later in Buffer B


// globals
float timeline;
float id, froot;

float map(vec3 p)
{
    vec3 q = p;
    float dist = 100.;
    float shape = 100.;
    
    // grid repeat
    float cell = 5.;
    vec3 pp = p+cell/2.;
    id = hash13(floor(pp/cell));
    p = repeat(pp, cell);
    p += (hash31(floor(id*100.))*2.-1.)*cell/4.;
    
    // timing
    float t = timeline;
    float tt = fract(t);
    float end = mix(.2,.8,id);
    float anim = ss(-.1,end,tt);
    
    // shape parameters
    float r = 1.; // range 
    float f = 1.7; // falloff coeficient
    const float count = 6.; // iterations
    
    // angle
    float n = id*1000.+anim*2. + floor(t); 
    
    // going outter
    r *= pow(ss(0.01,end,tt), .2);
    
    // kaleidoscopic iterated function
    float a = 1.;
    for (float i = 0.; i < count; ++i)
    {
        // twist faster and faster
        p.xy *= rot(n/a);
        p.yz *= rot(n/a);
        
        // fold
        p.x = abs(p.x)-r*a;
        
        // falloff
        a /= f;
    }
    
    // tube (the stems)
    dist = min(dist, max(abs(p.y)-.1*(1.-anim), length(p.xz)-.02));
    
    // sphere (the froots)
    froot = ss(end-.3,end,tt);
    float size = .1 * froot - (1.-froot) * .1;
    dist = min(dist, length(p)-size);
    
    // crop shell to avoid camera collision
    dist = max(dist, -length(q)+cell/2.);
    
    return dist * .8;
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    // 2d coordinates
    vec2 uv = fragCoord/R.xy;
    vec2 p = (fragCoord-R.xy/2.)/R.y;
    
    // noise
    vec3 rng = hash(uvec3(fragCoord, FRAMECOUNT));
    vec3 blu = texture(image14854, fragCoord/1024.).rgb;
    vec2 blur = normalize(blu.xy*2.-1.) * blu.z;
    
    // buffer
    vec4 frame = texture(BuffA, uv);
    
    // timeline used by map for animation
    timeline = rng.x*.001+TIME/delay;

    // 3d coordinates
    vec3 pos = vec3(0,0,0);
    vec3 ray = normalize(vec3(p, -2.));
    vec2 angle = (hash21(floor(TIME/delay))*2.-1.)*.5;
    ray.yz *= rot(angle.x);
    ray.xz *= rot(angle.y);
    
    // raymarch
    float total = 0.;
    float steps = 0.;
    float dof = 0.;
    const float count = 80.;
    for (steps = count; steps > 0.; --steps) {
        float dist = map(pos);
        if (dist < .001 || total > max_dist) break;
        
        // dithering
        dist *= 0.9 + 0.1 * rng.z;
        
        // dof far blur
        dof += 0.0001 * ss(15.,30.,total);
        ray.xy += blur * dof;
        
        total += dist;
        pos += ray * dist;
    }

    // draw if closer
    float depth = frame.a;
    bool closer = total < depth || depth < .001;
    if (total < max_dist && closer)
    {
        // data pack
        fragColor = vec4(
            steps/count, // used to apply shadow
            fract(timeline), // used to offset tint
            floor(id*100.) + froot * .99, // used as material
            total); // depth buffer
    }
    else
    {
        // keep previous result
        fragColor = frame;
    }
    
    // clear between transition
    fragColor *= step(.01, fract(TIME/delay));
	return fragColor; 
 } 


			//******** BuffB Code Begins ********


// Frooty by Leon Denise 2023-01-16

// compute ambient occlusion and normal from Buffer A depth

vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord/RENDERSIZE.xy;
    vec4 color = texture(BuffA, uv);
    
    // ambient occlusion
    float ao = 0.;
    if (color.a < 20. && color.a > 0.5)
    {
        const float count = 6.;
        for (float f = 0.; f < count; ++f)
        {
            vec3 rng = hash(uvec3(fragCoord, float(FRAMECOUNT)+196.*f));
            vec2 offset = 20.*pow(rng.z, 4.)*normalize(rng.xy*2.-1.)/R.xy;
            float total = texture(BuffA, uv+offset).a;
            if (total < 20. && total > 0.5)
                ao += 4.*abs(color.a - total);
        }
    }
    float frame = texture(BuffB, uv).a;
    fragColor.a = mix(frame, 1.-clamp(ao, 0., 1.), .1);
    
    // normal
    vec3 unit = vec3(1./R.xy, 0);
    float w = T(uv+unit.xz);
    float e = T(uv-unit.xz);
    float n = T(uv+unit.zy);
    float s = T(uv-unit.zy);
    bool edge = w * e * n * s < .001;
    if (color.a > .001 && color.a < max_dist && !edge)
    {
        vec3 normal = normalize(vec3(w-e, n-s, color.r*color.r*.05));
        fragColor.rgb = normal;
    }
    else
    {
        fragColor.rgb = vec3(0,0,1);
    }
	return fragColor; 
 } 



// Frooty by Leon Denise 2023-01-16

// a floral iteration over 
// Baroque Fractal Pattern https://www.shadertoy.com/view/flcBD4
// Taste of Noise 7 https://www.shadertoy.com/view/NddSWs

// Buffer A draws shape trails with a depth and id buffer
// Buffer B calculates ambient occlusion and normal
// Image does the coloring and lighting

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec3 color = vec3(0);
    
    // 2d coordinates
    vec2 uv = fragCoord/RENDERSIZE.xy;
    vec2 p = (fragCoord-R.xy/2.)/R.y;
    
    // pixel data from Buffer A (steps, timeline, id, depth)
    vec4 data = texture(BuffA, uv);
    float total = data.a;
    
    if (total > .001 && total < max_dist)
    {
        // data from Buffer A
        float shade = data.r;
        float time = data.g;
        float mat = floor(data.b);
        float froot = fract(data.b);
        
        // data from Buffer B (normal, ao)
        data = texture(BuffB, uv);
        vec3 normal = data.xyz;
        float ao = data.a;
        
        // palettes
        // Inigo Quilez https://iquilezles.org/articles/palettes
        vec3 tintCold = .5+.5*cos(vec3(2,3,1)*5.7+4.-time*4.);
        vec3 tintWarm = .5+.5*cos(vec3(1,2,3)*4.9+mat*1.+uv.y*2.+5.);
        color = mix(tintCold, tintWarm, froot);
        
        // lighting
        color += pow(dot(normal, N(0,1,1))*.5+.5, 4.)*froot;
        //color += .2*(1.-pow(abs(dot(normal, N(0,0,1))), .5));
        
        // shadow
        color *= shade;
        color *= ao * .5 + .5;
        color *= ss(100.,0.,total);
    }
    else
    {
        // background
        color = vec3(0.039,0.110,0.322) * smoothstep(2., -2., length(uv-.5));
    }
    
    // fade transition
    float t = fract(TIME/delay);
    color *= ss(.0,.1,t) * ss(1.,.9,t);
    
    // vignette
    color *= ss(1.5, .5, length(p));
    
    fragColor = vec4(color, 1);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderPassB();
	}
	if(PASSINDEX == 2){
		return renderMainImage();
	}
}