vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** Common Code Begins ********


#define repeat(p,r) (mod(p,r)-r/2.)
mat2 rot(float a) { return mat2(cos(a),-sin(a),sin(a),cos(a)); }
vec3 lookAt (vec3 from, vec3 at, vec2 uv, float fov)
{
  vec3 z = normalize(at-from);
  vec3 x = normalize(cross(z, vec3(0,1,0)));
  vec3 y = normalize(cross(x, z));
  return normalize(z * fov + uv.x * x + uv.y * y);
}

// Dave Hoskins
// https://www.shadertoy.com/view/4djSRW
vec3 hash33(vec3 p3) {
	p3 = fract(p3 * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yxz+33.33);
    return fract((p3.xxy + p3.yxx)*p3.zyx);
}

			//******** BuffA Code Begins ********


float rng, material;
vec3 coord;

// signed distance function
float map(vec3 p)
{
    float dist = 100.;
    float shape = 200.;
    float t = 206.+smoothTimeC*.1;
    float anim = fract(t);
    float index = floor(t);
    float signal = rng*sin(smoothTimeC*50.)*pow(anim,40.)*.1;
    t = pow(anim,0.1)+index+signal;
    float a = 1.;
    vec3 e = vec3(.01,0.05,0.02);
    const float count = 5.;
    for (float i = 0.; i < count; ++i) {
        
        p.xz = abs(p.xz)-.15*a;
        p.xz *= rot(t/a+i);
        p.yz *= rot(t/a+i);
        p = p - clamp(p, -e*a, e*a);
        dist = min(dist, length(p.xy)-.01*a);
        a /= 1.5;
    }
    
    coord = p;
    
    shape = max((dist - .002), p.z+.005);
    material = shape < dist ? 1. : 0.;
    dist = min(dist, shape);
    
    return dist;
}

void coloring (inout vec3 color, in vec3 pos, in vec3 normal, in vec3 ray, in vec2 uv, in float shade)
{
    // Inigo Quilez color palette
    // https://iquilezles.org/www/articles/palettes/palettes.htm
    vec3 tint = .5+.5*cos(vec3(0,.3,.6)*6.283+smoothTimeB*.012+uv.y*2.);
    vec3 rf = reflect(ray, normal);
    
    if (material == 0.)
    {
        float top = dot(rf, vec3(0,1,0))*.5+.5;
        float glow = dot(normal, ray)*.5+.5;
        color = vec3(0.9) * pow(dot(normal, -normalize(pos))*.5+.5, 0.5);
        color += vec3(.5)*clamp(top,0.,1.);
        color += tint*glow;
        color *= shade;
    }
    else
    {
        float top = dot(rf, vec3(0,1,0))*.5+.5;
        float front = dot(normal, vec3(0,0,1))*.5+.5;
        color = tint;
        color *= smoothstep(.01,.0,sin(coord.z*200.));
        color += vec3(.5)*top;
        color += vec3(.2)*clamp(top,0.,1.);
        color *= shade*front;
    }
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord-RENDERSIZE.xy/2.)/RENDERSIZE.y;
    vec3 color = vec3(.2)*smoothstep(2.,.5,length(uv));
    material = 0.;
    
    // coordinates
    vec3 pos = vec3(0,0,1.2);
    vec3 at = vec3(0);
    pos.xz *= rot(cos(smoothTimeC*.1)*.2);
    pos.zy *= rot(sin(smoothTimeC*.2)*.1);
    vec3 ray = lookAt(pos, at, uv, 1.);
    
    // noise
    vec3 blue = texture(image30, fragCoord/2048.).xyz;
    vec3 white = hash33(vec3(fragCoord, FRAMECOUNT));
    rng = white.x;
    
    // start ahead
    pos += ray * white.z * .2;
    
    // blur edges
    float dof = .0125*smoothstep(.125, 2., length(uv));
    ray.xy += vec2(cos(blue.x*6.28),sin(blue.x*6.28))*white.z*dof*0.01;
    
    // raymarch
    float maxDist = 50.;
    const float count = 100.;
    float steps = 0.;
    float total = 0.;
    for (steps = count; steps > 0.; --steps) {
        float dist = map(pos);
        if (dist < total/RENDERSIZE.y || total > maxDist) break;
        dist *= 0.9+0.1*blue.z;
        ray += white * total*.01;
        pos += ray * dist;
        total += dist;
    }
    
    // coloring
    float shade = steps/count;
    if (shade > .001 && total < maxDist) {
        // NuSan
        // https://www.shadertoy.com/view/3sBGzV
        vec2 noff = vec2(.001,0);
        vec3 normal = normalize(map(pos)-vec3(map(pos-noff.xyy), map(pos-noff.yxy), map(pos-noff.yyx)));
        coloring(color, pos, normal, ray, uv, shade);
    }
    
    fragColor = vec4(color, 1);
	return fragColor; 
 } 


			//******** BuffB Code Begins ********

// Temporal Anti Aliasing from:
// https://www.elopezr.com/temporal-aa-and-the-quest-for-the-holy-trail/

// but only the color clamping...

vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    vec3 color = texture(BuffA, uv).rgb;
    vec3 temporal = texture(BuffB, uv).rgb;
    vec3 minColor = vec3(9999.), maxColor = vec3(-9999.);
    for(int x = -1; x <= 1; ++x){
        for(int y = -1; y <= 1; ++y){
            vec3 c = texture(BuffA, uv + vec2(x, y) / RENDERSIZE.xy).rgb;
            minColor = min(minColor, c);
            maxColor = max(maxColor, c);
        }
    }
    temporal = clamp(temporal, minColor, maxColor);
    fragColor.rgb = mix(color, temporal, 0.9);
    fragColor.a = 1.0;
	return fragColor; 
 } 



// Meta Mecha Tech
// mechanical metamorphoses
// with the color strips a la Chris Foss 

// main code is in Buffer A
// Buffer B is a minimal temporal anti aliasing
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord/RENDERSIZE.xy;
    fragColor = texture(BuffB, uv);
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