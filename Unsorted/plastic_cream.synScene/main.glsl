

			//******** Common Code Begins ********


#define R RENDERSIZE.xy

// Dave Hoskins https://www.shadertoy.com/view/4djSRW
vec2 hash23(vec3 p3)
{
	p3 = fract(p3 * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx+33.33);
    return fract((p3.xx+p3.yz)*p3.zy);
}

float gyroid (vec3 seed)
{
    return dot(sin(seed),cos(seed.yzx));
}

float fbm (vec3 seed)
{
    float result = 0.;
    float a = .5;
    for (int i = 0; i < 3; ++i) {
        seed += result / 2.;
        result += gyroid(seed/a)*a;
        a /= 2.;
    }
    return result;
}

			//******** BuffA Code Begins ********


const float speed = 1.;
const float turbulences = 1.;
const float attraction = 2.0;

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    if (FRAMECOUNT < 1)
    {
        fragColor = vec4(hash23(vec3(fragCoord, 0.)), 0, 0);
        return fragColor;
    }

    // coordinates
    vec2 uv = fragCoord/R;
    vec2 mouse = (_mouse.xy - R.xy / 2.)/R.y;
    vec2 p = (fragCoord.xy - R.xy / 2.)/R.y;
    vec2 offset = vec2(0);
    float dist = length(p);
    
    vec4 buffer = texture(BuffA, uv);
    
    // turbulences
    float noise = fbm(vec3(p * 3., dist-smoothTimeC*.1*speed));
    noise = pow(abs(noise), 0.5);
    float angle = noise * 6.28;
    offset += turbulences * vec2(cos(angle), sin(angle));

    // attraction
    offset += (attraction*(1.0+basshits*3.)) * normalize(p) * sin(dist * 9. + smoothTime);
    
    float dt = 60. * 0.016667;
    
    // displace frame buffer
    vec4 frame = texture(BuffA, uv + dt * offset * speed / R);
    
    // edge spawn
    bool spawn = fragCoord.x < 1. || fragCoord.x > R.x-1.
              || fragCoord.y < 1. || fragCoord.y > R.y-1.
              || (_mouse.z > .0 && length(p-mouse) < 50./R.y);
    
    // spawn from noise
    vec2 rng = hash23(vec3(fragCoord, FRAMECOUNT));
    if (spawn) frame = vec4(step(0.5, rng.x),step(0.5, rng.y),0,0);
    
    // neighbor values
    vec2 neighbors = vec2(0);
    for (float x = -1.; x <= 1.; ++x)
    {
        for (float y = -1.; y <= 1.; ++y)
        {
            if (x == 0. && y == 0.) continue;
            neighbors += texture(BuffA, uv+vec2(x,y)/R).rg;
        }
    }
    
    // animation fade
    frame.r += 4.0 * (neighbors.r > 4.0 ? 1. : -1.) * 0.016667;
    frame.g += 4.0 * (neighbors.g > 4.0 ? 1. : -1.) * 0.016667;
    
    fragColor = vec4(clamp(frame.rg, 0., 1.), noise, 1.);//, frame.a + (neighbors.r + neighbors.g) * 0.016667 * .1);
	return fragColor; 
 } 



// Plastic Cream

// i was playing with Conway's game of life cellular automaton
// thinking i had an intuition for a smooth fading version
// with a lifetime gradient and a neighbor ratio

// found by accident that it can make reaction diffusion patterns
// got overwhelmed but couldn't posted another turing pattern feedback
// so i went a bit too far and now it's weird enough as i like it

// iteration from Wasp Blanket https://www.shadertoy.com/view/NlVBz1

#define T(uv) abs(texture(BuffA,uv).b)
#define N(v) normalize(v)

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord/RENDERSIZE.xy;
    vec4 frame = texture(BuffA, uv);
    
    // tints
    vec3 tint = .5 + .5 * cos(vec3(1,2,3)*5. + length(uv-.5)*3.);
    fragColor = vec4(frame.r * tint, 1);
    tint = .5 + .5 * cos(vec3(1,2,3)*5. + length(uv-.5)*2. + 3.);
    fragColor.rgb += frame.g * tint;
    
    // normal
    float height = 1.;
    vec3 unit = vec3(20./RENDERSIZE.xy, 0);
    vec3 normal = normalize(vec3(T(uv+unit.xz)-T(uv-unit.xz),
                                 T(uv-unit.zy)-T(uv+unit.zy),
                                 T(uv) * height));
    
    // light
    fragColor += vec4(.5) * clamp(dot(normal, N(vec3(-1,4,1))), 0., 1.);
    fragColor += vec4(.5) * pow(clamp(dot(normal, N(vec3(-1,1,2))), 0., 1.), 20.);
    
    // shadows
    fragColor -= vec4(0.5) * clamp(dot(normal, N(vec3(0,-4,1))), 0., 1.);
    fragColor *= pow(max(abs(frame.b)-.02,.0), .5);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderMainImage();
	}
}