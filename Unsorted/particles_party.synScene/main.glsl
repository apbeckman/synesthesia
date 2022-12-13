

			//******** Common Code Begins ********


float gyroid (vec3 seed) { return dot(sin(seed),cos(seed.yzx)); }
float fbm (vec3 seed)
{
    float result = 0., a = .5;
    for (int i = 0; i < 3; ++i, a/=2.)
    {
        result += (gyroid(seed/a))*a;
    }
    return result;
}

// Dave Hoskins https://www.shadertoy.com/view/4djSRW
float hash11(float p)
{
    p = fract(p * .1031);
    p *= p + 33.33;
    p *= p + p;
    return fract(p);
}
vec2 hash21(float p)
{
	vec3 p3 = fract(vec3(p) * vec3(.1031, .1030, .0973));
	p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.xx+p3.yz)*p3.zy);

}
vec4 hash41(float p)
{
	vec4 p4 = fract(vec4(p) * vec4(.1031, .1030, .0973, .1099));
    p4 += dot(p4, p4.wzxy+33.33);
    return fract((p4.xxyz+p4.yzzw)*p4.zywx);
}

			//******** BuffA Code Begins ********


const float count = 100.;

float speed = .3;
float friction = 3.;
float fade = 0.1;
float thin = 0.02;

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    // coordinates
    vec2 uv = fragCoord/RENDERSIZE.xy;
    vec2 p = (2.*fragCoord-RENDERSIZE.xy)/RENDERSIZE.y;
    
    // buffer
    vec4 frame = texture(BuffA, uv);
    
    // pixels are data (bottom left)
    if (fragCoord.y < 1. && fragCoord.x < count)
    {
        float id = fragCoord.x;
            
        // init
        if (FRAMECOUNT <= 1)
        {
            // random position and velocity
            frame = vec4(hash41(id)*2.-1.);
            frame.zw *= .01;
        }
        else
        {
            // coordinates
            vec2 aspect = vec2(RENDERSIZE.x/RENDERSIZE.y, 1.);
            vec2 p = frame.xy;
            vec2 offset = vec2(0);
            vec2 target = vec2(0);
            
            // respawn
            float t = smoothTime * 10.;
            float idd = id+floor(t) * count;
            if (hash11(idd) > .95 && fract(t) < .1)
            {
                frame = hash41(idd)*2.-1.;
                frame.xy *= aspect;
                frame.zw *= .01;
                fragColor = frame;
                return fragColor;
            }
            
            // interaction
            if (_mouse.z > 0.)
            {
                target = (2.*_mouse.xy-RENDERSIZE.xy)/RENDERSIZE.y;
            }
            
            // curl
            float noise = fbm(vec3(p, length(p) + smoothTimeC));
            float a = noise * 6.28;
            offset += vec2(cos(a), sin(a));
            
            // target
            offset += normalize(target.xy-p) * 2. * length(target.xy-p);
            
            // jitter
            offset += (hash21(id)*2.-1.)*(.5+.5*sin(smoothTimeC));
            
            // inertia
            vec2 velocity = frame.zw;
            velocity = velocity * (1.-friction*0.016667) + offset * speed * 0.016667;
            
            // apply
            frame.xy += velocity;
            frame.zw = velocity;
        }
    }
    
    // pixels are colors
    else
    {
        float matID = 0.;
        float dist = 100.;
        float dither = texture(image14854, fragCoord/1920.).r;

        for (float i = 0.; i < count; ++i)
        {
            // iterate pixel data
            vec4 data = texelFetch(BuffA, ivec2(i,0), 0);
            
            // circle shape (jitter blending with previous pos)
            vec2 pos = data.xy - data.zw * dither;
            float shape = length(pos-p);
            matID = shape < dist ? i : matID;
            dist = min(dist, shape);
        }

        // grayscale
        float shade = smoothstep(thin,.0,dist);

        // buffer
        frame.r = max(frame.r - fade, shade);
        
        // material layer
        if (dist < thin) frame.g = matID;
    }
    
    fragColor = frame;
	return fragColor; 
 } 



// Particles Party

// simple colorful particles system
// will try voronoi tracking next time

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec3 color = vec3(0);
    
    // coordinates
    vec2 uv = fragCoord/RENDERSIZE.xy;
    float rng = texture(image14854, fragCoord/1024.).r;
    vec2 aspect = vec2(RENDERSIZE.y/RENDERSIZE.x, 1.);
    
    // data
    vec4 data = texture(BuffA, uv);
    float shade = data.r;
    float mat = data.g;
    
    // rainbow
    color = .5+.5*cos(vec3(1,2,3)*4.9 + mat);
    
    // light
    vec3 un = vec3(0.005*aspect, 0);
    #define T(un) texture(BuffA, uv+un).r
    vec3 normal = normalize(vec3(T(un.xz)-T(-un.xz),T(un.zy)-T(-un.zy), .5));
    float d = dot(normal, normalize(vec3(0,-2,1)))*.5+.5;
    color += pow(d, 10.);
    
    // shadow
    color *= smoothstep(.0,.01,shade);

    fragColor = vec4(color, 1);
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