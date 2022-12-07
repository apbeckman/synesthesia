

			//******** BuffA Code Begins ********

// main reaction-diffusion loop

// actually the diffusion is realized as a separated two-pass Gaussian blur kernel and is stored in buffer C

#define pi2_inv 0.159154943091895335768883763372

vec2 complex_mul(vec2 factorA, vec2 factorB){
    return vec2( factorA.x*factorB.x - factorA.y*factorB.y, factorA.x*factorB.y + factorA.y*factorB.x);
}

vec2 spiralzoom(vec2 domain, vec2 center, float n, float spiral_factor, float zoom_factor, vec2 pos){
    vec2 uv = domain - center;
    float d = length(uv);
    return vec2( atan(uv.y, uv.x)*n*pi2_inv + d*spiral_factor, -log(d)*zoom_factor) + pos;
}

vec2 complex_div(vec2 numerator, vec2 denominator){
    return vec2( numerator.x*denominator.x + numerator.y*denominator.y,
                numerator.y*denominator.x - numerator.x*denominator.y)/
        vec2(denominator.x*denominator.x + denominator.y*denominator.y);
}

float circle(vec2 uv, vec2 aspect, float scale){
    return clamp( 1. - length((uv-0.5)*aspect*scale), 0., 1.);
}

float sigmoid(float x) {
    return 2./(1. + exp2(-x)) - 1.;
}

float smoothcircle(vec2 uv, vec2 aspect, float radius, float ramp){
    return 0.5 - sigmoid( ( length( (uv - 0.5) * aspect) - radius) * ramp) * 0.5;
}

float conetip(vec2 uv, vec2 pos, float size, float min)
{
    vec2 aspect = vec2(1.,RENDERSIZE.y/RENDERSIZE.x);
    return max( min, 1. - length((uv - pos) * aspect / size) );
}

float warpFilter(vec2 uv, vec2 pos, float size, float ramp)
{
    return 0.5 + sigmoid( conetip(uv, pos, size, -16.) * ramp) * 0.5;
}

vec2 vortex_warp(vec2 uv, vec2 pos, float size, float ramp, vec2 rot)
{
    vec2 aspect = vec2(1.,RENDERSIZE.y/RENDERSIZE.x);

    vec2 pos_correct = 0.5 + (pos - 0.5);
    vec2 rot_uv = pos_correct + complex_mul((uv - pos_correct)*aspect, rot)/aspect;
    float _filter = warpFilter(uv, pos_correct, size, ramp);
    return mix(uv, rot_uv, _filter);
}

vec2 vortex_pair_warp(vec2 uv, vec2 pos, vec2 vel)
{
    vec2 aspect = vec2(1.,RENDERSIZE.y/RENDERSIZE.x);
    float ramp = 5.;

    float d = 0.2;

    float l = length(vel);
    vec2 p1 = pos;
    vec2 p2 = pos;

    if(l > 0.){
        vec2 normal = normalize(vel.yx * vec2(-1., 1.))/aspect;
        p1 = pos - normal * d / 2.;
        p2 = pos + normal * d / 2.;
    }

    float w = l / d * 2.;

    // two overlapping rotations that would annihilate when they were not displaced.
    vec2 circle1 = vortex_warp(uv, p1, d, ramp, vec2(cos(w),sin(w)));
    vec2 circle2 = vortex_warp(uv, p2, d, ramp, vec2(cos(-w),sin(-w)));
    return (circle1 + circle2) / 2.;
}

vec2 mouseDelta(){
    vec2 pixelSize = 1. / RENDERSIZE.xy;
    float eighth = 1./8.;
    vec4 oldMouse = texture(BuffD, vec2(7.5 * eighth, 2.5 * eighth));
    vec4 nowMouse = vec4(_mouse.xy / RENDERSIZE.xy, _mouse.zw / RENDERSIZE.xy);
    if(oldMouse.z > pixelSize.x && oldMouse.w > pixelSize.y && 
       nowMouse.z > pixelSize.x && nowMouse.w > pixelSize.y)
    {
        return nowMouse.xy - oldMouse.xy;
    }
    return vec2(0.);
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    vec2 pixelSize = 1. / RENDERSIZE.xy;
    

    vec2 mouseV = mouseDelta();
    vec2 aspect = vec2(1.,RENDERSIZE.y/RENDERSIZE.x);
    uv = vortex_pair_warp(uv, _mouse.xy*pixelSize, mouseV*aspect*1.4);

    vec4 blur1 = texture(BuffC, uv);
    
    vec4 noise = texture(image30, fragCoord.xy / RENDERSIZE.xy + fract(vec2(42,56)*TIME));

    // get the gradients from the blurred image
	vec2 d = pixelSize*4.;
	vec4 dx = (texture(BuffC, fract(uv + vec2(1,0)*d)) - texture(BuffC, fract(uv - vec2(1,0)*d))) * 0.5;
	vec4 dy = (texture(BuffC, fract(uv + vec2(0,1)*d)) - texture(BuffC, fract(uv - vec2(0,1)*d))) * 0.5;
    
    vec2 uv_red = uv + vec2(dx.x, dy.x)*pixelSize*8.; // add some diffusive expansion
    
    float new_red = texture(BuffA, fract(uv_red)).x + (noise.x - 0.5) * 0.0025 - 0.002; // stochastic decay
	new_red -= (texture(BuffC, fract(uv_red + (noise.xy-0.5)*pixelSize)).x -
				texture(BuffA, fract(uv_red + (noise.xy-0.5)*pixelSize))).x * 0.047; // reaction-diffusion
        
    if(FRAMECOUNT<10)
    {
        fragColor = noise; 
    }
    else
    {
        fragColor.x = clamp(new_red, 0., 1.);
    }

//    fragColor = noise; // need a restart?
	return fragColor; 
 } 


			//******** BuffB Code Begins ********

// horizontal Gaussian blur pass

vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 pixelSize = 1./ RENDERSIZE.xy;
    vec2 uv = fragCoord.xy * pixelSize;
    
    float h = pixelSize.x;
	vec4 sum = vec4(0.0);
	sum += texture(BuffA, fract(vec2(uv.x - 4.0*h, uv.y)) ) * 0.05;
	sum += texture(BuffA, fract(vec2(uv.x - 3.0*h, uv.y)) ) * 0.09;
	sum += texture(BuffA, fract(vec2(uv.x - 2.0*h, uv.y)) ) * 0.12;
	sum += texture(BuffA, fract(vec2(uv.x - 1.0*h, uv.y)) ) * 0.15;
	sum += texture(BuffA, fract(vec2(uv.x + 0.0*h, uv.y)) ) * 0.16;
	sum += texture(BuffA, fract(vec2(uv.x + 1.0*h, uv.y)) ) * 0.15;
	sum += texture(BuffA, fract(vec2(uv.x + 2.0*h, uv.y)) ) * 0.12;
	sum += texture(BuffA, fract(vec2(uv.x + 3.0*h, uv.y)) ) * 0.09;
	sum += texture(BuffA, fract(vec2(uv.x + 4.0*h, uv.y)) ) * 0.05;
    
    fragColor.xyz = sum.xyz/0.98; // normalize
	fragColor.a = 1.;
	return fragColor; 
 } 


			//******** BuffC Code Begins ********

// vertical Gaussian blur pass

vec4 renderPassC() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 pixelSize = 1./ RENDERSIZE.xy;
    vec2 uv = fragCoord.xy * pixelSize;

    float v = pixelSize.y;
	vec4 sum = vec4(0.0);
	sum += texture(BuffB, fract(vec2(uv.x, uv.y - 4.0*v)) ) * 0.05;
	sum += texture(BuffB, fract(vec2(uv.x, uv.y - 3.0*v)) ) * 0.09;
	sum += texture(BuffB, fract(vec2(uv.x, uv.y - 2.0*v)) ) * 0.12;
	sum += texture(BuffB, fract(vec2(uv.x, uv.y - 1.0*v)) ) * 0.15;
	sum += texture(BuffB, fract(vec2(uv.x, uv.y + 0.0*v)) ) * 0.16;
	sum += texture(BuffB, fract(vec2(uv.x, uv.y + 1.0*v)) ) * 0.15;
	sum += texture(BuffB, fract(vec2(uv.x, uv.y + 2.0*v)) ) * 0.12;
	sum += texture(BuffB, fract(vec2(uv.x, uv.y + 3.0*v)) ) * 0.09;
	sum += texture(BuffB, fract(vec2(uv.x, uv.y + 4.0*v)) ) * 0.05;
    
    fragColor.xyz = sum.xyz/0.98; // normalize
	fragColor.a = 1.;
	return fragColor; 
 } 


			//******** BuffD Code Begins ********

// not used (yet), but hooray for 8 channel feedback

vec4 renderPassD() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    vec2 pixelSize = 1. / RENDERSIZE.xy;
    float eighth = 1./8.;
    if(uv.x > 7.*eighth && uv.x < 8.*eighth && uv.y > 2.*eighth && uv.y < 3.*eighth)
    {
        fragColor = vec4(_mouse.xy / RENDERSIZE.xy, _mouse.zw / RENDERSIZE.xy);
    }
	return fragColor; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    vec2 pixelSize = 1. / RENDERSIZE.xy;
    vec2 aspect = vec2(1.,RENDERSIZE.y/RENDERSIZE.x);

    vec4 noise = texture(image30, fragCoord.xy / RENDERSIZE.xy + fract(vec2(42,56)*TIME));
    
	vec2 lightSize=vec2(4.);

    // get the gradients from the blurred image
	vec2 d = pixelSize*2.;
	vec4 dx = (texture(BuffD, uv + vec2(1,0)*d) - texture(BuffD, uv - vec2(1,0)*d))*0.5;
	vec4 dy = (texture(BuffD, uv + vec2(0,1)*d) - texture(BuffD, uv - vec2(0,1)*d))*0.5;

	// add the pixel gradients
	d = pixelSize*1.;
	dx += texture(BuffA, uv + vec2(1,0)*d) - texture(BuffA, uv - vec2(1,0)*d);
	dy += texture(BuffA, uv + vec2(0,1)*d) - texture(BuffA, uv - vec2(0,1)*d);
    vec2 circ= vec2(cos(smoothTimeB), sin(smoothTimeB));
	vec2 displacement = vec2(dx.x,dy.x)*lightSize; // using only the red gradient as displacement vector
	float light = pow(max(1.-distance(0.5+(uv-0.5)*aspect*lightSize + displacement,0.5+(_mouse.xy*pixelSize-0.5)*aspect*circ+lightSize),0.),4.);

	// recolor the red channel
	vec4 rd = vec4(texture(BuffA,uv+vec2(dx.x,dy.x)*pixelSize*8.).x)*vec4(0.7,1.5,2.0,1.0)-vec4(0.3,1.0,1.0,1.0);

    // and add the light map
    fragColor = mix(rd,vec4(8.0,6.,2.,1.), light*0.75*vec4(1.-texture(BuffA,uv+vec2(dx.x,dy.x)*pixelSize*8.).x)); 
	
	//fragColor = texture(BuffA, uv); // bypass    
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
		return renderPassC();
	}
	if(PASSINDEX == 3){
		return renderPassD();
	}
	if(PASSINDEX == 4){
		return renderMainImage();
	}
}