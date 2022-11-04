

			//******** Common Code Begins ********

vec2 R;
vec4 M;
int I;
#define N 12.
#define R3D vec3(R/N,N*N)
#define e(a) mat2(cos(a),-sin(a),sin(a),cos(a))
#define d2(U) ((U).xy+vec2(mod(floor((U).z),N),floor(floor((U).z)/N))*R/N)
#define d3(u) vec3(mod(u,R/N),floor(u/R*N).x+floor(u/R*N).y*N)
#define _3D  vec3 U = d3(u)
#define Sampler vec4 T(vec3 U) {return mix(texture(iChannel0,d2(vec3(U.xy,floor(U.z)))/R),texture(iChannel0,d2(vec3(U.xy, ceil(U.z)))/R),fract(U.z));}
#define A(U) texture(cha,d2(mod(U,R3D))/R)
#define B(U) texture(chb,d2(mod(U,R3D))/R)
#define Main void mainImage (out vec4 Q, in vec2 U)
float signe (float x) {return atan(100.*x);}
void prog (vec3 U, out vec4 a, out vec4 b, sampler2D cha, sampler2D chb) {
	
    a = vec4(0); b = vec4(0);
    for (int x = -1; x <= 1; x++)
    for (int y = -1; y <= 1; y++)
    for (int z = -1; z <= 1; z++)
    {
        vec3 u = vec3(x,y,z);
    	vec4 aa = A(U+u), bb = B(U+u);
        aa.xyz += bb.xyz;
        #define q 1.1
		vec3 w1 = clamp(aa.xyz-0.5*q,U - 0.5,U + 0.5),
             w2 = clamp(aa.xyz+0.5*q,U - 0.5,U + 0.5);
        float m = (w2.x-w1.x)*(w2.y-w1.y)*(w2.z-w1.z)/(q*q*q);
        aa.xyz = 0.5*(w1+w2);
        a.xyz += aa.xyz*aa.w*m;
        b.xyz += bb.xyz*aa.w*m;
        a.w += aa.w*m;
    }
    if (a.w>0.) {
        a.xyz/=a.w;
        b.xyz/=a.w;
    }
}
void prog2 (vec3 U, out vec4 a, out vec4 b, sampler2D cha, sampler2D chb) {
	
    a = A(U); b = B(U);
    vec3 f = vec3(0); float m = 0.;
    for (int x = -1; x <= 1; x++)
    for (int y = -1; y <= 1; y++)
    for (int z = -1; z <= 1; z++)
    {
        vec3 u = vec3(x,y,z);
        float l = length(u);
        if (l>0.) {
    		vec4 aa = A(U+u), bb = B(U+u);
            f += 1e-2*(aa.w*(1.-.2*aa.w))*u/l;
            m += aa.w;
        }
    }
    if (m>0.) b.xyz += f/m;
    
    
    // Boundaries:
   	b.xyz -= 1e-3*signe(a.w)*(a.xyz-0.5*R3D)*sin(1e-5*float(I));

    
    if (I<1||U.x<1.||R3D.x-U.x<1.||R3D.y-U.y<1.||R3D.x-U.x<1.||U.z<1.||R3D.z-U.z<1.) {
    	a = vec4(U,0);
        b = vec4(0);
        if (length(U-0.5*R3D) < 0.2*R3D.y) a.w = 20.;
    }
}

			//******** BuffA Code Begins ********

#define pi2_inv 0.159154943091895335768883763372

vec2 lower_left(vec2 uv)
{
    return fract(uv * 0.5);
}

vec2 lower_right(vec2 uv)
{
    return fract((uv - vec2(1, 0.)) * 0.5);
}

vec2 upper_left(vec2 uv)
{
    return fract((uv - vec2(0., 1)) * 0.5);
}

vec2 upper_right(vec2 uv)
{
    return fract((uv - 1.) * 0.5);
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

vec4 BlurA(vec2 uv, int level)
{
    if(level <= 0)
    {
        return texture(BuffA, fract(uv));
    }

    uv = upper_left(uv);
    for(int depth = 1; depth < 8; depth++)
    {
        if(depth >= level)
        {
            break;
        }
        uv = lower_right(uv);
    }

    return texture(BuffD, uv);
}
	
vec2 GradientA(vec2 uv, vec2 d, vec4 selector, int level){
	vec4 dX = 0.5*BlurA(uv + vec2(1.,0.)*d, level) - 0.5*BlurA(uv - vec2(1.,0.)*d, level)+highhits*0.0125;
	vec4 dY = 0.5*BlurA(uv + vec2(0.,1.)*d, level) - 0.5*BlurA(uv - vec2(0.,1.)*d, level);
	return vec2( dot(dX, selector), dot(dY, selector) );
}

vec2 rot90(vec2 vector){
	return vector.yx*vec2(1,-1);
}

vec2 complex_mul(vec2 factorA, vec2 factorB){
    return vec2( factorA.x*factorB.x - factorA.y*factorB.y, factorA.x*factorB.y + factorA.y*factorB.x);//add control for succ
}

vec2 spiralzoom(vec2 domain, vec2 center, float n, float spiral_factor, float zoom_factor, vec2 pos){
    vec2 uv = domain - center;
    float d = length(uv);
    return vec2( atan(uv.y, uv.x)*n*pi2_inv + d*(spiral_factor*(0.75+basshits)), -log(d)*(zoom_factor*(0.75+syn_MidLevel))) + pos;
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
    float filterv = warpFilter(uv, pos_correct, size, ramp);
    return mix(uv, rot_uv, filterv);
}

vec2 vortex_pair_warp(vec2 uv, vec2 pos, vec2 vel)
{
    vec2 aspect = vec2(1.,RENDERSIZE.y/RENDERSIZE.x);
    float ramp = 4.;

    float d = 0.125;

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

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    vec4 noise = (texture(image30, fragCoord.xy / RENDERSIZE.xy + fract(vec2(42,56)*TIME))-0.5)*2.;

    if(FRAMECOUNT<10)
    {
        fragColor = noise;
        return fragColor;
    }

    vec2 mouseV = mouseDelta();
    if(length(mouseV)==0.){
        fragColor = BlurA(uv, 0);
        //return;
    }
    vec2 aspect = vec2(1.,RENDERSIZE.y/RENDERSIZE.x);
    vec2 pixelSize = 1. / RENDERSIZE.xy;

    uv = vortex_pair_warp(uv, _mouse.xy*pixelSize, mouseV*aspect*1.4);
    
    // expansion
    vec2 gradientLookupDistance = pixelSize*3.;
    float expansionFactor = 1.;
    
    // reaction-diffusion  
    float differentialFactor = 12./256.;
    float increment = - 3./256.;
    float noiseFactor = 2./256.;
    
    // rock-paper-scissor
    float feedBack = 6./(256.*(1.25+basshits));
    float feedForward = 6./(256.*(1.25-basshits));

	fragColor.r = BlurA(uv + GradientA(uv, gradientLookupDistance, vec4(4.,0.,-2.,0.), 1)*pixelSize*expansionFactor, 0).r;
	fragColor.g = BlurA(uv + GradientA(uv, gradientLookupDistance, vec4(0.,4.,0.,-2.), 1)*pixelSize*expansionFactor, 0).g;
	fragColor.b = BlurA(uv + GradientA(uv, gradientLookupDistance, vec4(-2.,0.,4.,0.), 1)*pixelSize*expansionFactor, 0).b;
    fragColor.a = BlurA(uv + GradientA(uv, gradientLookupDistance, vec4(0.,-2.,0.,4.), 1)*pixelSize*expansionFactor, 0).a;

   	fragColor += (BlurA(uv, 1) - BlurA(uv, 2))*(differentialFactor-basshits/100);

    fragColor += increment + noise * noiseFactor;

    fragColor -= fragColor.argb * feedBack;
    fragColor += fragColor.gbar * feedForward;
    
    fragColor = clamp(fragColor, 0., 1.);

//    fragColor = noise; // reset
	return fragColor; 
 } 


			//******** BuffB Code Begins ********

#define pi2_inv 0.159154943091895335768883763372
/*
vec2 lower_left(vec2 uv)
{
    return fract(uv * 0.5);
}

vec2 lower_right(vec2 uv)
{
    return fract((uv - vec2(1, 0.)) * 0.5);
}

vec2 upper_left(vec2 uv)
{
    return fract((uv - vec2(0., 1)) * 0.5);
}

vec2 upper_right(vec2 uv)
{
    return fract((uv - 1.) * 0.5);
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

vec4 BlurA(vec2 uv, int level)
{
    if(level <= 0)
    {
        return texture(iChannel0, fract(uv));
    }

    uv = upper_left(uv);
    for(int depth = 1; depth < 8; depth++)
    {
        if(depth >= level)
        {
            break;
        }
        uv = lower_right(uv);
    }

    return texture(BuffD, uv);
}
*/
vec4 BlurB(vec2 uv, int level)
{
    if(level <= 0)
    {
        return texture(BuffB, fract(uv));
    }

    uv = lower_left(uv);
    for(int depth = 1; depth < 8; depth++)
    {
        if(depth >= level)
        {
            break;
        }
        uv = lower_right(uv);
    }

    return texture(BuffD, uv);
}
/*
vec2 GradientA(vec2 uv, vec2 d, vec4 selector, int level){
    vec4 dX = 0.5*BlurA(uv + vec2(1.,0.)*d, level) - 0.5*BlurA(uv - vec2(1.,0.)*d, level);
    vec4 dY = 0.5*BlurA(uv + vec2(0.,1.)*d, level) - 0.5*BlurA(uv - vec2(0.,1.)*d, level);
    return vec2( dot(dX, selector), dot(dY, selector) );
}

vec2 rot90(vec2 vector){
    return vector.yx*vec2(1,-1);
}

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
    float filterv = warpFilter(uv, pos_correct, size, ramp);
    return mix(uv, rot_uv, filterv);
}

vec2 vortex_pair_warp(vec2 uv, vec2 pos, vec2 vel)
{
    vec2 aspect = vec2(1.,RENDERSIZE.y/RENDERSIZE.x);
    float ramp = 4.;

    float d = 0.125;

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
*/

vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    vec4 noise = texture(image30, fragCoord.xy / RENDERSIZE.xy + fract(vec2(42,56)*TIME));

    if(FRAMECOUNT<10)
    {
        fragColor = noise;
        return  fragColor;
    }


    uv = 0.5 + (uv - 0.5)*0.99;
    vec2 pixelSize = 1./RENDERSIZE.xy;
    vec2 mouseV = mouseDelta();
    vec2 aspect = vec2(1.,RENDERSIZE.y/RENDERSIZE.x);
    uv = vortex_pair_warp(uv, _mouse.xy*pixelSize, mouseV*aspect*1.4);

    float time = (smoothTime)/40.;
    uv = uv + vec2(sin(time*0.6 + uv.x*2. +1.) - sin(time*0.214 + uv.y*2. +1.), sin(time*0.168 + uv.x*2. +1.) - sin(time*0.115 +uv.y*2. +1.))*pixelSize*1.5;

    fragColor = BlurB(uv, 0);
    fragColor += ((BlurB(uv, 1) - BlurB(uv, 2))*0.5 + (noise-0.5) * 0.004)+highhits*0.5; 

    fragColor = clamp(fragColor, 0., 1.);

    //fragColor = noise; // reset
	return fragColor; 
 } 


			//******** BuffC Code Begins ********

// resolution reduction and horizontal blur
/*
vec2 lower_left(vec2 uv)
{
    return fract(uv * 0.5);
}

vec2 lower_right(vec2 uv)
{
    return fract((uv - vec2(1, 0.)) * 0.5);
}

vec2 upper_left(vec2 uv)
{
    return fract((uv - vec2(0., 1)) * 0.5);
}

vec2 upper_right(vec2 uv)
{
    return fract((uv - 1.) * 0.5);
}
*/
vec4 blur_horizontal(sampler2D channel, vec2 uv, float scale)
{
    float h = scale / RENDERSIZE.x;
    vec4 sum = vec4(0.0);

    sum += texture(channel, fract(vec2(uv.x - 4.0*h, uv.y)) ) * 0.05;
    sum += texture(channel, fract(vec2(uv.x - 3.0*h, uv.y)) ) * 0.09;
    sum += texture(channel, fract(vec2(uv.x - 2.0*h, uv.y)) ) * 0.12;
    sum += texture(channel, fract(vec2(uv.x - 1.0*h, uv.y)) ) * 0.15;
    sum += texture(channel, fract(vec2(uv.x + 0.0*h, uv.y)) ) * 0.16;
    sum += texture(channel, fract(vec2(uv.x + 1.0*h, uv.y)) ) * 0.15;
    sum += texture(channel, fract(vec2(uv.x + 2.0*h, uv.y)) ) * 0.12;
    sum += texture(channel, fract(vec2(uv.x + 3.0*h, uv.y)) ) * 0.09;
    sum += texture(channel, fract(vec2(uv.x + 4.0*h, uv.y)) ) * 0.05;

    return sum/0.98; // normalize
}

vec4 blur_horizontal_left_column(vec2 uv, int depth)
{
    float h = pow(2., float(depth)) / RENDERSIZE.x;    
    vec2 uv1, uv2, uv3, uv4, uv5, uv6, uv7, uv8, uv9;

    uv1 = fract(vec2(uv.x - 4.0 * h, uv.y) * 2.);
    uv2 = fract(vec2(uv.x - 3.0 * h, uv.y) * 2.);
    uv3 = fract(vec2(uv.x - 2.0 * h, uv.y) * 2.);
    uv4 = fract(vec2(uv.x - 1.0 * h, uv.y) * 2.);
    uv5 = fract(vec2(uv.x + 0.0 * h, uv.y) * 2.);
    uv6 = fract(vec2(uv.x + 1.0 * h, uv.y) * 2.);
    uv7 = fract(vec2(uv.x + 2.0 * h, uv.y) * 2.);
    uv8 = fract(vec2(uv.x + 3.0 * h, uv.y) * 2.);
    uv9 = fract(vec2(uv.x + 4.0 * h, uv.y) * 2.);

    if(uv.y > 0.5)
    {
        uv1 = upper_left(uv1);
        uv2 = upper_left(uv2);
        uv3 = upper_left(uv3);
        uv4 = upper_left(uv4);
        uv5 = upper_left(uv5);
        uv6 = upper_left(uv6);
        uv7 = upper_left(uv7);
        uv8 = upper_left(uv8);
        uv9 = upper_left(uv9);
    }
    else{
        uv1 = lower_left(uv1);
        uv2 = lower_left(uv2);
        uv3 = lower_left(uv3);
        uv4 = lower_left(uv4);
        uv5 = lower_left(uv5);
        uv6 = lower_left(uv6);
        uv7 = lower_left(uv7);
        uv8 = lower_left(uv8);
        uv9 = lower_left(uv9);
    }

    for(int level = 0; level < 8; level++)
    {
        if(level >= depth)
        {
            break;
        }

        uv1 = lower_right(uv1);
        uv2 = lower_right(uv2);
        uv3 = lower_right(uv3);
        uv4 = lower_right(uv4);
        uv5 = lower_right(uv5);
        uv6 = lower_right(uv6);
        uv7 = lower_right(uv7);
        uv8 = lower_right(uv8);
        uv9 = lower_right(uv9);
    }

    vec4 sum = vec4(0.0);

    sum += texture(BuffD, uv1) * 0.05;
    sum += texture(BuffD, uv2) * 0.09;
    sum += texture(BuffD, uv3) * 0.12;
    sum += texture(BuffD, uv4) * 0.15;
    sum += texture(BuffD, uv5) * 0.16;
    sum += texture(BuffD, uv6) * 0.15;
    sum += texture(BuffD, uv7) * 0.12;
    sum += texture(BuffD, uv8) * 0.09;
    sum += texture(BuffD, uv9) * 0.05;

    return sum/0.98; // normalize
}

vec4 renderPassC() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord.xy / RENDERSIZE.xy;

    if(uv.x < 0.5)
    {
        vec2 uv_half = fract(uv*2.);
        if(uv.y > 0.5)
        {
            fragColor = blur_horizontal(BuffA, uv_half, 1.);
        }
        else
        {
            fragColor = blur_horizontal(BuffB, uv_half, 1.);
        }
    }
    else
    {
        for(int level = 0; level < 8; level++)
        {
            if((uv.x > 0.5 && uv.y > 0.5) || (uv.x <= 0.5))
            {
                break;
            }
            vec2 uv_half = fract(uv*2.);
            fragColor = blur_horizontal_left_column(uv_half, level);
            uv = uv_half;
        }
    }
	return fragColor; 
 } 


			//******** BuffD Code Begins ********

// vertical blur (second pass)
/*
vec2 lower_left(vec2 uv)
{
    return fract(uv * 0.5);
}

vec2 lower_right(vec2 uv)
{
    return fract((uv - vec2(1, 0.)) * 0.5);
}

vec2 upper_left(vec2 uv)
{
    return fract((uv - vec2(0., 1)) * 0.5);
}

vec2 upper_right(vec2 uv)
{
    return fract((uv - 1.) * 0.5);
}
*/
vec4 blur_vertical_upper_left(sampler2D channel, vec2 uv)
{
    float v = 1. / RENDERSIZE.y;
    vec4 sum = vec4(0.0);
    sum += texture(channel, upper_left(vec2(uv.x, uv.y - 4.0*v)) ) * 0.05;
    sum += texture(channel, upper_left(vec2(uv.x, uv.y - 3.0*v)) ) * 0.09;
    sum += texture(channel, upper_left(vec2(uv.x, uv.y - 2.0*v)) ) * 0.12;
    sum += texture(channel, upper_left(vec2(uv.x, uv.y - 1.0*v)) ) * 0.15;
    sum += texture(channel, upper_left(vec2(uv.x, uv.y + 0.0*v)) ) * 0.16;
    sum += texture(channel, upper_left(vec2(uv.x, uv.y + 1.0*v)) ) * 0.15;
    sum += texture(channel, upper_left(vec2(uv.x, uv.y + 2.0*v)) ) * 0.12;
    sum += texture(channel, upper_left(vec2(uv.x, uv.y + 3.0*v)) ) * 0.09;
    sum += texture(channel, upper_left(vec2(uv.x, uv.y + 4.0*v)) ) * 0.05;
    return sum/0.98; // normalize
}

vec4 blur_vertical_lower_left(sampler2D channel, vec2 uv)
{
    float v = 1. / RENDERSIZE.y;
    vec4 sum = vec4(0.0);
    sum += texture(channel, lower_left(vec2(uv.x, uv.y - 4.0*v)) ) * 0.05;
    sum += texture(channel, lower_left(vec2(uv.x, uv.y - 3.0*v)) ) * 0.09;
    sum += texture(channel, lower_left(vec2(uv.x, uv.y - 2.0*v)) ) * 0.12;
    sum += texture(channel, lower_left(vec2(uv.x, uv.y - 1.0*v)) ) * 0.15;
    sum += texture(channel, lower_left(vec2(uv.x, uv.y + 0.0*v)) ) * 0.16;
    sum += texture(channel, lower_left(vec2(uv.x, uv.y + 1.0*v)) ) * 0.15;
    sum += texture(channel, lower_left(vec2(uv.x, uv.y + 2.0*v)) ) * 0.12;
    sum += texture(channel, lower_left(vec2(uv.x, uv.y + 3.0*v)) ) * 0.09;
    sum += texture(channel, lower_left(vec2(uv.x, uv.y + 4.0*v)) ) * 0.05;
    return sum/0.98; // normalize
}

vec4 blur_vertical_left_column(vec2 uv, int depth)
{
    float v = pow(2., float(depth)) / RENDERSIZE.y;

    vec2 uv1, uv2, uv3, uv4, uv5, uv6, uv7, uv8, uv9;

    uv1 = fract(vec2(uv.x, uv.y - 4.0*v) * 2.);
    uv2 = fract(vec2(uv.x, uv.y - 3.0*v) * 2.);
    uv3 = fract(vec2(uv.x, uv.y - 2.0*v) * 2.);
    uv4 = fract(vec2(uv.x, uv.y - 1.0*v) * 2.);
    uv5 = fract(vec2(uv.x, uv.y + 0.0*v) * 2.);
    uv6 = fract(vec2(uv.x, uv.y + 1.0*v) * 2.);
    uv7 = fract(vec2(uv.x, uv.y + 2.0*v) * 2.);
    uv8 = fract(vec2(uv.x, uv.y + 3.0*v) * 2.);
    uv9 = fract(vec2(uv.x, uv.y + 4.0*v) * 2.);

    if(uv.y > 0.5)
    {
        uv1 = upper_left(uv1);
        uv2 = upper_left(uv2);
        uv3 = upper_left(uv3);
        uv4 = upper_left(uv4);
        uv5 = upper_left(uv5);
        uv6 = upper_left(uv6);
        uv7 = upper_left(uv7);
        uv8 = upper_left(uv8);
        uv9 = upper_left(uv9);
    }
    else{
        uv1 = lower_left(uv1);
        uv2 = lower_left(uv2);
        uv3 = lower_left(uv3);
        uv4 = lower_left(uv4);
        uv5 = lower_left(uv5);
        uv6 = lower_left(uv6);
        uv7 = lower_left(uv7);
        uv8 = lower_left(uv8);
        uv9 = lower_left(uv9);
    }

    for(int level = 0; level < 8; level++)
    {
        if(level > depth)
        {
            break;
        }

        uv1 = lower_right(uv1);
        uv2 = lower_right(uv2);
        uv3 = lower_right(uv3);
        uv4 = lower_right(uv4);
        uv5 = lower_right(uv5);
        uv6 = lower_right(uv6);
        uv7 = lower_right(uv7);
        uv8 = lower_right(uv8);
        uv9 = lower_right(uv9);
    }

    vec4 sum = vec4(0.0);

    sum += texture(BuffC, uv1) * 0.05;
    sum += texture(BuffC, uv2) * 0.09;
    sum += texture(BuffC, uv3) * 0.12;
    sum += texture(BuffC, uv4) * 0.15;
    sum += texture(BuffC, uv5) * 0.16;
    sum += texture(BuffC, uv6) * 0.15;
    sum += texture(BuffC, uv7) * 0.12;
    sum += texture(BuffC, uv8) * 0.09;
    sum += texture(BuffC, uv9) * 0.05;

    return sum/0.98; // normalize
}

vec4 renderPassD() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    vec2 uv_orig = uv;
    vec2 uv_half = fract(uv*2.);
    if(uv.x < 0.5)
    {
        if(uv.y > 0.5)
        {
            fragColor = blur_vertical_upper_left(BuffC, uv_half);
        }
        else
        {
            fragColor = blur_vertical_lower_left(BuffC, uv_half);
        }
    }
    else
    {
        for(int level = 0; level < 8; level++)
        {
            if((uv.x > 0.5 && uv.y >= 0.5) || (uv.x < 0.5))
            {
                break;
            }
            vec2 uv_half = fract(uv*2.);
            fragColor = blur_vertical_left_column(uv_half, level);
            uv = uv_half;
        }  
    }
    uv = uv_orig;
    float eighth = 1./8.;
    if(uv.x > 7.*eighth && uv.x < 8.*eighth && uv.y > 2.*eighth && uv.y < 3.*eighth)
    {
        fragColor = vec4(_mouse.xy / RENDERSIZE.xy, _mouse.zw / RENDERSIZE.xy);
    }
	return fragColor; 
 } 

/*
vec2 lower_left(vec2 uv)
{
    return fract(uv * 0.5);
}

vec2 lower_right(vec2 uv)
{
    return fract((uv - vec2(1, 0.)) * 0.5);
}

vec2 upper_left(vec2 uv)
{
    return fract((uv - vec2(0., 1)) * 0.5);
}

vec2 upper_right(vec2 uv)
{
    return fract((uv - 1.) * 0.5);
}

vec4 BlurA(vec2 uv, int level)
{
    if(level <= 0)
    {
        return texture(BuffA, fract(uv));
    }

    uv = upper_left(uv);
    for(int depth = 1; depth < 8; depth++)
    {
        if(depth >= level)
        {
            break;
        }
        uv = lower_right(uv);
    }

    return texture(BuffD, uv);
}

vec4 BlurB(vec2 uv, int level)
{
    if(level <= 0)
    {
        return texture(BuffB, fract(uv));
    }

    uv = lower_left(uv);
    for(int depth = 1; depth < 8; depth++)
    {
        if(depth >= level)
        {
            break;
        }
        uv = lower_right(uv);
    }

    return texture(BuffD, uv);
}

vec2 GradientA(vec2 uv, vec2 d, vec4 selector, int level){
    vec4 dX = 0.5*BlurA(uv + vec2(1.,0.)*d, level) - 0.5*BlurA(uv - vec2(1.,0.)*d, level);
    vec4 dY = 0.5*BlurA(uv + vec2(0.,1.)*d, level) - 0.5*BlurA(uv - vec2(0.,1.)*d, level);
    return vec2( dot(dX, selector), dot(dY, selector) );
}
*/
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    vec2 pixelSize = 1. / RENDERSIZE.xy;
    vec2 aspect = vec2(1.,RENDERSIZE.y/RENDERSIZE.x);

    vec2 d = pixelSize*2.;
    vec4 dx = (BlurA(uv + vec2(1,0)*d, 1) - BlurA(uv - vec2(1,0)*d, 1))*0.5;
    vec4 dy = (BlurA(uv + vec2(0,1)*d, 1) - BlurA(uv - vec2(0,1)*d, 1))*0.5;

    d = pixelSize*1.;
    dx += BlurA(uv + vec2(1,0)*d, 0) - BlurA(uv - vec2(1,0)*d, 0);
    dy += BlurA(uv + vec2(0,1)*d, 0) - BlurA(uv - vec2(0,1)*d, 0);
    vec2 lightSize=vec2(0.5);

    fragColor = BlurA(uv+vec2(dx.x,dy.x)*pixelSize*8., 0).x * vec4(0.7,1.66,2.0,1.0) - vec4(0.3,1.0,1.0,1.0);
    fragColor = mix(fragColor,vec4(8.0,6.,2.,1.), BlurA(uv + vec2(dx.x,dy.x)*lightSize, 3).y*0.4*0.75*vec4(1.-BlurA(uv+vec2(dx.x,dy.x)*pixelSize*8., 0).x)); 
    fragColor = mix(fragColor, vec4(0.1,0.,0.4,0.), BlurA(uv, 1).a*length(GradientA(uv, pixelSize*2., vec4(0.,0.,0.,1.), 0))*5.);
    fragColor = mix(fragColor, vec4(1.25,1.35,1.4,0.), BlurA(uv, 0).x*BlurA(uv + GradientA(uv, pixelSize*2.5, vec4(-256.,32.,-128.,32.), 1)*pixelSize, 2).y);
    fragColor = mix(fragColor, vec4(0.25,0.75,1.,0.), BlurA(uv, 1).x*length(GradientA(uv+GradientA(uv, pixelSize*2., vec4(0.,0.,128.,0.), 1)*pixelSize, pixelSize*2., vec4(0.,0.,0.,1.), 0))*5.);
    fragColor = mix(fragColor, vec4(1.,1.25,1.5,0.), 0.5*(1.-BlurA(uv, 0)*1.).a*length(GradientA(uv+GradientA(uv, pixelSize*2., vec4(0.,128.,0.,0.), 1)*pixelSize, pixelSize*1.5, vec4(0.,0.,16.,0.), 0)));

    //    fragColor = BlurA(uv, 0); // simple bypass
    //    fragColor = BlurB(uv, 0); // simple bypass
    //    fragColor = texture(BuffD, uv); // raw Gaussian pyramid

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