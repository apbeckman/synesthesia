
// ****************** PASS 0 ***********************
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


vec4 BlurB(vec2 uv, int level)
{
    if(level <= 0)
    {
        return texture(buffB2, fract(uv));
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

    return texture(buffD, uv);
}

vec2 GradientB(vec2 uv, vec2 d, vec4 selector, int level){
    vec4 dX = 0.5*BlurB(uv + vec2(1.,0.)*d, level) - 0.5*BlurB(uv - vec2(1.,0.)*d, level);
    vec4 dY = 0.5*BlurB(uv + vec2(0.,1.)*d, level) - 0.5*BlurB(uv - vec2(0.,1.)*d, level);
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
    float f1lter = warpFilter(uv, pos_correct, size, ramp);
    return mix(uv, rot_uv, f1lter);
}

vec2 vortex_pair_warp(vec2 uv, vec2 pos, vec2 vel)
{
    vec2 aspect = vec2(1.,RENDERSIZE.y/RENDERSIZE.x);
    float ramp = 4.;

    float d = 0.425;

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

vec2 kaleidoscope(vec2 uvIn, float n) {
  vec2 uv = uvIn;
  float angle = PI/n;
  
  float r = length(uv);
  float a = atan(uv.y, uv.x)/angle;
  
  a = mix(fract(a), 1.0 - fract(a), mod(floor(a), 2.0))*angle;
  
  return vec2(cos(a), sin(a))*r;
}

vec2 autoSliceTransform(vec2 uv)
{
    vec2 diagPos = _rotate(_uvc, PI*0.25);
    float gridderX = floor(mod(diagPos.x*30.0,2.0));
    float gridderY = floor(mod(diagPos.y*20.0,2.0));
    float tinyGridX = floor(mod(diagPos.x*40.0,2.0));
    float tinyGridY = floor(mod(diagPos.y*50.0,2.0));
    int beatMode = int(mod(syn_BeatTime,4.0));
    float singer = -1+2.0*floor(mod(syn_BeatTime*0.5,2.0));

    vec2 gridder;

    switch(beatMode)
    {
    case 0: gridder = gridderX*vec2(0.0,1.0); break;
    case 1: gridder = gridderY*vec2(1.0,0.0); break;
    case 2: gridder = tinyGridX*vec2(0.0,1.0); break;
    case 3: gridder = tinyGridY*vec2(1.0,0.0); break;
    }
// gridder = _rotate(gridder, PI*0.25);
    return uv+gridder*pow(syn_HighHits,2.0)*0.01*singer*pow(syn_Intensity*syn_HighPresence+syn_Intensity*syn_BassPresence,2.0);
}

// *************** PASS 0 *****************
// Turing Sim

void mainImage1( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    vec4 noise = texture(colornoise, fragCoord.xy / RENDERSIZE.xy + fract(vec2(42,56)*TIME));

    if(FRAMECOUNT<10)
    {
        fragColor = noise;
        return;
    }

    float zoomAmt = 0.999+(inOrOut*(1.0-abs(in_or_out))+in_or_out*2.0)*0.0075*pow(syn_BassLevel*0.75+syn_MidLevel*0.25,2.0);
    // float zoomAmt = 1.0;


    uv = 0.5 + (uv - 0.5)*zoomAmt;
    vec2 pixelSize = 1./RENDERSIZE.xy;
    vec2 aspect = vec2(1.,RENDERSIZE.y/RENDERSIZE.x);
    if (slicing > 0.5){
        uv = autoSliceTransform(uv);
    }

    //uv = uv - vec2(0.0,GradientB(uv, pixelSize, vec4(-128,-128.,-128.,-128.), 1)*syn_HighHits*0.0001);
    float logo = dot(_loadUserImage().rgb, vec3(1.0))*0.33;

    uv -= vec2(-_uvc.x*0.005*basshits, _uvc.y*0.005*highhits)*(1.0-logo)*stretching;

    vec4 old = BlurB(uv, 0).rgbr;
    fragColor = old;
    fragColor += ((BlurB(uv+vec2(0.0,length(fragColor)*drips*0.01), 1) - BlurB(uv, 2))*0.5 + (noise-0.5) * 0.004); 

    float fbmmer = _fbm(_uvc*10.0+smoothTimeC);

    if (circle_on > 0.5){
        vec2 posSurf = _toPolarTrue(_uvc);
        
        posSurf.t = posSurf.t*(20.0);

        float stepper = 1.0;

        // stepper = step(0.2,(distance(_toRect(paintPos*vec2(1.0,2*PI)), _toRect(posSurf*vec2(1.0,2*PI))))*stepper);
        stepper = smoothstep(0.0, 0.4,posSurf.r-circle_size);
        fragColor.ra = mix(fragColor.ra, vec2(0.5, length(fragColor)*2.0), stepper);
    }

    fragColor.ra = mix(fragColor.ra, vec2(0.7, length(fragColor)*2.0), _pulse(_uv.y, down_sweep*1.1-0.1, 0.05)*step(0.001, down_sweep));
    fragColor.ra = mix(fragColor.ra, vec2(0.7, length(fragColor)*2.0), _pulse(_uv.y, (1.0-up_sweep)*1.1, 0.05)*step(0.001, up_sweep));
    fragColor.ra = mix(fragColor.ra, vec2(0.7, length(fragColor)*2.0), _pulse(_uv.x, left_sweep*1.1-0.1, 0.05)*step(0.001, left_sweep));
    fragColor.ra = mix(fragColor.ra, vec2(0.7, length(fragColor)*2.0), _pulse(_uv.x, (1.0-right_sweep)*1.1, 0.05)*step(0.001, right_sweep));

    if (syn_MediaType > 0.5){

        fragColor.r -= fragColor.r*logo*syn_BassLevel;
        fragColor.b += fragColor.g*logo;
        fragColor.g = mix(fragColor.g, logo, syn_HighHits);

    }

    fragColor = clamp(fragColor, 0., 1.);

    // fragColor = noise; // reset
}



// *************** PASS 2 *****************
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

    return sum/0.998; // normalize
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

    sum += texture(buffD, uv1) * 0.05;
    sum += texture(buffD, uv2) * 0.09;
    sum += texture(buffD, uv3) * 0.12;
    sum += texture(buffD, uv4) * 0.15;
    sum += texture(buffD, uv5) * 0.16;
    sum += texture(buffD, uv6) * 0.15;
    sum += texture(buffD, uv7) * 0.12;
    sum += texture(buffD, uv8) * 0.09;
    sum += texture(buffD, uv9) * 0.05;

    return sum/0.98; // normalize
}

void mainImage2( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord.xy / RENDERSIZE.xy;

    if(uv.x < 0.5)
    {
        vec2 uv_half = fract(uv*2.);
        if(uv.y > 0.5)
        {
            fragColor = blur_horizontal(buffB2, uv_half, 1.);
        }
        else
        {
            fragColor = blur_horizontal(buffB2, uv_half, 1.);
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
}


// void mainImage2( out vec4 fragColor, in vec2 fragCoord )
// {
//     vec2 uv = fragCoord.xy / RENDERSIZE.xy;

//     if(uv.x < 0.5)
//     {
//         vec2 uv_half = fract(uv*2.);
//         // if(uv.y > 0.5)
//         // {
//         //     fragColor = blur_horizontal(buffB, uv_half, 1.);
//         // }
//         // else
//         // {
//             fragColor = blur_horizontal(buffB, uv_half, 1.);
//         // }
//     }
//     else
//     {
//         for(int level = 0; level < 8; level++)
//         {
//             if((uv.x > 0.25 && uv.y > 0.25) || (uv.x <= 0.5))
//             {
//                 break;
//             }
//             vec2 uv_half = fract(uv*2.);
//             fragColor = blur_horizontal_left_column(uv_half, level);
//             uv = uv_half;
//         }
//     }
// }

// *************** PASS 4 *****************

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

    sum += texture(buffC, uv1) * 0.05;
    sum += texture(buffC, uv2) * 0.09;
    sum += texture(buffC, uv3) * 0.12;
    sum += texture(buffC, uv4) * 0.15;
    sum += texture(buffC, uv5) * 0.16;
    sum += texture(buffC, uv6) * 0.15;
    sum += texture(buffC, uv7) * 0.12;
    sum += texture(buffC, uv8) * 0.09;
    sum += texture(buffC, uv9) * 0.05;

    return sum/0.98; // normalize
}

void mainImage3( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    vec2 uv_orig = uv;
    vec2 uv_half = fract(uv*2.);
    if(uv.x < 0.5)
    {
        if(uv.y > 0.5)
        {
            fragColor = blur_vertical_upper_left(buffC, uv_half);
        }
        else
        {
            fragColor = blur_vertical_lower_left(buffC, uv_half);
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
}

// *************** Main Image Pass *****************

void mainImage4( out vec4 fragColor, in vec2 fragCoord )
{

    vec2 uv = fragCoord.xy / RENDERSIZE.xy;

    uv /= (1.0+length(_uvc*PI)*(0.5+0.5*sin(smoothTimeC*0.135))*distort);
    // if (uv.x*mod(syn_BeatTime*5.0+syn_RandomOnBeat*2.0,8.0)*0.5+uv.y*mod(syn_BeatTime*13.0+syn_RandomOnBeat*1.7,8.0)*0.5 > 1.0){
    //     uv = 1.0-uv;
    // }

    // if (uv.x > 0.5){
    //     uv.x = 1.0-uv.x;
    // }
    vec2 pixelSize = 1. / RENDERSIZE.xy;
    vec2 aspect = vec2(1.,RENDERSIZE.y/RENDERSIZE.x);

    vec2 lightSize=vec2(0.5);

    // *** Color Regime 0 ***
    vec3 midHighCol = vec3(-1.0)*0.5;
    vec3 midCol = vec3(1.0)*2.0;
    vec3 bassCol1 = vec3(1.0)*2.0;    
    vec3 bassCol2 = vec3(-1.0);
    vec3 highCol1 = vec3(-1.0);
    vec3 highCol2 = vec3(1.0);

    float mixer = color_palette;
    mixer = smoothstep(0.25, 0.75, clamp(mixer, 0.0, 1.0));

    // *** Color Regime 1 ***
    midHighCol = mix(midHighCol, _normalizeRGB(46, 9, 39), color_palette);
    bassCol1 = mix(bassCol1, vec3(0.7,0.3,0.1), color_palette);
    highCol1 = mix(highCol1, _normalizeRGB(4, 10, 111), color_palette);
    bassCol2 = mix(bassCol2, _normalizeRGB(255, 12, 0), color_palette);
    highCol2 = mix(highCol2, _normalizeRGB(255, 140, 255), color_palette);
    midCol = mix(midCol, _normalizeRGB(4, 117, 111), color_palette);

    mixer = color_palette - 1.0;
    mixer = smoothstep(0.25, 0.75, clamp(mixer, 0.0, 1.0));

    // *** Color Regime 2 ***
    midHighCol = mix(midHighCol, vec3(0.0,0.0,1.0)*0.5, mixer);
    bassCol1 = mix(bassCol1,vec3(0.8,0.4,0.0)*2.0, mixer);
    highCol1 = mix(highCol1,vec3(1.0,0.4,0.0)*2.0, mixer);
    bassCol2 = mix(bassCol2,vec3(0.8,0.3,0.0), mixer);
    highCol2 = mix(highCol2,vec3(0.0,0.7,1.0), mixer);
    midCol = mix(midCol,vec3(0.0,0.8,0.9), mixer);

    vec2 zoneMids = vec2(sin(smoothTimeC*0.23), cos(smoothTimeC*0.17));
    vec2 zoneMidHighs = vec2(sin(syn_MidHighTime*0.23), cos(syn_MidHighTime*0.17));
    vec2 zoneHighs = vec2(sin(smoothTimeB*0.23), cos(smoothTimeB*0.17));

    midHighCol *= mix(1.0,(0.4+0.7*syn_MidHighPresence*distance(zoneMidHighs, _uvc)*0.5), 1.0-no_shading);
    bassCol1 *= mix(1.0,(0.3+0.7*pow(syn_BassLevel,2.0)), 1.0-no_shading);
    bassCol2 *= mix(1.0,(0.2+0.9*syn_BassPresence*1.0*syn_Intensity), 1.0-no_shading);
    highCol1 *= mix(1.0,(0.2+0.9*syn_HighPresence*syn_Intensity), 1.0-no_shading);
    midCol *= mix(1.0,(0.2+0.9*distance(zoneMids, _uvc)*syn_MidPresence), 1.0-no_shading);

    if (flashing > 0.5){
        highCol2 *= (0.1+0.9*syn_HighLevel*(0.5+syn_Intensity*0.5+syn_HighHits*0.25)*distance(zoneHighs, _uvc));
    } else {
        highCol2 *= ((0.5+syn_Presence*0.5)*distance(zoneHighs, _uvc));
    }
    
    float logo = dot(_loadUserImage().rgb, vec3(1.0))/3.0;
    bassCol2 *= (1.0-logo);
    bassCol1 *= (1.0-logo);


    vec3 finalColor = vec3(0.0);
    finalColor = mix(finalColor, midCol*1.5, clamp(BlurB(uv + GradientB(uv, pixelSize, vec4(-256.,32.,-128.,32.)*(1.0-logo), 0)*pixelSize, 0).y,0.0,1.0));
    finalColor = mix(finalColor, bassCol1, clamp(abs(GradientB(uv, pixelSize*2.0, vec4(2.0,0.0,0.0,0.0),0).y),0.0,1.0));
    finalColor = mix(finalColor, midHighCol, clamp(BlurB(uv, 4).a,0.0,1.0));
    finalColor = mix(finalColor, bassCol2, clamp(BlurB(uv, 0).x*BlurB(uv + GradientB(uv, pixelSize*2.5, vec4(-256.,32.,-128.,32.)*0.2, 1)*pixelSize, 1).x,0.0,1.0));
    finalColor = mix(finalColor, highCol1, clamp(dot(GradientB(uv, pixelSize*1., vec4(0., 10., 0., 0.), 1), vec2(sin(syn_HighTime*1.0), cos(syn_HighTime*1.0))),0.0,1.0));
    finalColor = mix(finalColor, highCol2*2.0, clamp(BlurB(uv, 1).b*dot(GradientB(uv, pixelSize*2.0, vec4(1.,0.,0.,0.), 0), vec2(sin(syn_HighTime*2.0), cos(syn_HighTime*2.0)))*5.,0.0,1.0));

    finalColor *= mix(1.0,(0.6+(sin(_uvc.x*5.0+TIME*0.1)*_uv.y+cos(_uvc.x*4.7-TIME*0.47)+cos(_uvc.x*11.0+TIME*0.89)*(1.0-_uv.y))/6.0),1.0-no_shading);
    finalColor += pow(finalColor,vec3(2.0))*(0.5+syn_Presence*0.5);

    vec2 vigPos = _uv*(1.0 - uv.yx);   //vec2(1.0)- uv.yx; -> 1.-u.yx; Thanks FabriceNeyret !
    float vig = vigPos.x*vigPos.y * 15.0; // multiply with sth for intensity
    vig = pow(vig, 0.25); // change pow for modifying the extend of the  vignette

    finalColor = clamp(finalColor, 0.0, 1.0);

    // finalColor = _rgb2hsv(finalColor);
    // finalColor.g += BlurB(uv, 3).r+TIME*0.1;
    // finalColor = _hsv2rgb(finalColor);

    fragColor = vec4(finalColor,1.0);
    fragColor *= vig;

       // fragColor = vec4(BlurB(uv, 0)); // simple bypass
       // fragColor -= vec4(0.5,0.5,0.5,0.0)*fragColor.a;
       // fragColor = vec4(BlurB(uv, 0).a); // simple bypass
       // fragColor = texture(buffD, _uv); // raw Gaussian pyramid
}

vec4 renderMain(void)
{
  if (PASSINDEX == 0.0){
    //buffB
    vec4 fragColor;
    mainImage1(fragColor, gl_FragCoord.xy);
    return fragColor;
  }else if (PASSINDEX == 1.0){
    return texture(buffB, _uv);
  }
  else if (PASSINDEX == 2.0){
    //buffC
    vec4 fragColor;
    mainImage2(fragColor, gl_FragCoord.xy);
    return fragColor;
  }
  else if (PASSINDEX == 3.0){
    //buffD
    vec4 fragColor;
    mainImage3(fragColor, gl_FragCoord.xy);
    return fragColor;
  }
  else if (PASSINDEX == 4.0){
    //Image
    vec4 fragColor;
    mainImage4(fragColor, gl_FragCoord.xy);
    return fragColor;
  }
}
