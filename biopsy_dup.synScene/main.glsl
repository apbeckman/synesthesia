vec3 iResolution = vec3(RENDERSIZE.x, RENDERSIZE.y, 0.0);
vec3 iChannelResolution = iResolution;
float iGlobalTime = TIME*1.0;
int iFrame = int(FRAMECOUNT);
float time = TIME;
vec2 resolution = RENDERSIZE;
vec4 iMouse = vec4(0.5);

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
        return texture(buffB, fract(uv));
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
    vec2 aspect = vec2(1.,iResolution.y/iResolution.x);
    return max( min, 1. - length((uv - pos) * aspect / size) );
}

float warpFilter(vec2 uv, vec2 pos, float size, float ramp)
{
    return 0.5 + sigmoid( conetip(uv, pos, size, -16.) * ramp) * 0.5;
}

vec2 vortex_warp(vec2 uv, vec2 pos, float size, float ramp, vec2 rot)
{
    vec2 aspect = vec2(1.,iResolution.y/iResolution.x);

    vec2 pos_correct = 0.5 + (pos - 0.5);
    vec2 rot_uv = pos_correct + complex_mul((uv - pos_correct)*aspect, rot)/aspect;
    float f1lter = warpFilter(uv, pos_correct, size, ramp);
    return mix(uv, rot_uv, f1lter);
}

vec2 vortex_pair_warp(vec2 uv, vec2 pos, vec2 vel)
{
    vec2 aspect = vec2(1.,iResolution.y/iResolution.x);
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
    return uv+gridder*pow(syn_HighHits,2.0)*0.01*singer*pow(syn_HighPresence+syn_BassPresence,2.0);
}

vec2 manualSliceTransform(vec2 uv)
{
    // vec4 data = texture(pass0, _uv);

    // vec2 direction = normalize(vec2(slice_pos.x, slice_pos.y)-vec2(sliceX, sliceY));

    // float angle = atan(direction.y/direction.x);
    vec2 direction = slicerDirection;
    float angle = slicerAngle;

    vec2 diagPos = _rotate(_uvc, angle+PI/2.0);

    vec2 flowDirection = -1.0*direction;
    float gridderX = floor(mod(diagPos.x*10.0,2.0));
    float gridderY = floor(mod(diagPos.y*10.0,2.0));
    float tinyGridX = floor(mod(diagPos.x*30.0,2.0));
    float tinyGridY = floor(mod(diagPos.y*30.0,2.0));
    int sliceMode = 3;

    float gridder;

    switch(sliceMode)
    {
    case 0: gridder = gridderX; break;
    case 1: gridder = gridderY; break;
    case 2: gridder = tinyGridX; break;
    case 3: gridder = tinyGridY; break;
    }
// gridder = _rotate(gridder, PI*0.25);
    return uv+direction*gridder*pow(slicerMotionAmt*0.01,1.0);
}

// *************** PASS 1 *****************
// Turing Sim

void mainImage1( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord.xy / iResolution.xy;
    vec4 noise = texture(colornoise, fragCoord.xy / iChannelResolution.xy + fract(vec2(42,56)*smoothTime));

    if(iFrame<10)
    {
        fragColor = noise;
        return;
    }

    float zoomAmt = 0.999+inOrOut*0.0075*pow(syn_BassLevel*(0.5+zoom),2.0 )*zoom;
    // float zoomAmt = 1.0;


    uv = 0.5 + (uv - 0.5)*zoomAmt;
    vec2 pixelSize = 1./iResolution.xy;
    vec2 mouseV = vec2(0.0, -(syn_BassLevel*0.25+syn_MidLevel*0.25)*0.01+syn_BPMTri2*0.001);
    vec2 aspect = vec2(1.,iResolution.y/iResolution.x);
    if (auto_slice > 0.5){
        uv = autoSliceTransform(uv);
    }
    else {
        uv = manualSliceTransform(uv);
    }


    float time = syn_Time;
    //uv = uv - vec2(0.0,GradientB(uv, pixelSize, vec4(-128,-128.,-128.,-128.), 1)*syn_HighHits*0.0001);


    fragColor = BlurB(uv, 0).rbbr;
    fragColor += ((BlurB(uv+vec2(0.0,length(fragColor)*drips*0.01), 1) - BlurB(uv, 2))*0.5 + (noise-0.5) * 0.004); 

    if (paint_on > 0.5){
        float fbmmer = _fbm(_uvc*10.0+TIME)*0.05;

        if (paint_invert < 0.5){
            if (distance(_uvc, vec2(paint_pos.x, paint_pos.y)*vec2(1.0,RENDERSIZE.y/RENDERSIZE.x)) < paint_size-fbmmer){
                fragColor = vec4(length(fragColor)*2.0);
            }
        } else {
            if (distance(_uvc, vec2(paint_pos.x, paint_pos.y)*vec2(1.0,RENDERSIZE.y/RENDERSIZE.x)) > paint_size-fbmmer){
                fragColor = vec4(length(fragColor)*2.0);
            }
        }
    }

    if (_exists(syn_UserImage)){
        fragColor = mix(fragColor, 1.0-fragColor,_loadUserImageAsMask());
    }


    fragColor = clamp(fragColor, 0., 1.);

    // fragColor = noise; // reset
}



// *************** PASS 2 *****************
vec4 blur_horizontal(sampler2D channel, vec2 uv, float scale)
{
    float h = scale / iResolution.x;
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
    float h = pow(2., float(depth)) / iResolution.x;    
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

    for(int level = 0; level < 10; level++)
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
    vec2 uv = fragCoord.xy / iResolution.xy;

    if(uv.x < 0.5)
    {
        vec2 uv_half = fract(uv*2.);
        if(uv.y > 0.5)
        {
            fragColor = blur_horizontal(buffB, uv_half, 1.);
        }
        else
        {
            fragColor = blur_horizontal(buffB, uv_half, 1.);
        }
    }
    else
    {
        for(int level = 0; level < 10; level++)
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
//     vec2 uv = fragCoord.xy / iResolution.xy;

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
    float v = 1. / iResolution.y;
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
    float v = 1. / iResolution.y;
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
    float v = pow(2., float(depth)) / iResolution.y;

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

    for(int level = 0; level < 10; level++)
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
    vec2 uv = fragCoord.xy / iResolution.xy;
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
        for(int level = 0; level < 10; level++)
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
        fragColor = vec4(iMouse.xy / iResolution.xy, iMouse.zw / iResolution.xy);
    }
}

// *************** Main Image Pass *****************

void mainImage4( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord.xy / iResolution.xy;

    // if (uv.x*mod(syn_BeatTime*5.0+syn_RandomOnBeat*2.0,8.0)*0.5+uv.y*mod(syn_BeatTime*13.0+syn_RandomOnBeat*1.7,8.0)*0.5 > 1.0){
    //     uv = 1.0-uv;
    // }

    // if (uv.x > 0.5){
    //     uv.x = 1.0-uv.x;
    // }
    vec2 pixelSize = 2. / iResolution.xy;
    vec2 aspect = vec2(1.,iResolution.y/iResolution.x);

    vec2 d = pixelSize*(1.+mode*2.0);
    vec4 dx = (BlurB(uv + vec2(1,0)*d, 1) - BlurB(uv - vec2(1,0)*d, 1))*0.5;
    vec4 dy = (BlurB(uv + vec2(0,1)*d, 1) - BlurB(uv - vec2(0,1)*d, 1))*0.5;

    d = pixelSize*1.;
    dx += BlurB(uv + vec2(1,0)*d, 0) - BlurB(uv - vec2(1,0)*d, 0);
    dy += BlurB(uv + vec2(0,1)*d, 0) - BlurB(uv - vec2(0,1)*d, 0);
    vec2 lightSize=vec2(0.5);

    //for B
    //midHighCol is general color?
    //bassCol1 is a color in between coral
    //highCol1 does almost nothing. creeps in on the edges.
    //bassCol2 is general color "on" the life spots
    //highCol2 does almost nothing when sim is running. creeps in from edges.
    //midCol is squiggles layer... Its red and white tho. 

    //for A
    //midHighCol is sheen near the wavefronts. like "newSpawn"
    //bassCol1 is big wavefronts and smooth shading
    //highCol1 is wavefronts
    //bassCol2 crisp wavefronts and smooth wavefronts, of a diff color probably.
    //highCol2 sort of ghostly messy wavefronts
    //midCol sort of a networky layer. Looks brainular.


    vec3 midHighCol = vec3(1.0)*4.0;
    vec3 bassCol1 = vec3(1.0)*2.0;
    vec3 highCol1 = vec3(1.0)*2.0;
    vec3 bassCol2 = vec3(1.0)*2.0;
    vec3 highCol2 = vec3(1.0)*2.0;
    vec3 midCol = vec3(1.0)*2.0;

    float cReg4 = 1.0;

    // *** Color Regime 4 ***
    midHighCol = mix(midHighCol, _normalizeRGB(46, 9, 39)*4.0, cReg4);
    bassCol1 = mix(bassCol1, _normalizeRGB(217, 0, 0)*2.0, cReg4);
    highCol1 = mix(highCol1, _normalizeRGB(4, 200, 111)*2.0, cReg4);
    bassCol2 = mix(bassCol2, _normalizeRGB(100, 45, 0)*2.0, cReg4);
    highCol2 = mix(highCol2, _normalizeRGB(255, 140, 0)*2.0, cReg4);
    midCol = mix(midCol, _normalizeRGB(4, 117, 111)*2.0, cReg4);

    vec2 zoneMids = vec2(sin(smoothTimeC*0.2), cos(smoothTimeC*0.15));
    vec2 zoneMidHighs = vec2(sin(smoothTimeC*0.2), cos(smoothTimeC*0.15));
    vec2 zoneHighs = vec2(sin(smoothTimeB*0.2), cos(smoothTimeB*0.15));

    midHighCol *= (0.1+0.9*syn_MidHighPresence*distance(zoneMidHighs, _uvc)*0.8);
    bassCol1 *= (0.3+0.9*(syn_BassPresence)*1.25);
    highCol1 *= (0.3+0.9*syn_HighPresence);
    bassCol2 *= (0.3+0.9*syn_BassPresence*1.5);
    if (flashing > 0.5){
        highCol2 *= (0.1+0.9*(pow((syn_HighHits*0.25+syn_HighLevel*0.25), 2.0)*0.75)*(0.5+syn_Presence*0.5)*distance(zoneHighs, _uvc));
    } else {
        highCol2 *= (0.2*(0.5+syn_Presence*0.5)*distance(zoneHighs, _uvc));
    }
    midCol *= (0.1+0.9*distance(zoneMids, _uvc)*syn_MidPresence);

    float lowIntensity = 1.0-syn_Presence;

    vec3 finalColor = vec3(0.0);
    finalColor = BlurB(uv+vec2(dx.x,dy.x)*pixelSize*8., 0).x * midHighCol;
    finalColor = mix(finalColor, bassCol1, BlurB(uv + vec2(dx.x,dy.x)*lightSize, 3).y*0.4*0.75*vec3(1.-BlurB(uv+vec2(dx.x,dy.x)*pixelSize*8., 0).x));
    finalColor = mix(finalColor, highCol1, BlurB(uv, 1).a*length(GradientB(uv, pixelSize*2., vec4(0.,0.,0., 1.), 0))*5.);
    finalColor = mix(finalColor, bassCol2, BlurB(uv, 0).x*BlurB(uv + GradientB(uv, pixelSize*2.5, vec4(-256.,32.,-128.,32.), 1)*pixelSize, 2).y);
    finalColor = mix(finalColor, highCol2, BlurB(uv, 1).x*length(GradientB(uv, pixelSize*2., vec4(0.,0.,0.,5.), 0))*5.);
    finalColor = mix(finalColor, midCol, 0.5*(1.-BlurB(uv, 0)*1.).g*length(GradientB(uv+GradientB(uv, pixelSize*2., vec4(0.,128.,0.,0.), 1)*pixelSize, pixelSize*1.5, vec4(0.,0.,16.,0.), 0)));

    finalColor *= 0.5+(sin(_uvc.x*5.0+smoothTimeB*0.25)*_uv.y+cos(_uvc.x*4.7-smoothTimeB*0.6)+cos(_uvc.x*11.0+smoothTimeB*0.89)*(1.0-_uv.y))/6.0;

    fragColor = color_mix-vec4(finalColor, 1.0);
    if (dot(fragColor.rgb, vec3(1.0))<0.0){
        fragColor = -fragColor;
    }
    fragColor = clamp(fragColor, 0.0, 1.0);
    //fragColor = pow(fragColor, vec4(2.5));
    fragColor *= (1.0-_uv.y*(2.0-length(fragColor)));
    // fragColor *= (1.1-_fbm(_uvc*20.0+TIME)*0.2);
    // fragColor *= 0.0+step(length(_uvc), 0.7-length(fragColor)*0.2);

    //    fragColor = BlurB(uv, 0); // simple bypass
    //    fragColor = BlurB(uv, 0); // simple bypass
       // fragColor = texture(buffD, _uv); // raw Gaussian pyramid

}

void mainImage5( out vec4 fragColor, in vec2 fragCoord )
{
    vec3 finalColor = texture(coloredSim, _uv).rgb;
    fragColor = vec4(finalColor, 1.0);
}


vec4 renderMain(void)
{
  if (PASSINDEX == 0.0){
    //buffB
    vec4 fragColor;
    mainImage1(fragColor, gl_FragCoord.xy);
    return fragColor;
  }
  else if (PASSINDEX == 1.0){
    //buffC
    vec4 fragColor;
    mainImage2(fragColor, gl_FragCoord.xy);
    return fragColor;
  }
  else if (PASSINDEX == 2.0){
    //buffD
    vec4 fragColor;
    mainImage3(fragColor, gl_FragCoord.xy);
    return fragColor;
  }
  else if (PASSINDEX == 3.0){
    //Image
    vec4 fragColor;
    mainImage4(fragColor, gl_FragCoord.xy);
    return fragColor;
  }
  else if (PASSINDEX == 4.0){
    //Image
    vec4 fragColor;
    mainImage5(fragColor, gl_FragCoord.xy);
    return fragColor;
  }
}
